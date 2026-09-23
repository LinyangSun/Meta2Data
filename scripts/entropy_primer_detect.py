#!/usr/bin/env python3
"""
Four-state (C/D/T/Q) primer detection and removal
(with mixed R1/R2 orientation support).

Use all quality-filtered reads from the selected sample and the first 20 bases
of every read to build the frequency matrix. C/D/T/Q represent
one/two/three/four bases with at least 10% read support and contribute fold
factors 1/2/3/4.

Match the caller-selected direction-specific primer database first and trim at
the matched endpoint. Only when no database primer matches, calculate fold
across the first 20 positions. Fold strictly below 16 identifies an unknown
primer and triggers a fixed 20-bp trim; fold 16 or greater leaves reads
unchanged.
"""
import sys
import os
import gzip
import math
import shutil
import argparse
import json
from collections import Counter
from itertools import zip_longest

import parameters
import read_layout


DETECTION_WINDOW = 20
UNKNOWN_PRIMER_FOLD_THRESHOLD = 16


# ===========================================================================
# Low-level FASTQ helpers (avoid BioPython for speed in trimming step)
# ===========================================================================

def _open_fq(filepath, mode="rt"):
    """Open FASTQ, transparently handling gzip."""
    filepath = str(filepath)
    if filepath.lower().endswith(".gz"):
        return gzip.open(filepath, mode)
    return open(filepath, mode.replace("t", ""))


def _iter_fastq(handle):
    """Yield (header, seq, plus, qual) tuples from a FASTQ handle."""
    while True:
        header = handle.readline()
        if not header:
            break
        header = header.rstrip("\n")
        seq = handle.readline().rstrip("\n")
        plus = handle.readline().rstrip("\n")
        qual = handle.readline().rstrip("\n")
        if not header.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
            raise ValueError("Malformed FASTQ record: missing header/plus line or unequal sequence and quality lengths")
        yield header, seq, plus, qual


# ===========================================================================
# Step 1: Read & pre-filter
# ===========================================================================

def _is_degraded_quality(filepath, n_reads=1000):
    """
    Check if a FASTQ file has degraded (dummy) quality scores.
    Degraded = only 1 unique quality character across sampled reads,
    meaning scores are placeholder values (e.g., all '#' from SRA/ENA
    for old 454 data). Quality filtering is meaningless for such data.
    """
    qual_chars = set()
    count = 0
    with _open_fq(filepath) as fh:
        for _header, _seq, _plus, qual in _iter_fastq(fh):
            qual_chars.update(qual)
            count += 1
            if count >= n_reads:
                break
    return len(qual_chars) <= 1


def _avg_quality(qual_str):
    """Average Phred quality score (Phred+33 encoding)."""
    if not qual_str:
        return 0.0
    return sum(ord(c) - 33 for c in qual_str) / len(qual_str)


def _sequence_complexity(seq):
    """Linguistic complexity = unique 2-mers / 16."""
    if len(seq) < 2:
        return 0.0
    kmers = set()
    for i in range(len(seq) - 1):
        kmers.add(seq[i:i + 2].upper())
    return len(kmers) / 16.0


def _read_entropy(seq):
    """
    Per-read Shannon entropy of base composition.
    H = -sum(p * log2(p)) for ACGT frequencies within the read.
    Max = 2.0 (equal ACGT), min = 0.0 (homopolymer).
    Filters out compositionally biased reads (e.g., poly-A, AT-rich artifacts).
    """
    if len(seq) < 2:
        return 0.0
    counts = Counter(seq.upper())
    total = sum(counts.get(b, 0) for b in 'ACGT')
    if total == 0:
        return 0.0
    h = 0.0
    for b in 'ACGT':
        c = counts.get(b, 0)
        if c > 0:
            p = c / total
            h -= p * math.log2(p)
    return h


def _filter_reason(seq, qual, min_len=50, min_avg_qual=20,
                   min_complexity=0.3, min_entropy=1.0, skip_qual=False):
    """Return the first failed filter, shared by independent and paired reads."""
    if len(seq) < min_len:
        return 'length'
    if not skip_qual and _avg_quality(qual) < min_avg_qual:
        return 'quality'
    if _sequence_complexity(seq) < min_complexity:
        return 'complexity'
    if _read_entropy(seq) < min_entropy:
        return 'entropy'
    return None


def read_and_filter(filepath, min_len=50,
                    min_avg_qual=20, min_complexity=0.3,
                    min_entropy=1.0, skip_qual=False):
    """
    Step 1: Read all reads from a FASTQ file and apply quality filters.
    Filters:
      - Length >= min_len
      - Average Phred quality >= min_avg_qual (skipped if skip_qual=True)
      - K-mer complexity >= min_complexity (unique 2-mers / 16)
      - Shannon entropy >= min_entropy (base composition diversity)
    Returns list of (seq_string, qual_string) tuples.
    """
    reads = []
    total = 0
    discarded = Counter()

    with _open_fq(filepath) as fh:
        for header, seq, plus, qual in _iter_fastq(fh):
            total += 1

            reason = _filter_reason(seq, qual, min_len, min_avg_qual,
                                    min_complexity, min_entropy, skip_qual)
            if reason:
                discarded[reason] += 1
                continue

            reads.append((seq, qual))

    print(f"  Scanned {total} reads, kept {len(reads)}", file=sys.stderr)
    if skip_qual:
        print(f"  Quality filtering: SKIPPED (degraded/dummy scores)",
              file=sys.stderr)
    print(f"  Discarded: {discarded['length']} (length<{min_len}), "
          f"{discarded['quality']} (avgQ<{min_avg_qual}), "
          f"{discarded['complexity']} (complexity<{min_complexity}), "
          f"{discarded['entropy']} (entropy<{min_entropy})", file=sys.stderr)

    if len(reads) < 500:
        print(f"  WARNING: Only {len(reads)} reads passed filters (< 500)",
              file=sys.stderr)

    return reads


# ===========================================================================
# Step 2: Position-wise base frequency matrix
# ===========================================================================

def build_frequency_matrix(reads, num_positions=DETECTION_WINDOW):
    """
    Build a num_positions x 4 frequency matrix.
    Column order: A=0, C=1, G=2, T=3.
    """
    base_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    counts = [[0, 0, 0, 0] for _ in range(num_positions)]

    for seq, _ in reads:
        upper = seq.upper()
        for i in range(min(num_positions, len(upper))):
            b = upper[i]
            if b in base_idx:
                counts[i][base_idx[b]] += 1

    matrix = []
    for i in range(num_positions):
        total = sum(counts[i])
        if total > 0:
            matrix.append([c / total for c in counts[i]])
        else:
            matrix.append([0.25, 0.25, 0.25, 0.25])

    return matrix


# ===========================================================================
# Phase 1: Mixed-orientation detection
# ===========================================================================

def detect_mixed_orientation(freq_matrix, reads, num_check=15):
    """
    Detect mixed R1/R2 orientation from the position-wise frequency matrix.

    Mixed orientation creates bimodal base distributions in the primer region:
    at each position, the forward primer base and reverse primer base produce
    two distinct frequency peaks.

    Criterion per position (i in 0..num_check-1):
        Sort frequencies desc → f1, f2, f3, f4
        Bimodal if: f1+f2 > 0.85  AND  f2 > 0.15  AND  f1 < 0.85

    If >50% of checked positions are bimodal → mixed orientation.

    Returns None if single orientation, or a dict:
        best_pos       – position with the clearest bimodal split
        majority_base  – dominant base at best_pos
        minority_base  – second base at best_pos
        majority_reads – list of (seq, qual) in majority group
        minority_reads – list of (seq, qual) in minority group
        ratio          – fraction of minority reads
    """
    bimodal_count = 0
    best_pos = -1
    best_f2 = 0.0

    n_check = min(num_check, len(freq_matrix))
    for i in range(n_check):
        freqs = sorted(freq_matrix[i], reverse=True)
        f1, f2 = freqs[0], freqs[1]
        if f1 + f2 > 0.85 and f2 > 0.15 and f1 < 0.85:
            bimodal_count += 1
            if f2 > best_f2:
                best_f2 = f2
                best_pos = i

    bimodal_ratio = bimodal_count / n_check
    print(f"\n  [Mixed-orientation check] "
          f"Bimodal positions: {bimodal_count}/{n_check} "
          f"(ratio={bimodal_ratio:.2f})", file=sys.stderr)

    if bimodal_ratio <= 0.5:
        print(f"  -> Single orientation", file=sys.stderr)
        return None

    # Identify the two dominant bases at best_pos
    idx_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    freqs_at_best = freq_matrix[best_pos]
    sorted_idx = sorted(range(4), key=lambda x: freqs_at_best[x],
                        reverse=True)
    majority_base = idx_base[sorted_idx[0]]
    minority_base = idx_base[sorted_idx[1]]

    # Split sampled reads by the base at best_pos
    majority_reads = []
    minority_reads = []
    for seq, qual in reads:
        if best_pos < len(seq) and seq[best_pos].upper() == minority_base:
            minority_reads.append((seq, qual))
        else:
            majority_reads.append((seq, qual))

    ratio = len(minority_reads) / max(1, len(reads))
    print(f"  -> MIXED ORIENTATION DETECTED", file=sys.stderr)
    print(f"  Best split position: {best_pos} "
          f"(majority='{majority_base}' "
          f"[{freqs_at_best[sorted_idx[0]]:.2f}], "
          f"minority='{minority_base}' "
          f"[{freqs_at_best[sorted_idx[1]]:.2f}])", file=sys.stderr)
    print(f"  Majority: {len(majority_reads)} reads, "
          f"Minority: {len(minority_reads)} reads "
          f"(flip ratio={ratio:.2f})", file=sys.stderr)

    return {
        'best_pos': best_pos,
        'majority_base': majority_base,
        'minority_base': minority_base,
        'majority_reads': majority_reads,
        'minority_reads': minority_reads,
        'ratio': ratio,
    }


# ===========================================================================
# Step 3: Four-state position classification (C / D / T / Q)
# ===========================================================================

BASES = ['A', 'C', 'G', 'T']
STATE_BY_SUPPORTED_BASE_COUNT = {1: 'C', 2: 'D', 3: 'T', 4: 'Q'}
STATE_FOLD_FACTORS = {'C': 1, 'D': 2, 'T': 3, 'Q': 4}
MIN_SUPPORTED_BASE_FREQUENCY = 0.10


def classify_position(freqs, min_supported_frequency=MIN_SUPPORTED_BASE_FREQUENCY):
    """Classify by the number of bases with meaningful read support.

    C/D/T/Q mean one/two/three/four supported bases. Frequencies below 10%
    are ignored so ordinary sequencing errors do not create false alleles.
    """
    supported = sum(f >= min_supported_frequency for f in freqs)
    supported = max(1, min(4, supported))
    return STATE_BY_SUPPORTED_BASE_COUNT[supported]


def classify_all_positions(freq_matrix, min_supported_frequency=MIN_SUPPORTED_BASE_FREQUENCY):
    """Classify every position in the frequency matrix as C/D/T/Q."""
    return [classify_position(freqs, min_supported_frequency) for freqs in freq_matrix]


# ===========================================================================
# Step 4: Fold calculation
# ===========================================================================

def calculate_primer_fold(states, boundary):
    """Calculate C/D/T/Q fold only from position 1 through the boundary."""
    if boundary < 0 or boundary > len(states):
        raise ValueError("boundary is outside the state sequence")
    fold = 1
    for state in states[:boundary]:
        try:
            fold *= STATE_FOLD_FACTORS[state]
        except KeyError as error:
            raise ValueError(f"Unknown primer state: {state!r}") from error
    return fold


# ===========================================================================
# Step 5: Consensus sequence
# ===========================================================================

def build_consensus_cdv(freq_matrix, states, primer_length):
    """
    Build consensus sequence for the detected primer region.

    C positions → dominant base (A/C/G/T)
    D/T/Q positions → N (degenerate)
    """
    consensus = []
    for i in range(primer_length):
        if states[i] == 'C':
            max_idx = freq_matrix[i].index(max(freq_matrix[i]))
            consensus.append(BASES[max_idx])
        else:  # D, T, or Q
            consensus.append('N')
    return ''.join(consensus)


# ===========================================================================
# Known 16S primer database (primary detection path)
# ===========================================================================
# A database match is accepted first and determines the trim endpoint. The
# first-20/fold-<16 fixed-20-bp rule runs only when this search has no match.

_IUPAC = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'},
    'W': {'A', 'T'}, 'K': {'G', 'T'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'},
}

_COMPLEMENT = str.maketrans(
    'ACGTRYSWKMBDHVNacgtryswkmbdhvn',
    'TGCAYRSWMKVHDBNtgcayrswmkvhdbn')


def reverse_complement_iupac(seq):
    """
    Reverse complement a sequence, supporting IUPAC degenerate bases.

    IUPAC complement rules:
      A<->T, C<->G, R(AG)<->Y(CT), S(GC)<->S(GC), W(AT)<->W(AT),
      K(GT)<->M(AC), B(CGT)<->V(ACG), D(AGT)<->H(ACT), N<->N

    Used for R-end primer matching: R2 reads start with the reverse primer
    oligo sequence (5'->3'), but some databases store reverse primers in
    the reference orientation (reverse complement of the oligo). This
    function allows matching in both orientations.
    """
    return seq.translate(_COMPLEMENT)[::-1]

# ---------------------------------------------------------------------------
# Primer database loading (CONS_F.fas / CONS_R.fas)
# ---------------------------------------------------------------------------

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_DOCS_DIR = os.path.join(os.path.dirname(_SCRIPT_DIR), "docs")


def load_primer_fasta(filepath):
    """Load primers from a FASTA file. Returns list of (name, sequence)."""
    primers = []
    name = None
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                name = line[1:]
            elif name is not None:
                primers.append((name, line.upper()))
                name = None
    return primers


PRIMERS_F = load_primer_fasta(os.path.join(_DOCS_DIR, "CONS_F.fas"))
PRIMERS_R = load_primer_fasta(os.path.join(_DOCS_DIR, "CONS_R.fas"))


def _iupac_compatible(base1, base2):
    """Check if two IUPAC bases are compatible (share at least one base)."""
    set1 = _IUPAC.get(base1.upper(), {base1.upper()})
    set2 = _IUPAC.get(base2.upper(), {base2.upper()})
    return bool(set1 & set2)


def match_primer_database(consensus_window, database, min_identity=0.85,
                          min_informative_fraction=0.50):
    """
    Sliding-window match of the four-state consensus against a primer database.

    database: list of (name, sequence) tuples, e.g. from load_primer_fasta().
              Use PRIMERS_F for forward/R1 reads, PRIMERS_R for reverse/R2.

    For each database primer (and its reverse complement), slides a window
    across the detection-window consensus to find the best matching position.

    The trim position is calculated as: offset + primer_length
    This means everything from position 0 to the end of the matched primer
    region is trimmed, regardless of prefix/suffix mismatches:
      - Database lacks prefix → prefix is trimmed together, no problem
      - Database lacks suffix → a few residual bases stay in read, harmless
      - Exact match → perfect

    N positions in the consensus (from D/T/Q states) are skipped
    as uninformative. At least 50% of positions must be informative for a
    valid match.

    Both orientations are tested so the match works regardless of how
    R-end primers are stored (oligo 5'->3' or reference orientation).

    Returns (primer_name, trim_position, identity) or (None, 0, 0.0).
    """
    if not consensus_window:
        return None, 0, 0.0

    best_name = None
    best_trim = 0
    best_identity = 0.0

    for name, db_seq in database:
        # Try both orientations: original and reverse complement
        candidates = [
            (db_seq, ""),
            (reverse_complement_iupac(db_seq), "_RC"),
        ]
        for seq, suffix in candidates:
            L = len(seq)
            if L > len(consensus_window):
                continue

            # Slide window across consensus
            for offset in range(len(consensus_window) - L + 1):
                segment = consensus_window[offset:offset + L]
                informative = 0
                matches = 0
                for i in range(L):
                    if segment[i] == 'N':
                        continue  # uninformative position, skip
                    informative += 1
                    if _iupac_compatible(segment[i], seq[i]):
                        matches += 1

                if informative < L * min_informative_fraction:
                    continue  # not enough informative positions

                identity = matches / informative
                if identity >= min_identity and identity > best_identity:
                    best_name = f"{name}{suffix}" if suffix else name
                    best_trim = offset + L
                    best_identity = identity

    return best_name, best_trim, best_identity


# ===========================================================================
# Step 6: Trimming helpers
# ===========================================================================

def _out_mode(out_path):
    """Write mode for an output FASTQ ('wt' for gzip, 'w' otherwise)."""
    return "wt" if str(out_path).lower().endswith(".gz") else "w"


def trim_single_file(in_path, out_path, trim_len):
    """Trim first trim_len bases/quals from every read in a FASTQ file."""
    out_mode = _out_mode(out_path)
    count = 0
    with _open_fq(in_path) as fin, _open_fq(out_path, out_mode) as fout:
        for header, seq, plus, qual in _iter_fastq(fin):
            fout.write(f"{header}\n{seq[trim_len:]}\n+\n{qual[trim_len:]}\n")
            count += 1
    return count


def _paired_records(first, second):
    """Read mates together and report mismatched read counts."""
    for one, two in zip_longest(_iter_fastq(first), _iter_fastq(second)):
        if one is None or two is None:
            raise ValueError("Paired FASTQ files contain different numbers of reads")
        yield one, two


def trim_paired_files(r1_in, r2_in, r1_out, r2_out, r1_trim, r2_trim):
    """Trim PE files in lockstep to maintain read pairing."""
    r1_mode = _out_mode(r1_out)
    r2_mode = _out_mode(r2_out)
    count = 0
    with _open_fq(r1_in) as f1i, _open_fq(r2_in) as f2i, \
         _open_fq(r1_out, r1_mode) as f1o, _open_fq(r2_out, r2_mode) as f2o:
        for (h1, s1, p1, q1), (h2, s2, p2, q2) in _paired_records(f1i, f2i):
            f1o.write(f"{h1}\n{s1[r1_trim:]}\n+\n{q1[r1_trim:]}\n")
            f2o.write(f"{h2}\n{s2[r2_trim:]}\n+\n{q2[r2_trim:]}\n")
            count += 1
    return count


def _reverse_orientation(seq, orientation):
    """Use the same minority-versus-other split for detection and trimming."""
    position = orientation['best_pos']
    minority = position < len(seq) and seq[position].upper() == orientation['minority_base']
    return minority if orientation['forward_is_majority'] else not minority


def trim_paired_files_mixed(r1_in, r2_in, r1_out, r2_out, trims, orientation):
    """Trim each original end/orientation before normalizing mate direction."""
    count = swapped = 0
    with _open_fq(r1_in) as f1i, _open_fq(r2_in) as f2i, \
         _open_fq(r1_out, _out_mode(r1_out)) as f1o, \
         _open_fq(r2_out, _out_mode(r2_out)) as f2o:
        for (h1, s1, p1, q1), (h2, s2, p2, q2) in _paired_records(f1i, f2i):
            count += 1
            if _reverse_orientation(s1, orientation):
                swapped += 1
                first_trim, second_trim = trims['r2']['forward'], trims['r1']['reverse']
                f1o.write(f"{h1}\n{s2[first_trim:]}\n+\n{q2[first_trim:]}\n")
                f2o.write(f"{h2}\n{s1[second_trim:]}\n+\n{q1[second_trim:]}\n")
            else:
                first_trim, second_trim = trims['r1']['forward'], trims['r2']['reverse']
                f1o.write(f"{h1}\n{s1[first_trim:]}\n+\n{q1[first_trim:]}\n")
                f2o.write(f"{h2}\n{s2[second_trim:]}\n+\n{q2[second_trim:]}\n")
    return count, swapped


def copy_file(src, dst):
    """Copy a file without modification."""
    shutil.copy2(src, dst)


# ===========================================================================
# Detection and dataset processing
# ===========================================================================

def find_files(input_dir):
    """Discover every sample with the shared pairing rule; select the first."""
    rows = read_layout.discover(input_dir)
    return read_layout.layout(rows), rows[0]['r1'], rows[0]['r2'] or None


def _find_pe_pairs(input_dir):
    rows = read_layout.discover(input_dir)
    if read_layout.layout(rows) != 'PE':
        raise ValueError('Expected paired FASTQ files')
    return ((row['r1'], row['r2']) for row in rows)


def detect_for_reads(reads, label, database, settings=None):
    """Accept database matches first; otherwise apply the strict b1 fold rule."""
    settings = parameters.current() if settings is None else settings
    window = settings['primer.window']
    threshold = settings['primer.fold_threshold']
    unknown_length = settings['primer.unknown_trim_length']
    result = dict(
        detected=False, primer_length=0, consensus='', primer_name='none',
        primer_fold=None, candidate_boundary=0, degenerate_count=0,
        degenerate_density=0.0, fallback_applied=False,
        database_match_identity=0.0, reads_used=len(reads),
        detection_window=window, fold_threshold=threshold,
        support_frequency=settings['primer.support_frequency'],
        database_identity_threshold=settings['primer.database_identity'],
        informative_fraction=settings['primer.informative_fraction'],
        states='', decision='no_valid_reads',
        message='Primer detection failed: no reads passed filters',
    )
    if reads:
        matrix = build_frequency_matrix(reads, window)
        states = classify_all_positions(matrix, settings['primer.support_frequency'])
        consensus = build_consensus_cdv(matrix, states, window)
        name, endpoint, identity = match_primer_database(
            consensus, database, settings['primer.database_identity'],
            settings['primer.informative_fraction'])
        result.update(states=''.join(states), consensus=consensus)
        if name:
            result.update(
                detected=True, primer_length=endpoint, primer_name=name,
                consensus=consensus[:endpoint], database_match_identity=identity,
                decision='known_database_match',
                message=f'Known primer: {name}; trim {endpoint} bp',
            )
        else:
            fold = calculate_primer_fold(states, window)
            degenerate = sum(state != 'C' for state in states)
            result.update(
                primer_fold=fold, candidate_boundary=window,
                degenerate_count=degenerate, degenerate_density=degenerate / window,
                fallback_applied=True,
            )
            if fold < threshold:
                result.update(
                    detected=True, primer_length=unknown_length, primer_name='unknown',
                    decision='unknown_primer',
                    message=f'Unknown primer: fold {fold} < {threshold}; trim {unknown_length} bp',
                )
            else:
                result.update(
                    decision='no_primer',
                    message=f'No primer: fold {fold} >= {threshold}; no trimming',
                )
    print(f"  {label}: {result['message']}", file=sys.stderr)
    return result


def _filter_options(filepath, settings):
    return dict(min_len=settings['primer.min_length'],
                min_avg_qual=settings['primer.min_average_quality'],
                min_complexity=settings['primer.min_complexity'],
                min_entropy=settings['primer.min_entropy'],
                skip_qual=_is_degraded_quality(filepath))


def _filtered_reads(filepath, settings):
    return read_and_filter(filepath, **_filter_options(filepath, settings))


def _detect_mixed_r2(rows, orientation, settings):
    """Group actual mates using R1 orientation, retaining R2's own filters."""
    groups = {'forward': [], 'reverse': []}
    for row in rows:
        options = _filter_options(row['r2'], settings)
        with _open_fq(row['r1']) as first, _open_fq(row['r2']) as second:
            for (_, r1_seq, _, _), (_, seq, _, qual) in _paired_records(first, second):
                if _filter_reason(seq, qual, **options) is None:
                    group = 'forward' if _reverse_orientation(r1_seq, orientation) else 'reverse'
                    groups[group].append((seq, qual))
    return {group: detect_for_reads(groups[group], 'R2/' + group, database, settings)
            for group, database in [('forward', PRIMERS_F), ('reverse', PRIMERS_R)]}


def detect_for_file(filepath, label, database, skip_qual=False, settings=None):
    """Detect from one file, using the same settings as dataset detection."""
    settings = parameters.current() if settings is None else settings
    reads = read_and_filter(
        filepath, min_len=settings['primer.min_length'],
        min_avg_qual=settings['primer.min_average_quality'],
        min_complexity=settings['primer.min_complexity'],
        min_entropy=settings['primer.min_entropy'], skip_qual=skip_qual)
    return detect_for_reads(reads, label, database, settings)


def _primer_entry(result, trim_length=0):
    """Keep detected length separate from the bases actually removed."""
    return {
        'name': result['primer_name'],
        'consensus': result['consensus'],
        'length': result['primer_length'],
        'trim_length': trim_length,
        'detected': result['detected'],
        'fold': result['primer_fold'],
        'decision': result['decision'],
        'reason': result['message'],
        'states': result['states'],
        'reads_used': result['reads_used'],
        'detection_window': result['detection_window'],
        'fold_threshold': result['fold_threshold'],
        'support_frequency': result['support_frequency'],
        'fallback_applied': result['fallback_applied'],
        'database_match_identity': result['database_match_identity'],
        'database_identity_threshold': result['database_identity_threshold'],
        'informative_fraction': result['informative_fraction'],
    }


def trim_single_file_mixed(in_path, out_path, forward_trim, reverse_trim, orientation):
    """Trim SE reads by the orientation group used during detection."""
    count = 0
    with _open_fq(in_path) as source, _open_fq(out_path, _out_mode(out_path)) as dest:
        for header, seq, plus, quality in _iter_fastq(source):
            reverse = _reverse_orientation(seq, orientation)
            trim = reverse_trim if reverse else forward_trim
            dest.write(f'{header}\n{seq[trim:]}\n{plus}\n{quality[trim:]}\n')
            count += 1
    return count


def _detect_orientation_groups(mixed, settings):
    """Use known primers to identify orientation, regardless of group size."""
    forward = detect_for_reads(mixed['majority_reads'], 'Majority/forward', PRIMERS_F, settings)
    reverse = detect_for_reads(mixed['minority_reads'], 'Minority/reverse', PRIMERS_R, settings)
    reverse_base = mixed['minority_base']
    if not all(result['decision'] == 'known_database_match' for result in (forward, reverse)):
        other_forward = detect_for_reads(mixed['minority_reads'], 'Minority/forward', PRIMERS_F, settings)
        other_reverse = detect_for_reads(mixed['majority_reads'], 'Majority/reverse', PRIMERS_R, settings)

        def score(first, second):
            return (sum(result['decision'] == 'known_database_match' for result in (first, second)),
                    first['database_match_identity'] + second['database_match_identity'])

        if score(other_forward, other_reverse) > score(forward, reverse):
            forward, reverse = other_forward, other_reverse
            reverse_base = mixed['majority_base']
    return forward, reverse, reverse_base


def _save_info(output_dir, info):
    path = os.path.join(output_dir, 'primer_info.json')
    with open(path, 'w') as handle:
        json.dump(info, handle, indent=2)
        handle.write('\n')
    print(f'  Primer info saved to: {path}', file=sys.stderr)


def process_dataset(input_dir, output_dir, settings, detect_only=False,
                    mixed_orientation=False):
    """Detect once per dataset, and decide whether to skip before writing reads."""
    info = {'mode': None, 'layout': None, 'status': 'failed', 'reason': '',
            'detect_only': detect_only,
            'skip_unknown_primers': settings['primer.skip_unknown'],
            'parameters': parameters.nested({k: v for k, v in settings.items()
                                             if k.startswith('primer.')})['primer']}
    try:
        if os.path.realpath(input_dir) == os.path.realpath(output_dir):
            raise ValueError('Primer input and output directories must be different')
        rows = read_layout.discover(input_dir)
        mode = read_layout.layout(rows)
        info.update(mode=mode, layout=mode, detection_sample=rows[0]['sample'])
        selected = [row for row in rows if row['sample'] == rows[0]['sample']]
        reads = {'r1': [], 'r2': []}
        for row in selected:
            reads['r1'].extend(_filtered_reads(row['r1'], settings))

        mixed = None
        if reads['r1'] and (mode == 'PE' or detect_only or mixed_orientation):
            mixed = detect_mixed_orientation(
                build_frequency_matrix(reads['r1'], settings['primer.window']),
                reads['r1'])
        if mixed:
            info['mode'] = mode + '_mixed'
            info['orientation'] = {key: mixed[key] for key in
                                   ('best_pos', 'majority_base', 'minority_base', 'ratio')}
            forward, reverse, reverse_base = _detect_orientation_groups(mixed, settings)
            info['orientation']['reverse_base'] = reverse_base
            info['orientation']['forward_is_majority'] = reverse_base == mixed['minority_base']
        else:
            if mode == 'PE':
                for row in selected:
                    reads['r2'].extend(_filtered_reads(row['r2'], settings))
            forward = detect_for_reads(reads['r1'], 'R1' if mode == 'PE' else 'SE',
                                       PRIMERS_F, settings)
            reverse = (detect_for_reads(reads['r2'], 'R2', PRIMERS_R, settings)
                       if mode == 'PE' else None)
        info['forward_primer'] = _primer_entry(forward)
        if reverse is not None:
            info['reverse_primer'] = _primer_entry(reverse)
        results = [forward] + ([reverse] if reverse is not None else [])
        if mixed and mode == 'PE':
            r2_results = _detect_mixed_r2(selected, info['orientation'], settings)
            info['original_ends'] = {
                'r1': {'forward': info['forward_primer'], 'reverse': info['reverse_primer']},
                'r2': {group: _primer_entry(result) for group, result in r2_results.items()},
            }
            results.extend(r2_results.values())
        if any(result['decision'] == 'no_valid_reads' for result in results):
            info['reason'] = 'Primer detection failed: no valid reads in a required end/orientation'
            _save_info(output_dir, info)
            print(f"ERROR: {info['reason']}", file=sys.stderr)
            return 1
        if settings['primer.skip_unknown'] and any(
                result['decision'] == 'unknown_primer' for result in results):
            info.update(status='skipped', reason='Unknown primer detected; --skip-unknown-primers enabled')
            _save_info(output_dir, info)
            print(f"SKIP: {info['reason']}", file=sys.stderr)
            return 98

        forward_trim = 0 if detect_only else forward['primer_length']
        reverse_trim = 0 if detect_only or reverse is None else reverse['primer_length']
        if mixed and mode == 'PE':
            trims = {end: {group: 0 if detect_only else entry['length']
                           for group, entry in groups.items()}
                     for end, groups in info['original_ends'].items()}
        for row in rows:
            first = row['r1']
            first_out = os.path.join(output_dir, os.path.basename(first))
            if detect_only or (forward_trim == 0 and reverse_trim == 0 and mixed is None):
                for direction in ('r1', 'r2') if mode == 'PE' else ('r1',):
                    source = row[direction]
                    copy_file(source, os.path.join(output_dir, os.path.basename(source)))
            elif mode == 'PE':
                second = row['r2']
                second_out = os.path.join(output_dir, os.path.basename(second))
                if mixed:
                    trim_paired_files_mixed(first, second, first_out, second_out,
                                            trims, info['orientation'])
                else:
                    trim_paired_files(first, second, first_out, second_out,
                                      forward_trim, reverse_trim)
            elif mixed:
                trim_single_file_mixed(first, first_out, forward_trim, reverse_trim,
                                       info['orientation'])
            else:
                trim_single_file(first, first_out, forward_trim)
        info['forward_primer']['trim_length'] = forward_trim
        if reverse is not None:
            info['reverse_primer']['trim_length'] = reverse_trim
        if mixed and mode == 'PE':
            for end, groups in info['original_ends'].items():
                for group, entry in groups.items():
                    entry['trim_length'] = trims[end][group]
        info.update(status='completed', reason='Detection only; reads copied unchanged' if detect_only
                    else 'Primer decisions applied to dataset')
        _save_info(output_dir, info)
        return 0
    except (ValueError, OSError) as error:
        info.update(status='failed', reason=str(error))
        _save_info(output_dir, info)
        print(f'ERROR: {error}', file=sys.stderr)
        return 1


def main():
    parser = argparse.ArgumentParser(
        description='Database-first primer detection with the b1 four-state fallback')
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('--detect-only', action='store_true',
                        help='Detect primers and copy every FASTQ unchanged; also check SE orientations')
    parser.add_argument('--mixed-orientation', action='store_true',
                        help='Check SE orientations and trim each group separately (PacBio vsearch)')
    parser.add_argument('--skip-unknown-primers', action='store_true', default=None,
                        help='Skip the dataset when any end/orientation has an unknown primer')
    parser.add_argument('--parameter', help='JSON file overriding built-in parameter defaults')
    args = parser.parse_args()
    try:
        settings = parameters.load(args.parameter) if args.parameter else parameters.current()
    except (ValueError, OSError) as error:
        parser.exit(2, f'Parameter error: {error}\n')
    if args.skip_unknown_primers is not None:
        settings['primer.skip_unknown'] = args.skip_unknown_primers
    os.makedirs(args.output, exist_ok=True)
    return process_dataset(os.path.abspath(args.input), os.path.abspath(args.output),
                           settings, args.detect_only, args.mixed_orientation)


if __name__ == '__main__':
    sys.exit(main())
