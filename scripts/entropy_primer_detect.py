#!/usr/bin/env python3
"""
Three-state (C/D/V) Primer Detection and Removal
(with mixed R1/R2 orientation support)

Principle:
  Analogous to FastQC "Per Base Sequence Content": build a position-wise
  base frequency matrix and classify each position into one of three states:
    C (Conserved) — one dominant base (f1 ≥ 0.80)
    D (Degenerate) — two dominant bases (f1+f2 ≥ 0.80, f2 ≥ 0.10)
    V (Variable)  — no clear dominant pattern (biological region)
  Primer = longest [CD]* prefix from position 0, terminated at first V.

Algorithm:
  Step 0:   Scan directory, pick first sample (PE or SE).
  Step 1:   Read & pre-filter ALL R1 reads (length, quality, complexity).
  Step 2:   Build 60×4 position-wise base frequency matrix.
  Phase 1:  Mixed-orientation detection (PE only).
            Check first 15 positions for bimodal base distributions.
            If >50% bimodal → mixed orientation detected.

  --- Branch A (mixed orientation) ---
  Phase 2:  Split sampled reads into majority/minority groups by the base
            at the best bimodal position. Run CDV detection (Steps 3-5)
            on each group independently to determine forward/reverse primer
            lengths.
  Phase 3:  Trim ALL files with per-read orientation correction.
            For each read pair, classify by a single-character check at
            best_pos. Flipped reads have R1↔R2 swapped before trimming.

  --- Branch B (single orientation) ---
  Step 3:   Classify each position as C / D / V.
  Step 4:   Find primer boundary (longest [CD]* prefix, with noise tolerance).
  Step 5:   Build consensus (C→dominant base, D→N) and report.
  Step 6:   Trim ALL files using detected primer length.

  Constraints:
    - Minimum primer length: 10 bp (V before pos 10 → no primer).
    - Maximum primer length: 30 bp.
    - Noise tolerance: up to 2 consecutive V positions surrounded by
      ≥3 non-V positions on each side are treated as degenerate primer
      bases (e.g., IUPAC N/V/W) and skipped.
"""

import sys
import os
import gzip
import glob
import math
import shutil
import argparse
import json
from collections import Counter


# ===========================================================================
# Low-level FASTQ helpers (avoid BioPython for speed in trimming step)
# ===========================================================================

def _open_fq(filepath, mode="rt"):
    """Open FASTQ, transparently handling gzip."""
    if filepath.endswith(".gz"):
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
        if not header:
            break
        yield header, seq, plus, qual


# ===========================================================================
# Step 1: Read & pre-filter
# ===========================================================================

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


def read_and_filter(filepath, min_len=50,
                    min_avg_qual=20, min_complexity=0.3,
                    min_entropy=1.0):
    """
    Step 1: Read ALL reads from a FASTQ file, applying quality filters.
    Uses the entire first sample for maximum species diversity in the
    frequency matrix, which improves primer boundary accuracy.
    Filters:
      - Length >= min_len
      - Average Phred quality >= min_avg_qual
      - K-mer complexity >= min_complexity (unique 2-mers / 16)
      - Shannon entropy >= min_entropy (base composition diversity)
    Returns list of (seq_string, qual_string) tuples.
    """
    reads = []
    total = 0
    disc_len = disc_qual = disc_cplx = disc_entropy = 0

    with _open_fq(filepath) as fh:
        for header, seq, plus, qual in _iter_fastq(fh):
            total += 1

            if len(seq) < min_len:
                disc_len += 1
                continue
            if _avg_quality(qual) < min_avg_qual:
                disc_qual += 1
                continue
            if _sequence_complexity(seq) < min_complexity:
                disc_cplx += 1
                continue
            if _read_entropy(seq) < min_entropy:
                disc_entropy += 1
                continue

            reads.append((seq, qual))

    print(f"  Scanned {total} reads, kept {len(reads)}", file=sys.stderr)
    print(f"  Discarded: {disc_len} (length<{min_len}), "
          f"{disc_qual} (avgQ<{min_avg_qual}), "
          f"{disc_cplx} (complexity<{min_complexity}), "
          f"{disc_entropy} (entropy<{min_entropy})", file=sys.stderr)

    if len(reads) < 500:
        print(f"  WARNING: Only {len(reads)} reads passed filters (< 500)",
              file=sys.stderr)

    return reads


# ===========================================================================
# Step 2: Position-wise base frequency matrix
# ===========================================================================

def build_frequency_matrix(reads, num_positions=60):
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
        majority_base  – dominant base at best_pos (forward-primer reads)
        minority_base  – second base at best_pos (reverse-primer reads)
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
# Step 3: Three-state position classification (C / D / V)
# ===========================================================================

BASES = ['A', 'C', 'G', 'T']


def classify_position(freqs):
    """
    Classify a single position based on its base frequency distribution.

    C (Conserved):  one dominant base, f1 >= 0.80
    D (Degenerate): two dominant bases, f1+f2 >= 0.80 and f2 >= 0.10
    V (Variable):   no clear dominant pattern (biological region)

    Analogous to FastQC "Per Base Sequence Content" interpretation.
    """
    sorted_f = sorted(freqs, reverse=True)
    f1, f2 = sorted_f[0], sorted_f[1]
    if f1 >= 0.80:
        return 'C'
    if f1 + f2 >= 0.80 and f2 >= 0.10:
        return 'D'
    return 'V'


def classify_all_positions(freq_matrix):
    """Classify every position in the frequency matrix as C, D, or V."""
    return [classify_position(freqs) for freqs in freq_matrix]


# ===========================================================================
# Step 4: Primer boundary detection
# ===========================================================================

def find_primer_boundary(states, min_len=10, max_len=30, max_v_run=2):
    """
    Primer = longest [CD]* prefix from position 0, terminated at first V.

    Constraints:
      - Minimum length 10 bp (V before pos 10 → no primer).
      - Maximum length 30 bp.
      - Noise tolerance: up to max_v_run (default 2) consecutive V positions
        are skipped if surrounded by >= 3 non-V positions on each side.
        This handles 3-way/4-way degenerate primer bases (e.g., N/V/W in
        IUPAC) that look variable by frequency alone.

    Returns primer_length (0 if no primer detected).
    """
    n = min(max_len, len(states))
    end = 0
    i = 0

    while i < n:
        if states[i] in ('C', 'D'):
            end = i + 1
            i += 1
        elif states[i] == 'V':
            # Count consecutive V's starting at position i
            v_start = i
            while i < len(states) and states[i] == 'V':
                i += 1
            v_count = i - v_start

            # Tolerance: up to max_v_run consecutive V's, if surrounded
            # by >= 3 non-V positions on each side
            if v_count <= max_v_run:
                before_ok = (v_start >= 3 and
                             all(s != 'V'
                                 for s in states[v_start - 3:v_start]))
                after_ok = (i + 2 < len(states) and
                            all(s != 'V'
                                for s in states[i:i + 3]))
                if before_ok and after_ok:
                    # Degenerate primer positions — treat as noise, skip
                    end = i
                    continue

            # Real variable region — primer ends before these V's
            break
        else:
            i += 1

    if end < min_len:
        return 0  # Too short to be a primer

    return end


# ===========================================================================
# Step 5: Consensus sequence
# ===========================================================================

def build_consensus_cdv(freq_matrix, states, primer_length):
    """
    Build consensus sequence for the detected primer region.

    C positions → dominant base (A/C/G/T)
    D positions → N (degenerate)
    Tolerated V → N
    """
    consensus = []
    for i in range(primer_length):
        if states[i] == 'C':
            max_idx = freq_matrix[i].index(max(freq_matrix[i]))
            consensus.append(BASES[max_idx])
        else:  # D or tolerated V
            consensus.append('N')
    return ''.join(consensus)


# ===========================================================================
# Known 16S primer database (boundary refinement)
# ===========================================================================
# CDV classification cannot distinguish primer-conserved from
# genomic-conserved positions. After CDV detects a primer exists,
# match the consensus against known 16S primers to get the exact length.

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


def match_primer_database(consensus_30bp, database, min_identity=0.85):
    """
    Sliding-window match of CDV consensus against a primer database.

    database: list of (name, sequence) tuples, e.g. from load_primer_fasta().
              Use PRIMERS_F for forward/R1 reads, PRIMERS_R for reverse/R2.

    For each database primer (and its reverse complement), slides a window
    across the 30bp consensus to find the best matching position.

    The trim position is calculated as: offset + primer_length
    This means everything from position 0 to the end of the matched primer
    region is trimmed, regardless of prefix/suffix mismatches:
      - Database lacks prefix → prefix is trimmed together, no problem
      - Database lacks suffix → a few residual bases stay in read, harmless
      - Exact match → perfect

    N positions in the consensus (from CDV D/V classification) are skipped
    as uninformative. At least 50% of positions must be informative for a
    valid match.

    Both orientations are tested so the match works regardless of how
    R-end primers are stored (oligo 5'->3' or reference orientation).

    Returns (primer_name, trim_position, identity) or (None, 0, 0.0).
    """
    if not consensus_30bp:
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
            if L > len(consensus_30bp):
                continue

            # Slide window across consensus
            for offset in range(len(consensus_30bp) - L + 1):
                segment = consensus_30bp[offset:offset + L]
                informative = 0
                matches = 0
                for i in range(L):
                    if segment[i] == 'N':
                        continue  # uninformative position, skip
                    informative += 1
                    if _iupac_compatible(segment[i], seq[i]):
                        matches += 1

                if informative < L * 0.5:
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

def trim_single_file(in_path, out_path, trim_len):
    """Trim first trim_len bases/quals from every read in a FASTQ file."""
    out_mode = "wt" if out_path.endswith(".gz") else "w"
    count = 0
    with _open_fq(in_path) as fin, _open_fq(out_path, out_mode) as fout:
        for header, seq, plus, qual in _iter_fastq(fin):
            fout.write(f"{header}\n{seq[trim_len:]}\n+\n{qual[trim_len:]}\n")
            count += 1
    return count


def trim_paired_files(r1_in, r2_in, r1_out, r2_out, r1_trim, r2_trim):
    """Trim PE files in lockstep to maintain read pairing."""
    r1_mode = "wt" if r1_out.endswith(".gz") else "w"
    r2_mode = "wt" if r2_out.endswith(".gz") else "w"
    count = 0
    with _open_fq(r1_in) as f1i, _open_fq(r2_in) as f2i, \
         _open_fq(r1_out, r1_mode) as f1o, _open_fq(r2_out, r2_mode) as f2o:
        for (h1, s1, p1, q1), (h2, s2, p2, q2) in zip(
                _iter_fastq(f1i), _iter_fastq(f2i)):
            f1o.write(f"{h1}\n{s1[r1_trim:]}\n+\n{q1[r1_trim:]}\n")
            f2o.write(f"{h2}\n{s2[r2_trim:]}\n+\n{q2[r2_trim:]}\n")
            count += 1
    return count


def trim_paired_files_mixed(r1_in, r2_in, r1_out, r2_out,
                            r1_trim, r2_trim, best_pos, majority_base):
    """
    Trim PE files with per-read orientation correction.

    For each read pair, checks R1's base at best_pos:
      - majority_base  -> normal:  R1 trimmed by r1_trim, R2 by r2_trim
      - other base     -> flipped: swap sequences, then trim

    This applies the sampled detection result to ALL reads in the file.
    Classification cost per read: one character comparison at best_pos.

    Returns (total_count, swapped_count).
    """
    r1_mode = "wt" if r1_out.endswith(".gz") else "w"
    r2_mode = "wt" if r2_out.endswith(".gz") else "w"
    count = 0
    swapped = 0
    with _open_fq(r1_in) as f1i, _open_fq(r2_in) as f2i, \
         _open_fq(r1_out, r1_mode) as f1o, \
         _open_fq(r2_out, r2_mode) as f2o:
        for (h1, s1, p1, q1), (h2, s2, p2, q2) in zip(
                _iter_fastq(f1i), _iter_fastq(f2i)):
            count += 1
            # Classify: if R1 base at best_pos != majority → flipped
            is_flipped = (best_pos < len(s1)
                          and s1[best_pos].upper() != majority_base)
            if is_flipped:
                # Swap: output R1 = input R2, output R2 = input R1
                swapped += 1
                f1o.write(f"{h1}\n{s2[r1_trim:]}\n+\n{q2[r1_trim:]}\n")
                f2o.write(f"{h2}\n{s1[r2_trim:]}\n+\n{q1[r2_trim:]}\n")
            else:
                f1o.write(f"{h1}\n{s1[r1_trim:]}\n+\n{q1[r1_trim:]}\n")
                f2o.write(f"{h2}\n{s2[r2_trim:]}\n+\n{q2[r2_trim:]}\n")
    return count, swapped


def copy_file(src, dst):
    """Copy a file without modification."""
    shutil.copy2(src, dst)


# ===========================================================================
# Step 0: File discovery
# ===========================================================================

def find_files(input_dir):
    """
    Identify SE/PE mode and pick the first sample.
    Returns: ("PE", first_r1, first_r2) or ("SE", first_file, None)
    """
    # PE: _R1/_R2
    for r1 in sorted(glob.glob(os.path.join(input_dir, "*_R1*.fastq*"))):
        r2 = r1.replace("_R1", "_R2")
        if os.path.exists(r2):
            return "PE", r1, r2

    # PE: _1/_2
    for r1 in sorted(glob.glob(os.path.join(input_dir, "*_1.fastq*"))):
        r2 = r1.replace("_1.fastq", "_2.fastq")
        if os.path.exists(r2):
            return "PE", r1, r2

    # SE: everything except R1/R2
    all_fq = sorted(glob.glob(os.path.join(input_dir, "*.fastq*")))
    se = [f for f in all_fq if "_R1" not in f and "_R2" not in f]
    if se:
        return "SE", se[0], None
    if all_fq:
        return "SE", all_fq[0], None

    return None, None, None


# ===========================================================================
# Detection pipeline (Steps 1-6) for one file
# ===========================================================================

def detect_for_reads(reads, label, database):
    """
    Run three-state (C/D/V) primer detection on pre-filtered reads.

    database: list of (name, sequence) tuples from load_primer_fasta().
              Use PRIMERS_F for forward/R1, PRIMERS_R for reverse/R2.

    Steps:
      2. Build 60×4 position-wise base frequency matrix.
      3. Classify each position as C (Conserved), D (Degenerate), V (Variable).
      4. Find primer boundary: longest [CD]* prefix with noise tolerance.
      5. Build consensus (C→dominant base, D/V→N), match against database.
    """
    print(f"\n{'=' * 60}", file=sys.stderr)
    print(f"  Analyzing {label}: {len(reads)} reads", file=sys.stderr)
    print(f"{'=' * 60}", file=sys.stderr)

    if not reads:
        print("  ERROR: No reads to analyze", file=sys.stderr)
        return dict(detected=False, primer_length=0, consensus="",
                    message="No primer detected (no reads)")

    # Step 2: Build frequency matrix
    print(f"\n[Step 2] Building frequency matrix (60 positions) ...",
          file=sys.stderr)
    freq_matrix = build_frequency_matrix(reads, num_positions=60)

    # Step 3: Classify each position as C/D/V
    print(f"\n[Step 3] Classifying positions (C/D/V) ...", file=sys.stderr)
    states = classify_all_positions(freq_matrix)
    state_str = ''.join(states)
    print(f"  Pos  0-29: {state_str[:30]}", file=sys.stderr)
    print(f"  Pos 30-59: {state_str[30:]}", file=sys.stderr)

    # Per-position detail for the first 30 positions (max primer range)
    for i in range(min(30, len(freq_matrix))):
        freqs = freq_matrix[i]
        paired = sorted(zip(BASES, freqs), key=lambda x: -x[1])
        top = '  '.join(f"{b}:{f:.2f}" for b, f in paired[:2])
        print(f"    [{i:2d}] {states[i]}  {top}", file=sys.stderr)

    # Step 4: Find primer boundary (CDV — used as fallback)
    print(f"\n[Step 4] Finding CDV primer boundary "
          f"(min=10bp, max=30bp) ...", file=sys.stderr)
    cdv_boundary = find_primer_boundary(states, min_len=10, max_len=30)
    print(f"  CDV boundary: {cdv_boundary} bp", file=sys.stderr)

    # D density: fraction of D (Degenerate) positions in the CDV region.
    # Real synthetic primers have ≤ 20% degenerate bases (max: 341F = 3/15).
    # Biological 16S sequence shows higher D density due to inter-species
    # variation.  If D density > 25%, the "primer" region is actually
    # genomic sequence → primer was already removed upstream.
    if cdv_boundary > 0:
        d_count = sum(1 for s in states[:cdv_boundary] if s == 'D')
        d_density = d_count / cdv_boundary
        print(f"  D density: {d_count}/{cdv_boundary} = {d_density:.1%}",
              file=sys.stderr)
    else:
        d_density = 0.0

    # Step 5: Build 30bp consensus, then match against primer database
    print(f"\n[Step 5] Database matching (sliding window) ...",
          file=sys.stderr)
    consensus_30 = build_consensus_cdv(freq_matrix, states, 30)
    print(f"  30bp consensus: {consensus_30}", file=sys.stderr)

    # Sliding window: tries both orientations (original + RC)
    db_name, db_trim, db_identity = match_primer_database(consensus_30,
                                                          database)

    if db_name:
        # Database match → use database-defined trim position
        is_rc = db_name.endswith("_RC")
        orient_msg = " (matched via reverse complement)" if is_rc else ""
        print(f"  Database match: {db_name} "
              f"(trim={db_trim}bp, identity={db_identity:.2f})"
              f"{orient_msg}", file=sys.stderr)
        if cdv_boundary > 0 and db_trim != cdv_boundary:
            print(f"  CDV boundary={cdv_boundary}bp -> "
                  f"overridden by database trim={db_trim}bp",
                  file=sys.stderr)
        primer_length = db_trim
        consensus = consensus_30[:db_trim]
        primer_name = db_name
        result = dict(detected=True, primer_length=primer_length,
                      consensus=consensus, primer_name=primer_name,
                      message=f"Primer detected: {primer_name} "
                              f"(trim={primer_length}bp)")
    elif cdv_boundary > 0 and d_density <= 0.25:
        # No database match, but D density consistent with a real primer
        # → fall back to CDV boundary
        print(f"  No database match, using CDV boundary: "
              f"{cdv_boundary}bp", file=sys.stderr)
        primer_length = cdv_boundary
        consensus = consensus_30[:cdv_boundary]
        result = dict(detected=True, primer_length=primer_length,
                      consensus=consensus, primer_name="unknown",
                      message=f"Primer detected: unknown "
                              f"(trim={primer_length}bp, CDV fallback)")
    elif cdv_boundary > 0 and d_density > 0.25:
        # No database match AND D density too high for a real primer.
        # Real primers have ≤ 20% degenerate bases; this region looks like
        # conserved 16S gene sequence → primer was already removed upstream.
        print(f"  No database match; D density {d_density:.0%} > 25% "
              f"→ not a primer (likely already trimmed)", file=sys.stderr)
        primer_length = 0
        result = dict(detected=False, primer_length=0, consensus="",
                      primer_name="none",
                      message=f"No primer detected "
                              f"(D density {d_density:.0%} > 25%)")
    else:
        # Neither database nor CDV detected a primer
        primer_length = 0
        result = dict(detected=False, primer_length=0, consensus="",
                      message="No primer detected")

    print(f"  1. Primer exists:     "
          f"{'Yes' if result['detected'] else 'No'}", file=sys.stderr)
    print(f"  2. Truncation length: {result['primer_length']} bp",
          file=sys.stderr)
    print(f"  3. Consensus:         "
          f"{result['consensus'] if result['consensus'] else 'N/A'}",
          file=sys.stderr)

    return result


def detect_for_file(filepath, label, database):
    """Run Steps 1-6 on a single FASTQ file. Returns result dict."""
    print(f"\n[Step 1] Reading and filtering reads from "
          f"{os.path.basename(filepath)} ...", file=sys.stderr)
    reads = read_and_filter(filepath)
    if not reads:
        print("  ERROR: No reads passed filters", file=sys.stderr)
        return dict(detected=False, primer_length=0, consensus="",
                    message="No primer detected (no reads)")
    return detect_for_reads(reads, label, database)


# ===========================================================================
# Main
# ===========================================================================

def _find_pe_pairs(input_dir):
    """Find all PE file pairs in input_dir. Yields (r1_path, r2_path)."""
    # Try _R1/_R2 pattern first
    found = False
    for r1 in sorted(glob.glob(os.path.join(input_dir, "*_R1*.fastq*"))):
        r2 = r1.replace("_R1", "_R2")
        if os.path.exists(r2):
            found = True
            yield r1, r2
    if found:
        return
    # Fallback: _1/_2 pattern
    for r1 in sorted(glob.glob(os.path.join(input_dir, "*_1.fastq*"))):
        r2 = r1.replace("_1.fastq", "_2.fastq")
        if os.path.exists(r2):
            yield r1, r2


def main():
    parser = argparse.ArgumentParser(
        description="Three-state (C/D/V) primer detection and removal")
    parser.add_argument("-i", "--input", required=True,
                        help="Input directory containing FASTQ files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for trimmed files")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60, file=sys.stderr)
    print("THREE-STATE (C/D/V) PRIMER DETECTION", file=sys.stderr)
    print("=" * 60, file=sys.stderr)

    # ------------------------------------------------------------------
    # Step 0: File discovery
    # ------------------------------------------------------------------
    print(f"\n[Step 0] Scanning directory: {input_dir}", file=sys.stderr)
    mode, first_file, second_file = find_files(input_dir)

    if mode is None:
        print("ERROR: No FASTQ files found in input directory",
              file=sys.stderr)
        sys.exit(1)

    print(f"  Mode: {mode}", file=sys.stderr)
    print(f"  First file:  {os.path.basename(first_file)}", file=sys.stderr)
    if second_file:
        print(f"  Second file: {os.path.basename(second_file)}",
              file=sys.stderr)

    # ------------------------------------------------------------------
    # Step 1: Read & filter first R1 sample
    # ------------------------------------------------------------------
    print(f"\n[Step 1] Reading and filtering reads from "
          f"{os.path.basename(first_file)} ...", file=sys.stderr)
    r1_reads = read_and_filter(first_file)
    if not r1_reads:
        print("  ERROR: No reads passed filters", file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 2: Build frequency matrix (used for both CDV classification
    #         and mixed-orientation detection)
    # ------------------------------------------------------------------
    print(f"\n[Step 2] Building frequency matrix (60 positions) ...",
          file=sys.stderr)
    freq_matrix = build_frequency_matrix(r1_reads, num_positions=60)

    # ------------------------------------------------------------------
    # Phase 1: Mixed-orientation detection (PE only)
    # ------------------------------------------------------------------
    mixed_info = None
    if mode == "PE":
        mixed_info = detect_mixed_orientation(freq_matrix, r1_reads)

    # ==================================================================
    # Branch A: Mixed orientation — detect from majority/minority,
    #           trim all files with per-read orientation correction
    # ==================================================================
    if mixed_info is not None:
        # Phase 2: CDV detection on each orientation group
        r1_result = detect_for_reads(mixed_info['majority_reads'],
                                     "R1 majority (forward primer)",
                                     PRIMERS_F)
        r2_result = detect_for_reads(mixed_info['minority_reads'],
                                     "R1 minority (reverse primer)",
                                     PRIMERS_R)

        r1_trim = r1_result['primer_length']
        r2_trim = r2_result['primer_length']

        # Summary
        print(f"\n{'=' * 60}", file=sys.stderr)
        print("DETECTION SUMMARY (Mixed Orientation)", file=sys.stderr)
        print(f"{'=' * 60}", file=sys.stderr)
        print(f"  Forward primer (majority R1): {r1_result['message']}",
              file=sys.stderr)
        print(f"  Reverse primer (minority R1): {r2_result['message']}",
              file=sys.stderr)
        print(f"  R1 trim: {r1_trim} bp, R2 trim: {r2_trim} bp",
              file=sys.stderr)

        # Phase 3: Trim all PE files with orientation correction
        print(f"\n{'=' * 60}", file=sys.stderr)
        print("TRIMMING (with orientation correction)", file=sys.stderr)
        print(f"{'=' * 60}", file=sys.stderr)

        best_pos = mixed_info['best_pos']
        majority_base = mixed_info['majority_base']
        total_pairs = 0
        total_swapped = 0
        pairs_found = False

        for r1, r2 in _find_pe_pairs(input_dir):
            pairs_found = True
            r1_out = os.path.join(output_dir, os.path.basename(r1))
            r2_out = os.path.join(output_dir, os.path.basename(r2))
            n, s = trim_paired_files_mixed(
                r1, r2, r1_out, r2_out,
                r1_trim, r2_trim, best_pos, majority_base)
            total_pairs += n
            total_swapped += s
            print(f"    {os.path.basename(r1)} + "
                  f"{os.path.basename(r2)}: "
                  f"{n} pairs ({s} swapped)", file=sys.stderr)

        if not pairs_found:
            print("  ERROR: No PE file pairs found for trimming",
                  file=sys.stderr)
            sys.exit(1)

        print(f"\n  Total: {total_pairs} pairs processed, "
              f"{total_swapped} orientation-corrected "
              f"({total_swapped/max(1,total_pairs)*100:.1f}%)",
              file=sys.stderr)

    # ==================================================================
    # Branch B: Normal single orientation
    # ==================================================================
    else:
        # Steps 3-6 on already-loaded R1 reads (avoid re-reading)
        r1_result = detect_for_reads(r1_reads,
                                     "R1" if mode == "PE" else "SE",
                                     PRIMERS_F)

        r2_result = None
        if mode == "PE" and second_file:
            r2_result = detect_for_file(second_file, "R2", PRIMERS_R)

        # Summary
        print(f"\n{'=' * 60}", file=sys.stderr)
        print("DETECTION SUMMARY", file=sys.stderr)
        print(f"{'=' * 60}", file=sys.stderr)

        if mode == "PE":
            print(f"  R1: {r1_result['message']}", file=sys.stderr)
            print(f"  R2: {r2_result['message']}", file=sys.stderr)
            r1_trim = r1_result['primer_length']
            r2_trim = r2_result['primer_length']
            any_detected = (r1_result['detected']
                            or r2_result['detected'])
        else:
            print(f"  SE: {r1_result['message']}", file=sys.stderr)
            r1_trim = r1_result['primer_length']
            r2_trim = 0
            any_detected = r1_result['detected']

        # Step 6: Trim the ENTIRE dataset
        print(f"\n{'=' * 60}", file=sys.stderr)
        print("TRIMMING", file=sys.stderr)
        print(f"{'=' * 60}", file=sys.stderr)

        if not any_detected:
            print("\n  No primers detected. Copying files unchanged ...",
                  file=sys.stderr)
            for f in glob.glob(os.path.join(input_dir, "*.fastq*")):
                copy_file(f, os.path.join(output_dir,
                                          os.path.basename(f)))
            print("  Done.", file=sys.stderr)
            return

        # --- SE trimming ---
        if mode == "SE":
            print(f"\n  Trimming all SE files: "
                  f"remove first {r1_trim} bp ...", file=sys.stderr)
            for f in sorted(glob.glob(os.path.join(input_dir,
                                                    "*.fastq*"))):
                out = os.path.join(output_dir, os.path.basename(f))
                n = trim_single_file(f, out, r1_trim)
                print(f"    {os.path.basename(f)}: {n} reads trimmed",
                      file=sys.stderr)

        # --- PE trimming ---
        else:
            print(f"\n  Trimming all PE pairs: "
                  f"R1 -{r1_trim} bp, R2 -{r2_trim} bp ...",
                  file=sys.stderr)
            pairs_found = False
            for r1, r2 in _find_pe_pairs(input_dir):
                pairs_found = True
                r1_out = os.path.join(output_dir, os.path.basename(r1))
                r2_out = os.path.join(output_dir, os.path.basename(r2))
                n = trim_paired_files(r1, r2, r1_out, r2_out,
                                      r1_trim, r2_trim)
                print(f"    {os.path.basename(r1)} + "
                      f"{os.path.basename(r2)}: {n} pairs trimmed",
                      file=sys.stderr)

            if not pairs_found:
                print("  ERROR: No PE file pairs found for trimming",
                      file=sys.stderr)
                sys.exit(1)

    # Save detected primer info to JSON for downstream tools (e.g. DADA2)
    primer_info = {}
    if mixed_info is not None:
        # Branch A: mixed orientation PE
        primer_info = {
            "mode": "PE_mixed",
            "forward_primer": {
                "name": r1_result.get('primer_name', 'unknown'),
                "consensus": r1_result.get('consensus', ''),
                "length": r1_result.get('primer_length', 0),
                "detected": r1_result.get('detected', False),
            },
            "reverse_primer": {
                "name": r2_result.get('primer_name', 'unknown'),
                "consensus": r2_result.get('consensus', ''),
                "length": r2_result.get('primer_length', 0),
                "detected": r2_result.get('detected', False),
            },
        }
    elif mode == "PE":
        primer_info = {
            "mode": "PE",
            "forward_primer": {
                "name": r1_result.get('primer_name', 'unknown'),
                "consensus": r1_result.get('consensus', ''),
                "length": r1_result.get('primer_length', 0),
                "detected": r1_result.get('detected', False),
            },
            "reverse_primer": {
                "name": r2_result.get('primer_name', 'unknown') if r2_result else 'unknown',
                "consensus": r2_result.get('consensus', '') if r2_result else '',
                "length": r2_result.get('primer_length', 0) if r2_result else 0,
                "detected": r2_result.get('detected', False) if r2_result else False,
            },
        }
    else:
        primer_info = {
            "mode": "SE",
            "forward_primer": {
                "name": r1_result.get('primer_name', 'unknown'),
                "consensus": r1_result.get('consensus', ''),
                "length": r1_result.get('primer_length', 0),
                "detected": r1_result.get('detected', False),
            },
        }

    primer_info_path = os.path.join(output_dir, "primer_info.json")
    with open(primer_info_path, 'w') as f:
        json.dump(primer_info, f, indent=2)
    print(f"\n  Primer info saved to: {primer_info_path}", file=sys.stderr)

    print(f"\n  Output directory: {output_dir}", file=sys.stderr)
    print("  Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
