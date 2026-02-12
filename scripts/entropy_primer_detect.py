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
  Step 1:   Read & pre-filter 2000 R1 reads (length, quality, complexity).
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
    - Noise tolerance: 1 isolated V surrounded by ≥3 non-V positions
      on each side is treated as sequencing noise and skipped.
"""

import sys
import os
import gzip
import glob
import math
import shutil
import argparse
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


def read_and_filter(filepath, max_reads=2000, min_len=50,
                    min_avg_qual=20, min_complexity=0.3,
                    min_entropy=1.0):
    """
    Step 1: Read FASTQ, apply quality filters, return up to max_reads reads.
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
            if len(reads) >= max_reads:
                break
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

def find_primer_boundary(states, min_len=10, max_len=30):
    """
    Primer = longest [CD]* prefix from position 0, terminated at first V.

    Constraints:
      - Minimum length 10 bp (V before pos 10 → no primer).
      - Maximum length 30 bp.
      - Noise tolerance: a *single* isolated V is skipped if it has
        >= 3 non-V positions on each side (likely sequencing noise).

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
            # Check if this is an isolated V (noise)
            is_single = (i + 1 < len(states) and states[i + 1] != 'V')
            before_ok = (i >= 3 and
                         all(s != 'V' for s in states[max(0, i - 3):i]))
            after_ok = (i + 3 < len(states) and
                        all(s != 'V' for s in states[i + 1:i + 4]))
            if is_single and before_ok and after_ok:
                # Isolated V — treat as noise, skip
                end = i + 1
                i += 1
            else:
                # Real variable region — primer ends here
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

def detect_for_reads(reads, label):
    """
    Run three-state (C/D/V) primer detection on pre-filtered reads.

    Steps:
      2. Build 60×4 position-wise base frequency matrix.
      3. Classify each position as C (Conserved), D (Degenerate), V (Variable).
      4. Find primer boundary: longest [CD]* prefix with noise tolerance.
      5. Build consensus (C→dominant base, D/V→N) and report.
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

    # Step 4: Find primer boundary
    print(f"\n[Step 4] Finding primer boundary "
          f"(min=10bp, max=30bp) ...", file=sys.stderr)
    primer_length = find_primer_boundary(states, min_len=10, max_len=30)

    # Step 5: Build consensus and report
    print(f"\n[Step 5] Result:", file=sys.stderr)
    if primer_length == 0:
        result = dict(detected=False, primer_length=0, consensus="",
                      message="No primer detected")
    else:
        consensus = build_consensus_cdv(freq_matrix, states, primer_length)
        result = dict(detected=True, primer_length=primer_length,
                      consensus=consensus,
                      message=f"Primer detected (length={primer_length}bp)")

    print(f"  1. Primer exists:     "
          f"{'Yes' if result['detected'] else 'No'}", file=sys.stderr)
    print(f"  2. Truncation length: {result['primer_length']} bp",
          file=sys.stderr)
    print(f"  3. Consensus:         "
          f"{result['consensus'] if result['consensus'] else 'N/A'}",
          file=sys.stderr)

    return result


def detect_for_file(filepath, label):
    """Run Steps 1-6 on a single FASTQ file. Returns result dict."""
    print(f"\n[Step 1] Reading and filtering reads from "
          f"{os.path.basename(filepath)} ...", file=sys.stderr)
    reads = read_and_filter(filepath)
    if not reads:
        print("  ERROR: No reads passed filters", file=sys.stderr)
        return dict(detected=False, primer_length=0, consensus="",
                    message="No primer detected (no reads)")
    return detect_for_reads(reads, label)


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
                                     "R1 majority (forward primer)")
        r2_result = detect_for_reads(mixed_info['minority_reads'],
                                     "R1 minority (reverse primer)")

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
                                     "R1" if mode == "PE" else "SE")

        r2_result = None
        if mode == "PE" and second_file:
            r2_result = detect_for_file(second_file, "R2")

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

    print(f"\n  Output directory: {output_dir}", file=sys.stderr)
    print("  Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
