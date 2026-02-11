#!/usr/bin/env python3
"""
Entropy-based 16S Primer Detection and Removal
(with mixed R1/R2 orientation support)

Principle:
  Primer positions are conserved across reads → low Shannon entropy.
  Biological positions are variable across reads → high Shannon entropy.
  The entropy jump point marks the primer→biology transition.

Algorithm:
  Step 0:   Scan directory, pick first sample (PE or SE).
  Step 1:   Read & pre-filter 2000 R1 reads (length, quality, complexity).
  Step 2:   Build 60×4 position-wise base frequency matrix.
  Phase 1:  Mixed-orientation detection (PE only).
            Check first 15 positions for bimodal base distributions.
            If >50% bimodal → mixed orientation detected.

  --- Branch A (mixed orientation) ---
  Phase 2:  Split sampled reads into majority/minority groups by the base
            at the best bimodal position. Run entropy detection (Steps 3-6)
            on each group independently to determine forward/reverse primer
            lengths.
  Phase 3:  Trim ALL files with per-read orientation correction.
            For each read pair, classify by a single-character check at
            best_pos. Flipped reads have R1↔R2 swapped before trimming.

  --- Branch B (single orientation) ---
  Step 3:   Compute Shannon entropy per position.
  Step 4:   Smooth with sliding window (w=3).
  Step 5:   Detect jump point (primer→biology transition).
  Step 6:   Classify result and extract consensus.
  Step 7:   Trim ALL files using detected primer length.
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


def read_and_filter(filepath, max_reads=2000, min_len=50,
                    min_avg_qual=20, min_complexity=0.3):
    """
    Step 1: Read FASTQ, apply quality filters, return up to max_reads reads.
    Returns list of (seq_string, qual_string) tuples.
    """
    reads = []
    total = 0
    disc_len = disc_qual = disc_cplx = 0

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

            reads.append((seq, qual))

    print(f"  Scanned {total} reads, kept {len(reads)}", file=sys.stderr)
    print(f"  Discarded: {disc_len} (length<{min_len}), "
          f"{disc_qual} (avgQ<{min_avg_qual}), "
          f"{disc_cplx} (complexity<{min_complexity})", file=sys.stderr)

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
# Step 3: Shannon entropy
# ===========================================================================

def compute_entropy(freq_matrix):
    """H(i) = -sum( p * log2(p) ) for each position."""
    entropy = []
    for freqs in freq_matrix:
        h = 0.0
        for p in freqs:
            if p > 0:
                h -= p * math.log2(p)
        entropy.append(h)
    return entropy


# ===========================================================================
# Step 4: Smoothing
# ===========================================================================

def smooth_entropy(entropy, window=3):
    """Sliding-window mean smoothing."""
    n = len(entropy)
    half = window // 2
    smoothed = []
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        smoothed.append(sum(entropy[lo:hi]) / (hi - lo))
    return smoothed


# ===========================================================================
# Step 5: Jump-point detection
# ===========================================================================

def detect_jump_point(smoothed):
    """
    baseline = mean(smoothed[40:60])
    Scan from position 0 rightward; first i satisfying BOTH:
      A) smoothed[i] > baseline * 0.75
      C) mean(smoothed[i:i+5]) > baseline * 0.7
    primer_length = i
    """
    baseline = sum(smoothed[40:60]) / 20.0

    if baseline < 0.1:
        print(f"  WARNING: Very low baseline entropy ({baseline:.4f})",
              file=sys.stderr)
        return 0, baseline

    for i in range(len(smoothed) - 5):
        cond_a = smoothed[i] > baseline * 0.75
        win_mean = sum(smoothed[i:i + 5]) / 5.0
        cond_c = win_mean > baseline * 0.7
        if cond_a and cond_c:
            return i, baseline

    # No jump found within analysis range
    return 0, baseline


# ===========================================================================
# Step 6: Classification & consensus
# ===========================================================================

IUPAC_MAP = {
    frozenset('A'): 'A', frozenset('C'): 'C',
    frozenset('G'): 'G', frozenset('T'): 'T',
    frozenset('AG'): 'R', frozenset('CT'): 'Y',
    frozenset('CG'): 'S', frozenset('AT'): 'W',
    frozenset('GT'): 'K', frozenset('AC'): 'M',
    frozenset('CGT'): 'B', frozenset('AGT'): 'D',
    frozenset('ACT'): 'H', frozenset('ACG'): 'V',
    frozenset('ACGT'): 'N',
}

BASES = ['A', 'C', 'G', 'T']


def _iupac_code(freqs):
    """IUPAC ambiguity code for one position's frequency vector."""
    max_freq = max(freqs)
    if max_freq >= 0.7:
        return BASES[freqs.index(max_freq)]
    present = frozenset(b for b, f in zip(BASES, freqs) if f > 0.1)
    if not present:
        return 'N'
    return IUPAC_MAP.get(present, 'N')


def classify_result(jump_point, baseline, freq_matrix):
    """Return a result dict describing primer status."""
    consensus = "".join(_iupac_code(freq_matrix[i])
                        for i in range(jump_point))

    if jump_point <= 1:
        return dict(detected=False, primer_length=0, confidence=None,
                    consensus="", message="No primer detected")
    elif 10 <= jump_point <= 35:
        return dict(detected=True, primer_length=jump_point,
                    confidence="high", consensus=consensus,
                    message=f"Primer detected (length={jump_point}bp, "
                            f"confidence=high)")
    elif jump_point > 35:
        return dict(detected=True, primer_length=jump_point,
                    confidence="low", consensus=consensus,
                    message=f"Adapter+primer detected "
                            f"(length={jump_point}bp, confidence=low)")
    else:  # 2-9
        return dict(detected=True, primer_length=jump_point,
                    confidence="low", consensus=consensus,
                    message=f"Short technical sequence "
                            f"(length={jump_point}bp, confidence=low)")


# ===========================================================================
# Step 7: Trimming helpers
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
    """Run Steps 2-6 on a pre-filtered list of (seq, qual) tuples."""
    print(f"\n{'=' * 60}", file=sys.stderr)
    print(f"  Analyzing {label}: {len(reads)} reads", file=sys.stderr)
    print(f"{'=' * 60}", file=sys.stderr)

    if not reads:
        print("  ERROR: No reads to analyze", file=sys.stderr)
        return classify_result(0, 0.0, [[0.25]*4]*60)

    # Step 2
    print(f"\n[Step 2] Building frequency matrix (60 positions) ...",
          file=sys.stderr)
    freq_matrix = build_frequency_matrix(reads, num_positions=60)

    # Step 3
    print(f"\n[Step 3] Computing Shannon entropy ...", file=sys.stderr)
    entropy = compute_entropy(freq_matrix)
    print(f"  Entropy[0:5]   = {['%.3f' % e for e in entropy[:5]]}",
          file=sys.stderr)
    print(f"  Entropy[25:30] = {['%.3f' % e for e in entropy[25:30]]}",
          file=sys.stderr)
    print(f"  Entropy[55:60] = {['%.3f' % e for e in entropy[55:60]]}",
          file=sys.stderr)

    # Step 4
    print(f"\n[Step 4] Smoothing (window=3) ...", file=sys.stderr)
    smoothed = smooth_entropy(entropy, window=3)

    # Step 5
    print(f"\n[Step 5] Detecting jump point ...", file=sys.stderr)
    jump_point, baseline = detect_jump_point(smoothed)
    print(f"  Baseline entropy (pos 40-60): {baseline:.4f}", file=sys.stderr)
    print(f"  Jump point: position {jump_point}", file=sys.stderr)

    # Step 6
    print(f"\n[Step 6] Classification ...", file=sys.stderr)
    result = classify_result(jump_point, baseline, freq_matrix)
    print(f"  {result['message']}", file=sys.stderr)
    if result['consensus']:
        print(f"  Consensus: {result['consensus']}", file=sys.stderr)

    return result


def detect_for_file(filepath, label):
    """Run Steps 1-6 on a single FASTQ file. Returns result dict."""
    print(f"\n[Step 1] Reading and filtering reads from "
          f"{os.path.basename(filepath)} ...", file=sys.stderr)
    reads = read_and_filter(filepath)
    if not reads:
        print("  ERROR: No reads passed filters", file=sys.stderr)
        return classify_result(0, 0.0, [[0.25]*4]*60)
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
        description="Entropy-based 16S primer detection and removal")
    parser.add_argument("-i", "--input", required=True,
                        help="Input directory containing FASTQ files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for trimmed files")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60, file=sys.stderr)
    print("ENTROPY-BASED PRIMER DETECTION", file=sys.stderr)
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
    # Step 2: Build frequency matrix (used for both entropy and
    #         mixed-orientation detection)
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
        # Phase 2: Entropy detection on each orientation group
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

        # Step 7: Trim the ENTIRE dataset
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
