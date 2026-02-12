#!/usr/bin/env python3
"""
Unused functions removed from py_16s.py (2026-02-12).

These functions were part of older pipeline approaches that have been
superseded by the CDV three-state primer detection (entropy_primer_detect.py)
and the current DADA2-based workflow.

Contents:
  - combine_stats_file:     Deblur stats merging (deblur no longer used)
  - add_prefix_to_file:     FASTA/BIOM prefix tool (never called in pipeline)
  - trim_pos_454:           454 platform trim calc (454 removed from pipeline)
  - trim_pos_dada2:         Placeholder, always returned None
  - detect_primers_16s:     BLAST-based primer detection (replaced by CDV)
  - _check_primer_blast:    Helper for detect_primers_16s
  - smart_trim_16s:         BBDuk/cutadapt primer detection (replaced by CDV)
"""

import sys
import subprocess
import os
import glob
import gzip
import statistics
from collections import Counter
import pandas as pd
from Bio import SeqIO


def combine_stats_file(filter_stats, deblur_stats, dataset_name):
    """Combine filter and deblur statistics"""
    filter_stats_df = pd.read_csv(filter_stats)
    deblur_stats_df = pd.read_csv(deblur_stats)
    merged_df = pd.merge(filter_stats_df, deblur_stats_df, on='sample-id', how='outer')
    out_path = os.path.join(os.path.dirname(filter_stats), f"{dataset_name}_process_stats.csv")
    merged_df.to_csv(out_path, index=False)


def add_prefix_to_file(in_fasta, in_table, prefix):
    """Add prefix to FASTA and BIOM files"""
    from Bio import SeqIO
    from biom import load_table, Table
    import h5py

    dir = os.path.dirname(in_fasta)
    out_fasta = os.path.join(dir, f"{prefix}-dna-sequence.fasta")

    modified_sequences = []
    for record in SeqIO.parse(in_fasta, "fasta"):
        record.id = prefix + record.id
        record.description = prefix + record.description
        modified_sequences.append(record)

    SeqIO.write(modified_sequences, out_fasta, "fasta")

    table = load_table(in_table)
    out_table = os.path.join(dir, f"{prefix}-feature-table.biom")
    observation_ids = table.ids(axis='observation')
    new_observation_ids = [prefix + obs_id for obs_id in observation_ids]
    new_data = table.matrix_data.tocsr()
    new_table = Table(new_data, new_observation_ids, table.ids())
    with h5py.File(out_table, 'w') as f:
        new_table.to_hdf5(f, 'representive sequence')


def trim_pos_454(file_path):
    """Calculate trim length for 454 sequencing"""
    with open(file_path, "r") as f:
        tsv = [line.rstrip("\n") for line in f]

    pos = [int(x) for x in tsv[0].split("\t")[1:]]
    counts = [float(x) for x in tsv[1].split("\t")[1:]]
    L = len(pos)

    threshold = 9000
    under_threshold = [i for i, v in enumerate(counts) if v < threshold]

    trim_length = under_threshold[0] if under_threshold else L

    if trim_length == 0:
        return (None,)

    print(f"{trim_length}")
    return (trim_length,)


def trim_pos_dada2(file_path):
    """Placeholder for DADA2 trim position calculation"""
    return (None, None)


def detect_primers_16s(input_path, tmp_path="/tmp", ref_path="./Meta2Data/docs/J01859.1.fna"):
    """
    Detect 16S primers using alignment-based approach

    Method:
    1. Sample 1000 reads from FASTQ file
    2. BLAST align reads to E. coli 16S reference (only plus strand)
    3. Find consensus alignment start position
    4. Calculate distance to nearest LEGAL_ANCHOR
    5. Determine if primers need trimming:
       - Distance < 10bp: Primers already removed
       - Distance >= 18bp: Primers present, calculate trim length

    Returns: (primers_found: bool, trim_length: int)
    """
    print("=" * 60, file=sys.stderr)
    print("16S PRIMER DETECTION (ALIGNMENT-BASED)", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print("\nStrategy: Align reads to reference -> Find consensus position -> Detect primers\n", file=sys.stderr)

    if not os.path.exists(ref_path):
        print(f"Reference not found: {ref_path}", file=sys.stderr)
        return False, 0

    END_ANCHORS = [28, 64, 357, 533, 799]
    print(f"Primer end positions (16S starts after): {END_ANCHORS}", file=sys.stderr)

    if not os.path.exists(f"{ref_path}.nhr"):
        subprocess.run(
            f"makeblastdb -in {ref_path} -dbtype nucl",
            shell=True, check=True, capture_output=True
        )

    fwd_files = sorted(glob.glob(os.path.join(input_path, "*_R1*.fastq*")))
    if not fwd_files:
        fwd_files = sorted(glob.glob(os.path.join(input_path, "*.fastq*")))
        if not fwd_files:
            return False, 0

    test_file = fwd_files[0]

    sampled_reads = []
    with (gzip.open(test_file, "rt") if test_file.endswith(".gz") else open(test_file, "r")) as h:
        for i, rec in enumerate(SeqIO.parse(h, "fastq")):
            if i >= 1000:
                break
            sampled_reads.append((f"read{i}", str(rec.seq)))

    if len(sampled_reads) < 100:
        return False, 0

    os.makedirs(tmp_path, exist_ok=True)
    query_fa = os.path.join(tmp_path, "primer_detect_reads.fa")
    with open(query_fa, 'w') as f:
        for read_id, seq in sampled_reads:
            f.write(f">{read_id}\n{seq}\n")

    blast_out = os.path.join(tmp_path, "primer_detect_alignment.txt")
    blast_cmd = (f"blastn -query {query_fa} -db {ref_path} "
                 f"-outfmt '6 qseqid qstart qend sstart send sstrand qlen' "
                 f"-out {blast_out} -task blastn-short -max_target_seqs 1 -evalue 1e-5")

    try:
        subprocess.run(blast_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        return False, 0

    alignment_starts = []
    query_starts = []

    with open(blast_out) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                qseqid, qstart, qend, sstart, send, sstrand, qlen = parts
                if sstrand == 'plus':
                    alignment_starts.append(int(sstart))
                    query_starts.append(int(qstart) - 1)

    if len(alignment_starts) < 50:
        return False, 0

    start_counter = Counter(alignment_starts)
    consensus_start, _ = start_counter.most_common(1)[0]
    nearest_anchor = min(END_ANCHORS, key=lambda a: abs(a - consensus_start))
    ref_distance = nearest_anchor - consensus_start

    if ref_distance < 10:
        return True, 0

    trim_candidates = []
    for align_start, query_start in zip(alignment_starts, query_starts):
        if abs(align_start - consensus_start) <= 5:
            trim_candidates.append(query_start + ref_distance)

    if not trim_candidates:
        trim_candidates = [qs + ref_distance for qs in query_starts]

    return True, int(statistics.median(trim_candidates))


def _check_primer_blast(seq, offset, ref_path, tmp_path, legal_anchors):
    """
    Helper function to BLAST a candidate primer sequence.
    Returns trim_length if valid primer found, None otherwise.
    """
    tmp_fa = os.path.join(tmp_path, f"primer_detect_offset{offset}.fa")
    with open(tmp_fa, "w") as f:
        f.write(f">candidate_offset{offset}\n{seq}")

    try:
        res = subprocess.check_output(
            f"blastn -query {tmp_fa} -db {ref_path} -outfmt '6 sstart' -task blastn-short",
            shell=True, stderr=subprocess.DEVNULL
        ).decode().strip()

        if res:
            sstart = int(res.split('\n')[0])
            for anchor in legal_anchors:
                if abs(sstart - anchor) < 20:
                    return offset + 20
    except subprocess.CalledProcessError:
        pass

    return None


def smart_trim_16s(input_path, output_path, tmp_path="/tmp", ref_path="./Meta2Data/docs/J01859.1.fna"):
    """Detect and trim 16S primers using BBDuk frequency + BLAST approach"""
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if not os.path.exists(ref_path):
        sys.exit(1)

    LEGAL_ANCHORS = [8, 27, 338, 341, 515, 518, 534, 785, 799, 806, 907, 926, 1046, 1099, 1100, 1391, 1492]

    if not os.path.exists(f"{ref_path}.nhr"):
        subprocess.run(
            f"makeblastdb -in {ref_path} -dbtype nucl",
            shell=True, check=True, capture_output=True
        )

    fwd_files = sorted(glob.glob(os.path.join(input_path, "*_R1*.fastq*")))
    if not fwd_files:
        sys.exit(1)

    test_file = fwd_files[0]

    prefixes = Counter()
    with (gzip.open(test_file, "rt") if test_file.endswith(".gz") else open(test_file, "r")) as h:
        for i, rec in enumerate(SeqIO.parse(h, "fastq")):
            if i >= 1000:
                break
            prefixes[str(rec.seq[:20])] += 1

    if not prefixes:
        sys.exit(1)

    top_seq, count = prefixes.most_common(1)[0]
    has_primers = False

    if count > 500:
        tmp_fa = os.path.join(tmp_path, "top_candidate.fa")
        with open(tmp_fa, "w") as f:
            f.write(f">top_candidate\n{top_seq}")

        try:
            res = subprocess.check_output(
                f"blastn -query {tmp_fa} -db {ref_path} -outfmt '6 sstart' -task blastn-short",
                shell=True, stderr=subprocess.DEVNULL
            ).decode().strip()

            if res:
                sstart = int(res.split('\n')[0])
                has_primers = any(abs(sstart - anchor) < 20 for anchor in LEGAL_ANCHORS)
        except subprocess.CalledProcessError:
            pass

    if has_primers:
        for r1 in fwd_files:
            r2 = r1.replace("_R1", "_R2")
            if not os.path.exists(r2):
                continue
            out1 = os.path.join(output_path, os.path.basename(r1))
            out2 = os.path.join(output_path, os.path.basename(r2))
            subprocess.run(
                f"cutadapt -u 20 -U 20 -o {out1} -p {out2} {r1} {r2} --quiet",
                shell=True, check=True
            )
    else:
        for f in glob.glob(os.path.join(input_path, "*.fastq*")):
            target = os.path.join(output_path, os.path.basename(f))
            if not os.path.exists(target):
                os.symlink(os.path.abspath(f), target)
