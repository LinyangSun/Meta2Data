#!/usr/bin/env python3

import sys
import subprocess
import os
import glob
import gzip
from collections import Counter
import pandas as pd
import numpy as np
from Bio import SeqIO


def check_and_install(module, module2):
    """Auto-install missing modules"""
    try:
        __import__(module)
    except ImportError:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", module2],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        __import__(module)


def GenerateDatasetsIDsFile(file_path, Bioproject, Data_SequencingPlatform=None, output_dir=None):
    """
    Generate datasets ID file from CSV

    Extracts unique BioProject IDs from metadata CSV.
    Platform parameter is optional since platform is detected dynamically downstream.

    Args:
        file_path: Path to metadata CSV
        Bioproject: Column name for BioProject ID
        Data_SequencingPlatform: (Optional) Column name for platform - if provided, outputs both columns
        output_dir: (Optional) Directory to write output file. If not provided, uses input file directory

    Returns:
        Array of dataset IDs
    """
    # Use output directory if provided, otherwise use input file's directory
    if output_dir:
        directory_path = output_dir
    else:
        directory_path = os.path.dirname(os.path.abspath(file_path))

    df = pd.read_csv(file_path)

    # Clean BioProject column
    if Bioproject in df.columns:
        df[Bioproject] = (
            df[Bioproject]
            .astype(str)
            .str.replace(' ', '', regex=True)
            .str.replace('\n', '', regex=True)
            .str.replace('\t', '', regex=True)
        )

    # If platform column is provided and exists, include it (backward compatibility)
    if Data_SequencingPlatform and Data_SequencingPlatform in df.columns:
        df[Data_SequencingPlatform] = (
            df[Data_SequencingPlatform]
            .astype(str)
            .str.replace(' ', '', regex=True)
            .str.replace('\n', '', regex=True)
            .str.replace('\t', '', regex=True)
        )
        # Output both columns
        df_pair = (
            df[[Bioproject, Data_SequencingPlatform]]
            .dropna(subset=[Bioproject])
            .drop_duplicates()
        )
    else:
        # Output only BioProject column (default behavior)
        # Platform is detected dynamically later via get_sequencing_platform()
        df_pair = (
            df[[Bioproject]]
            .dropna(subset=[Bioproject])
            .drop_duplicates()
        )

    datasets = df_pair.values.astype(object)
    out_path = f"{directory_path}/datasets_ID.txt"
    np.savetxt(out_path, datasets, fmt="%s")
    return datasets


def GenerateSRAsFile(file_path, Bioproject, SRA_Number, Biosample=None, output_dir=None):
    """
    Generate SRA files for each bioproject

    Args:
        file_path: Path to metadata CSV
        Bioproject: Column name for BioProject ID
        SRA_Number: Column name for SRA/Run accession
        Biosample: (Optional) Column name for BioSample ID - if not provided, uses SRA_Number for naming
        output_dir: (Optional) Directory to write output files
    """
    # Use output directory if provided, otherwise use input file's directory
    if output_dir:
        directory_path = output_dir
    else:
        directory_path = os.path.dirname(os.path.abspath(file_path))
    df = pd.read_csv(file_path)

    # Determine which columns to clean
    columns_to_clean = [Bioproject, SRA_Number]
    if Biosample and Biosample in df.columns:
        columns_to_clean.append(Biosample)

    for col in columns_to_clean:
        if col in df.columns:
            df[col] = (
                df[col]
                .astype(str)
                .str.replace(' ', '', regex=True)
                .str.replace('\n', '', regex=True)
                .str.replace('\t', '', regex=True)
            )

    # Create rename column: use Biosample if available, otherwise use SRA_Number
    if Biosample and Biosample in df.columns:
        df["rename"] = df[Bioproject] + '-' + df[Biosample]
    else:
        # Fallback: use SRA_Number (Run ID) when Biosample is not available
        df["rename"] = df[Bioproject] + '-' + df[SRA_Number]

    datasets = np.array(
        [str(x).strip().replace('\t', '') for x in df[Bioproject].dropna().unique()],
        dtype=object
    )

    for value in datasets:
        sub = df[df[Bioproject] == value].copy()
        out_df = pd.concat([
            sub[SRA_Number],
            sub["rename"]
        ], axis=1)

        folder_path = f"{directory_path}/{value}/"
        os.makedirs(folder_path, exist_ok=True)
        out_file = f'{folder_path}/{value}_sra.txt'
        out_df.to_csv(out_file, sep='\t', header=None, index=None)


def mk_manifest_SE(file_path):
    """Generate single-end manifest file"""
    df = pd.read_csv(file_path, sep='\t', header=None)
    dataset_name = os.path.basename(file_path).replace("-file.txt", "")
    df1 = pd.DataFrame({'sample-id': [], 'absolute-filepath': []})

    for i, row in df.iterrows():
        basenames = os.path.basename(row[0])
        path = os.path.dirname(row[0])
        sample = '.'.join(basenames.split('.')[:-1])
        df1.loc[i, 'sample-id'] = sample
        df1.loc[i, 'absolute-filepath'] = os.path.join(path, basenames)

    df_unique = df1.drop_duplicates(subset=['sample-id'])
    out_path = os.path.join(os.path.dirname(file_path), f"{dataset_name}_manifest.tsv")
    df_unique.to_csv(out_path, sep='\t', index=False)


def mk_manifest_PE(file_path):
    """Generate paired-end manifest file"""
    df = pd.read_csv(file_path, sep='\t', header=None)
    dataset_name = os.path.basename(file_path).replace("-file.txt", "")
    df1 = pd.DataFrame({
        'sample-id': [],
        'forward-absolute-filepath': [],
        'reverse-absolute-filepath': []
    })

    for i, row in df.iterrows():
        basenames = os.path.basename(row[0])
        sample = basenames.rsplit('_', 1)[0]
        path = os.path.dirname(row[0])
        df1.loc[i, 'sample-id'] = sample
        df1.loc[i, 'forward-absolute-filepath'] = os.path.join(path, f"{sample}_1.fastq")
        df1.loc[i, 'reverse-absolute-filepath'] = os.path.join(path, f"{sample}_2.fastq")

    df_unique = df1.drop_duplicates(subset=['sample-id'])
    out_path = os.path.join(os.path.dirname(file_path), f"{dataset_name}_manifest.tsv")
    df_unique.to_csv(out_path, sep='\t', index=False)


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


def trim_pos_deblur(file_path):
    """Calculate trim positions for Deblur"""
    with open(file_path, "r") as f:
        tsv = [line.rstrip("\n") for line in f]
    
    pos = [int(x) for x in tsv[0].split("\t")[1:]]
    L = len(pos)
    counts = [float(x) for x in tsv[1].split("\t")[1:]]
    
    threshold = min(9000, max(counts) * 0.7)
    under = [i for i, v in enumerate(counts) if v < threshold]
    outer1 = max(under[0] if under else (L - 1), 40)
    
    inner_region = range(0, 20)
    outer_region = range(max(0, outer1 - 20), outer1)
    
    rows = []
    for line in tsv[3:]:
        vals = [float(x) for x in line.split("\t")[1:][:L]]
        rows.append(vals)
    
    per_row = []
    for vals in rows:
        na_idx = {i for i, v in enumerate(vals) if v < 20}
        
        # Calculate start position
        inner_na = sorted(i for i in na_idx if i in inner_region)
        n_inner = len(inner_na)
        
        if n_inner == 0:
            start = 0
        elif n_inner >= 15:
            start = None
        else:
            s = set(inner_na)
            r = 0
            while r in s and r < 20:
                r += 1
            start = r if r > 0 else 0
        
        # Calculate end position
        outer_na = sorted(i for i in na_idx if i in outer_region)
        n_outer = len(outer_na)
        
        if n_outer == 0:
            end = outer1
        elif n_outer >= 15:
            end = None
        else:
            s = set(outer_na)
            r = 0
            while (outer1 - r) in s:
                r += 1
            end = outer1 - r - 1 if r > 0 else outer1
        
        per_row.append((start, end))
    
    starts = [s for s, e in per_row if s is not None]
    ends = [e for s, e in per_row if e is not None]
    
    if not starts or not ends:
        return (None, None)
    
    valid_rows = sum(1 for s, e in per_row if s is not None and e is not None)
    if valid_rows < len(rows) * 0.5:
        return (None, None)
    
    final_start = max(0, min(max(starts), L - 1)) + 1
    final_end = max(0, min(min(ends), L - 1)) - 1
    
    if final_start >= final_end:
        return (None, None)
    
    print(f"{final_start},{final_end}")
    return (final_start, final_end)


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
    Detect 16S primers in FASTQ files using sliding window approach

    Method:
    1. Sample first 1000 reads from FASTQ file
    2. Extract 20bp sequences at different offsets (positions)
    3. Find most common 20bp sequence at each offset (conserved = likely primer)
    4. BLAST against E. coli 16S reference to verify it's a real primer
    5. Check if BLAST hit matches known 16S primer binding sites

    Returns: (primers_found: bool, trim_length: int)
    """
    print("=" * 60, file=sys.stderr)
    print("16S PRIMER DETECTION METHOD", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print("\nPurpose: Detect if reads contain adapters/barcodes + 16S primers", file=sys.stderr)
    print("Strategy: Two-stage sliding window search with BLAST verification\n", file=sys.stderr)

    if not os.path.exists(ref_path):
        print(f"✗ Reference not found: {ref_path}", file=sys.stderr)
        return False, 0

    # Known 16S primer binding sites on E. coli 16S rRNA gene
    LEGAL_ANCHORS = [8, 27, 338, 341, 515, 518, 534, 785, 799, 806, 907, 926, 1046, 1099, 1100, 1391, 1492]
    print(f"Known 16S primer sites: {len(LEGAL_ANCHORS)} positions", file=sys.stderr)
    print(f"  (e.g., V3-V4: 341F/785R, V4: 515F/806R)\n", file=sys.stderr)

    # Two-stage search strategy
    COARSE_OFFSETS = [0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60]  # Fast initial search
    FINE_SEARCH_RANGE = 60  # If coarse fails, try every position

    print(f"Stage 1: Coarse search at {len(COARSE_OFFSETS)} key offsets", file=sys.stderr)
    print(f"  Offsets: {COARSE_OFFSETS}", file=sys.stderr)
    print(f"Stage 2: Fine search (0-{FINE_SEARCH_RANGE}bp) if Stage 1 fails\n", file=sys.stderr)

    # Build BLAST database if needed
    if not os.path.exists(f"{ref_path}.nhr"):
        print("Building BLAST database for 16S reference...", file=sys.stderr)
        subprocess.run(
            f"makeblastdb -in {ref_path} -dbtype nucl",
            shell=True, check=True, capture_output=True
        )
        print("✓ BLAST database ready\n", file=sys.stderr)

    # Find input files
    fwd_files = sorted(glob.glob(os.path.join(input_path, "*_R1*.fastq*")))
    if not fwd_files:
        # Try single-end pattern
        fwd_files = sorted(glob.glob(os.path.join(input_path, "*.fastq*")))
        if not fwd_files:
            print("✗ No FASTQ files found", file=sys.stderr)
            return False, 0

    test_file = fwd_files[0]
    print(f"Analyzing: {os.path.basename(test_file)}", file=sys.stderr)
    print(f"Sampling: First 1000 reads\n", file=sys.stderr)

    # Stage 1: Coarse search with fixed offsets
    print("-" * 60, file=sys.stderr)
    print("STAGE 1: Coarse Search (Fast)", file=sys.stderr)
    print("-" * 60, file=sys.stderr)
    offset_prefixes = {offset: Counter() for offset in COARSE_OFFSETS}

    with (gzip.open(test_file, "rt") if test_file.endswith(".gz") else open(test_file, "r")) as h:
        for i, rec in enumerate(SeqIO.parse(h, "fastq")):
            if i >= 1000:
                break
            for offset in COARSE_OFFSETS:
                if len(rec.seq) >= offset + 20:
                    offset_prefixes[offset][str(rec.seq[offset:offset+20])] += 1

    print("Extracting most common 20bp sequences at each offset...", file=sys.stderr)

    # Check each coarse offset for primer matches
    os.makedirs(tmp_path, exist_ok=True)

    for offset in COARSE_OFFSETS:
        if not offset_prefixes[offset]:
            continue

        top_seq, count = offset_prefixes[offset].most_common(1)[0]
        frequency_pct = (count / 1000) * 100

        if count > 500:  # High frequency suggests conserved primer sequence
            print(f"\nOffset {offset:2d}bp: Found conserved sequence (freq: {count}/1000 = {frequency_pct:.1f}%)", file=sys.stderr)
            print(f"  Sequence: {top_seq}", file=sys.stderr)
            print(f"  BLASTing against 16S reference...", file=sys.stderr)

            result = _check_primer_blast(top_seq, offset, ref_path, tmp_path, LEGAL_ANCHORS)
            if result is not None:
                print(f"  ✓ CONFIRMED: Matches known 16S primer site!", file=sys.stderr)
                print(f"  → Trim length: {result}bp (adapter/barcode={offset}bp + primer=20bp)", file=sys.stderr)
                print("=" * 60, file=sys.stderr)
                return True, result
            else:
                print(f"  ✗ No match to 16S primer sites", file=sys.stderr)

    print(f"\n✗ Stage 1: No primers detected at coarse offsets", file=sys.stderr)

    # Stage 2: Fine search if coarse search failed
    print("\n" + "-" * 60, file=sys.stderr)
    print("STAGE 2: Fine Search (Thorough)", file=sys.stderr)
    print("-" * 60, file=sys.stderr)
    print(f"Searching every position from 0-{FINE_SEARCH_RANGE}bp...\n", file=sys.stderr)

    fine_offset_prefixes = {offset: Counter() for offset in range(0, FINE_SEARCH_RANGE)}
    
    with (gzip.open(test_file, "rt") if test_file.endswith(".gz") else open(test_file, "r")) as h:
        for i, rec in enumerate(SeqIO.parse(h, "fastq")):
            if i >= 1000:
                break
            for offset in range(0, FINE_SEARCH_RANGE):
                if len(rec.seq) >= offset + 20:
                    fine_offset_prefixes[offset][str(rec.seq[offset:offset+20])] += 1

    # Check fine offsets
    candidates_found = 0
    for offset in range(0, FINE_SEARCH_RANGE):
        if not fine_offset_prefixes[offset]:
            continue

        top_seq, count = fine_offset_prefixes[offset].most_common(1)[0]
        frequency_pct = (count / 1000) * 100

        if count > 500:
            candidates_found += 1
            print(f"Offset {offset:2d}bp: Conserved sequence (freq: {count}/1000 = {frequency_pct:.1f}%)", file=sys.stderr)
            result = _check_primer_blast(top_seq, offset, ref_path, tmp_path, LEGAL_ANCHORS)
            if result is not None:
                print(f"  ✓ CONFIRMED: Matches known 16S primer site!", file=sys.stderr)
                print(f"  → Trim length: {result}bp", file=sys.stderr)
                print("=" * 60, file=sys.stderr)
                return True, result

    print(f"\n✗ Stage 2: Checked {FINE_SEARCH_RANGE} positions, found {candidates_found} candidates", file=sys.stderr)
    print("  None matched known 16S primer sites\n", file=sys.stderr)

    print("=" * 60, file=sys.stderr)
    print("CONCLUSION: No 16S primers detected", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print("\nInterpretation:", file=sys.stderr)
    print("  • Reads may already be primer-trimmed", file=sys.stderr)
    print("  • Or using non-standard primer sequences", file=sys.stderr)
    print("  • Or primer sequences are variable (not conserved)", file=sys.stderr)
    print("  → Pipeline will proceed without primer trimming\n", file=sys.stderr)

    return False, 0


def _check_primer_blast(seq, offset, ref_path, tmp_path, legal_anchors):
    """
    Helper function to BLAST a candidate primer sequence

    Process:
    1. Write candidate sequence to FASTA file
    2. BLAST against E. coli 16S reference
    3. Get alignment start position on reference
    4. Check if position matches known primer binding sites (±20bp tolerance)

    Returns trim_length if valid primer found, None otherwise
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
            print(f"    BLAST hit at position: {sstart}bp on 16S gene", file=sys.stderr)

            # Check if within 20bp of any known primer site
            for anchor in legal_anchors:
                distance = abs(sstart - anchor)
                if distance < 20:
                    print(f"    Match! Within {distance}bp of primer site {anchor}", file=sys.stderr)
                    # Trim length = offset (adapter/barcode) + 20 (primer)
                    trim_length = offset + 20
                    return trim_length

            print(f"    No match to known primer sites (closest: {min(abs(sstart - a) for a in legal_anchors)}bp away)", file=sys.stderr)
        else:
            print(f"    No BLAST hits found", file=sys.stderr)

    except subprocess.CalledProcessError as e:
        print(f"    BLAST error: {e}", file=sys.stderr)

    return None


def smart_trim_16s(input_path, output_path, tmp_path="/tmp", ref_path="./Meta2Data/docs/J01859.1.fna"):
    """Detect and trim 16S primers"""
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if not os.path.exists(ref_path):
        sys.exit(1)

    LEGAL_ANCHORS = [8, 27, 338, 341, 515, 518, 534, 785, 799, 806, 907, 926, 1046, 1099, 1100, 1391, 1492]

    # Build BLAST database if needed
    if not os.path.exists(f"{ref_path}.nhr"):
        subprocess.run(
            f"makeblastdb -in {ref_path} -dbtype nucl",
            shell=True, check=True, capture_output=True
        )

    fwd_files = sorted(glob.glob(os.path.join(input_path, "*_R1*.fastq*")))
    if not fwd_files:
        sys.exit(1)

    test_file = fwd_files[0]

    # Count prefix frequencies
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
        # Trim primers using cutadapt
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
        # Symlink files without trimming
        for f in glob.glob(os.path.join(input_path, "*.fastq*")):
            target = os.path.join(output_path, os.path.basename(f))
            if not os.path.exists(target):
                os.symlink(os.path.abspath(f), target)


def get_sequencing_platform(srr_id, bioproject_id=None):
    """
    Get sequencing platform from SRA accession

    Args:
        srr_id: SRA/Run accession (SRR/ERR/DRR/CRR)
        bioproject_id: Optional BioProject ID (required for CNCB/CRR accessions)

    Returns:
        Platform name (e.g., 'ILLUMINA', 'OXFORD_NANOPORE') or None

    Note:
        - Works for NCBI accessions (SRR/ERR/DRR) via Entrez API
        - Works for CNCB accessions (CRR) via CNCB GSA API
        - CNCB queries require bioproject_id parameter (e.g., PRJCA040882)
        - Returns None if platform cannot be determined
    """
    # Check if this is a CRR accession (CNCB/China)
    if srr_id and srr_id.startswith('CRR'):
        return _get_platform_from_cncb(srr_id, bioproject_id)

    # Handle NCBI accessions (SRR/ERR/DRR)
    return _get_platform_from_ncbi(srr_id)


def _get_platform_from_cncb(crr_id, bioproject_id=None):
    """
    Get sequencing platform from CNCB/GSA for CRR accession

    Args:
        crr_id: CRR accession (e.g., CRR1878501)
        bioproject_id: Optional BioProject ID (e.g., PRJCA040882)
                      CNCB API requires BioProject to query run info

    Returns:
        Platform name or None

    Note:
        The CNCB getRunInfo API requires a BioProject accession, not a run accession.
        If bioproject_id is not provided, will attempt to extract from context.
    """
    import requests
    from io import StringIO

    BASE_URL = "https://ngdc.cncb.ac.cn/gsa"
    HEADERS = {"User-Agent": "Mozilla/5.0"}

    # Determine what to query - prefer BioProject ID
    search_term = bioproject_id if bioproject_id else crr_id

    try:
        # Query CNCB GSA API for run information
        # Note: This API expects a BioProject ID, not a run ID
        url = f"{BASE_URL}/search/getRunInfo"
        data = f'searchTerm=%26quot%3B{search_term}%26quot%3BtotalDatas=9999%3BdownLoadCount=9999'

        print(f"  Querying CNCB API with: {search_term}", file=sys.stderr)

        resp = requests.post(
            url,
            data=data,
            headers={**HEADERS, "Content-Type": "application/x-www-form-urlencoded"},
            timeout=30
        )
        resp.raise_for_status()

        # Parse CSV response
        csv_content = resp.text
        if csv_content.count('\n') < 2:
            print(f"Warning: No data returned from CNCB for {search_term}", file=sys.stderr)
            if bioproject_id is None:
                print(f"  Hint: CNCB API requires BioProject ID, not run ID", file=sys.stderr)
                print(f"  Try passing bioproject_id parameter", file=sys.stderr)
            return None

        # Read as DataFrame
        df = pd.read_csv(StringIO(csv_content))

        if df.empty:
            print(f"Warning: Empty response from CNCB for {search_term}", file=sys.stderr)
            return None

        print(f"  ✓ Retrieved {len(df)} runs from CNCB", file=sys.stderr)

        # If we have a specific run ID, filter for it
        if crr_id and 'Run' in df.columns:
            run_df = df[df['Run'] == crr_id]
            if not run_df.empty:
                df = run_df
                print(f"  ✓ Found metadata for run {crr_id}", file=sys.stderr)
            else:
                print(f"Warning: Run {crr_id} not found in BioProject results", file=sys.stderr)
                # Use first run's platform as fallback
                print(f"  Using platform from first run as fallback", file=sys.stderr)

        # Look for platform information in common column names
        platform_columns = ['Platform', 'Instrument', 'Model', 'Sequencing Platform']

        for col in platform_columns:
            if col in df.columns:
                platform_value = df[col].iloc[0] if not df[col].empty else None
                if platform_value and str(platform_value) != 'nan':
                    # Normalize platform names to match NCBI format
                    platform_str = str(platform_value).upper()

                    print(f"  Platform from CNCB: {platform_str}", file=sys.stderr)

                    # Map common platform names
                    if 'ILLUMINA' in platform_str or 'HISEQ' in platform_str or 'NOVASEQ' in platform_str or 'MISEQ' in platform_str:
                        return 'ILLUMINA'
                    elif 'NANOPORE' in platform_str or 'MINION' in platform_str or 'PROMETHION' in platform_str:
                        return 'OXFORD_NANOPORE'
                    elif 'PACBIO' in platform_str or 'SEQUEL' in platform_str:
                        return 'PACBIO_SMRT'
                    elif 'ION' in platform_str or 'TORRENT' in platform_str:
                        return 'ION_TORRENT'
                    elif '454' in platform_str or 'ROCHE' in platform_str:
                        return 'LS454'
                    else:
                        # Return the original value if no mapping found
                        return platform_str

        print(f"Warning: Platform column not found in CNCB metadata", file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        return None

    except Exception as e:
        print(f"Warning: Failed to retrieve platform from CNCB for {search_term}: {str(e)}", file=sys.stderr)
        return None


def _get_platform_from_ncbi(srr_id):
    """
    Get sequencing platform from NCBI for SRR/ERR/DRR accession

    Args:
        srr_id: NCBI SRA accession (SRR/ERR/DRR)

    Returns:
        Platform name or None
    """
    from Bio import Entrez
    import xml.etree.ElementTree as ET

    Entrez.email = "your_email@example.com"

    try:
        search_handle = Entrez.esearch(db="sra", term=srr_id)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if not search_results['IdList']:
            print(f"Warning: No results found for {srr_id}", file=sys.stderr)
            return None

        uid = search_results['IdList'][0]

        fetch_handle = Entrez.efetch(db="sra", id=uid, retmode="xml")
        xml_data = fetch_handle.read()
        fetch_handle.close()

        root = ET.fromstring(xml_data)
        platform = root.find('.//PLATFORM')

        if platform is not None and len(platform) > 0:
            return platform[0].tag

        print(f"Warning: Platform not found in metadata for {srr_id}", file=sys.stderr)
        return None

    except Exception as e:
        print(f"Warning: Failed to retrieve platform for {srr_id}: {str(e)}", file=sys.stderr)
        return None


if __name__ == "__main__":
    import argparse
    
    required_modules = [
        ("pandas", "pandas"),
        ("numpy", "numpy"),
        ("argparse", "argparse"),
        ("Bio", "biopython"),
        ("biom", "biom-format")
    ]
    
    for module, install_name in required_modules:
        check_and_install(module, install_name)
    
    parser = argparse.ArgumentParser(description="Process bioinformatics files.")
    parser.add_argument(
        "function",
        choices=[
            "GenerateSRAsFile",
            "GenerateDatasetsIDsFile",
            "trim_pos_454",
            "mk_manifest_SE",
            "mk_manifest_PE",
            "trim_pos_dada2",
            "trim_pos_deblur",
            "combine_stats_file",
            "add_prefix_to_file",
            "detect_primers_16s",
            "smart_trim_16s",
            "get_sequencing_platform"
        ],
        help="The function to execute."
    )
    parser.add_argument("--FilePath", help="Path to input file")
    parser.add_argument("--SequencingPlatform", help="Sequencing Platform column name")
    parser.add_argument("--Forward", help="Forward column name")
    parser.add_argument("--Reverse", help="Reverse column name")
    parser.add_argument("--Datasets_ID", help="Datasets_ID column name")
    parser.add_argument("--Bioproject", help="Bioproject column name")
    parser.add_argument("--SRA_Number", help="SRA_Number column name")
    parser.add_argument("--Biosample", help="Biosample column name")
    parser.add_argument("--region", help="Region column name")
    parser.add_argument("--filter_stats", help="Path to filter stats file")
    parser.add_argument("--deblur_stats", help="Path to deblur stats file")
    parser.add_argument("--dataset_name", help="Dataset ID")
    parser.add_argument("--In_fasta", help="Input FASTA file")
    parser.add_argument("--In_table", help="Input table file")
    parser.add_argument("--Prefix", help="Prefix for output files")
    parser.add_argument("--input_path", help="Input directory with FASTQ files")
    parser.add_argument("--output_path", help="Output directory")
    parser.add_argument("--OutputDir", help="Output directory for generated files")
    parser.add_argument("--tmp_path", default="/tmp", help="Temporary directory")
    parser.add_argument("--ref_path", default="./Meta2Data/docs/J01859.1.fna", help="E. coli 16S reference")
    parser.add_argument("--srr_id", help="SRA accession number")
    parser.add_argument("--bioproject_id", help="BioProject ID (required for CNCB/CRR accessions)")

    args = parser.parse_args()
    
    # Execute functions
    if args.function == "mk_manifest_SE":
        mk_manifest_SE(args.FilePath)
        
    elif args.function == "mk_manifest_PE":
        mk_manifest_PE(args.FilePath)
        
    elif args.function == "combine_stats_file":
        combine_stats_file(args.filter_stats, args.deblur_stats, args.dataset_name)
        
    elif args.function == "trim_pos_dada2":
        trim_pos_dada2(args.FilePath)
        
    elif args.function == "trim_pos_454":
        trim_pos_454(args.FilePath)
        
    elif args.function == "trim_pos_deblur":
        trim_pos_deblur(args.FilePath)
        
    elif args.function == "GenerateDatasetsIDsFile":
        # SequencingPlatform is optional - if not provided, only outputs BioProject IDs
        # OutputDir is optional - if not provided, uses input file directory
        GenerateDatasetsIDsFile(args.FilePath, args.Bioproject, args.SequencingPlatform, args.OutputDir)
        
    elif args.function == "GenerateSRAsFile":
        # Biosample and OutputDir are optional
        # If Biosample not provided, uses SRA_Number for sample naming
        GenerateSRAsFile(args.FilePath, args.Bioproject, args.SRA_Number, args.Biosample, args.OutputDir)
        
    elif args.function == "add_prefix_to_file":
        add_prefix_to_file(args.In_fasta, args.In_table, args.Prefix)

    elif args.function == "detect_primers_16s":
        primers_found, trim_length = detect_primers_16s(
            args.input_path, 
            args.tmp_path, 
            args.ref_path
        )
        if primers_found:
            print(f"TRIM:{trim_length}", file=sys.stdout)
            sys.exit(0)  # Success - primers detected
        else:
            sys.exit(1)  # No primers found

    elif args.function == "smart_trim_16s":
        smart_trim_16s(args.input_path, args.output_path, args.tmp_path, args.ref_path)

    elif args.function == "get_sequencing_platform":
        platform = get_sequencing_platform(args.srr_id, args.bioproject_id)
        if platform:
            print(platform)
        else:
            sys.exit(1)