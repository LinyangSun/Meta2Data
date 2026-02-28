#!/usr/bin/env python3

import sys
import subprocess
import os
import glob
import gzip
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
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
        df["rename"] = df[Bioproject] + '_' + df[Biosample]
    else:
        # Fallback: use SRA_Number (Run ID) when Biosample is not available
        df["rename"] = df[Bioproject] + '_' + df[SRA_Number]

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


def subset_meta_for_test(file_path, Bioproject, SRA_Number, output_dir=None, n=2):
    """
    Subset metadata for test mode: keep only the first N SRA entries per BioProject.

    Args:
        file_path: Path to metadata CSV
        Bioproject: Column name for BioProject ID
        SRA_Number: Column name for SRA/Run accession
        output_dir: (Optional) Directory to write output file. If not provided, uses input file directory
        n: Number of SRA entries to keep per BioProject (default: 2)

    Returns:
        Path to the subset CSV file
    """
    if output_dir:
        directory_path = output_dir
    else:
        directory_path = os.path.dirname(os.path.abspath(file_path))

    df = pd.read_csv(file_path)

    if Bioproject not in df.columns:
        print(f"Error: Column '{Bioproject}' not found in {file_path}", file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        sys.exit(1)
    if SRA_Number not in df.columns:
        print(f"Error: Column '{SRA_Number}' not found in {file_path}", file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        sys.exit(1)

    # Keep first N SRA entries per BioProject
    subset = df.groupby(Bioproject, sort=False).head(n)

    out_path = os.path.join(directory_path, "test_subset_meta.csv")
    subset.to_csv(out_path, index=False)

    total_bioprojects = subset[Bioproject].nunique()
    total_sras = len(subset)
    print(f"Test subset created: {total_bioprojects} BioProjects, {total_sras} SRA entries", file=sys.stderr)
    print(f"Subset file: {out_path}", file=sys.stderr)

    return out_path


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
    """Generate paired-end manifest file from a list of actual file paths."""
    df = pd.read_csv(file_path, sep='\t', header=None)
    dataset_name = os.path.basename(file_path).replace("-file.txt", "")

    # Group actual file paths by sample ID and read direction
    samples = {}
    for _, row in df.iterrows():
        filepath = row[0].strip()
        basename = os.path.basename(filepath)
        # Extract sample name and suffix: e.g. "PROJ_SAMPLE_1.fastq.gz" → ("PROJ_SAMPLE", "1.fastq.gz")
        sample, suffix = basename.rsplit('_', 1)
        if suffix.startswith('1.fastq'):
            samples.setdefault(sample, {})['forward'] = filepath
        elif suffix.startswith('2.fastq'):
            samples.setdefault(sample, {})['reverse'] = filepath

    rows = []
    for sample in sorted(samples):
        paths = samples[sample]
        if 'forward' in paths and 'reverse' in paths:
            rows.append({
                'sample-id': sample,
                'forward-absolute-filepath': paths['forward'],
                'reverse-absolute-filepath': paths['reverse'],
            })

    df1 = pd.DataFrame(rows)
    out_path = os.path.join(os.path.dirname(file_path), f"{dataset_name}_manifest.tsv")
    df1.to_csv(out_path, sep='\t', index=False)


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

    Returns:
        Platform name or None

    Strategy:
        1. Try querying by CRR run ID directly using getRunInfoByCra endpoint
        2. If fails, try querying by BioProject ID
        3. Parse response and extract platform information

    Note:
        This is a copy/adaptation of MetaDL's CNCB query logic for AmpliconPIP use.
        According to iSeq updates (2024), API endpoint changed from getRunInfo to getRunInfoByCra.
    """
    import requests
    from io import StringIO

    BASE_URL = "https://ngdc.cncb.ac.cn/gsa"
    HEADERS = {"User-Agent": "Mozilla/5.0"}

    # Strategy 1: Try querying by run ID directly with getRunInfoByCra endpoint
    print(f"  Attempting CNCB query with run ID: {crr_id}", file=sys.stderr)

    try:
        # Try updated getRunInfoByCra endpoint with CRR ID
        url = f"{BASE_URL}/search/getRunInfoByCra"
        data = f'searchTerm=%26quot%3B{crr_id}%26quot%3BtotalDatas=9999%3BdownLoadCount=9999'

        resp = requests.post(
            url,
            data=data,
            headers={**HEADERS, "Content-Type": "application/x-www-form-urlencoded"},
            timeout=30
        )
        resp.raise_for_status()

        csv_content = resp.text
        if csv_content.count('\n') >= 2:
            # Successfully got data with run ID
            print(f"  ✓ getRunInfoByCra with CRR ID successful", file=sys.stderr)
            platform = _parse_cncb_platform_response(csv_content, crr_id)
            if platform:
                return platform
    except Exception as e:
        print(f"  getRunInfoByCra with CRR failed: {str(e)}", file=sys.stderr)

    # Strategy 2: If BioProject ID provided, try querying with it
    if bioproject_id:
        print(f"  Attempting CNCB query with BioProject: {bioproject_id}", file=sys.stderr)

        try:
            # Try with BioProject using getRunInfoByCra
            url = f"{BASE_URL}/search/getRunInfoByCra"
            data = f'searchTerm=%26quot%3B{bioproject_id}%26quot%3BtotalDatas=9999%3BdownLoadCount=9999'

            resp = requests.post(
                url,
                data=data,
                headers={**HEADERS, "Content-Type": "application/x-www-form-urlencoded"},
                timeout=30
            )
            resp.raise_for_status()

            csv_content = resp.text
            if csv_content.count('\n') >= 2:
                print(f"  ✓ getRunInfoByCra with BioProject successful", file=sys.stderr)
                platform = _parse_cncb_platform_response(csv_content, crr_id)
                if platform:
                    return platform
        except Exception as e:
            print(f"  getRunInfoByCra with BioProject failed: {str(e)}", file=sys.stderr)

        # Strategy 3: Fallback to old getRunInfo endpoint with BioProject
        print(f"  Trying legacy getRunInfo endpoint with BioProject", file=sys.stderr)

        try:
            url = f"{BASE_URL}/search/getRunInfo"
            data = f'searchTerm=%26quot%3B{bioproject_id}%26quot%3BtotalDatas=9999%3BdownLoadCount=9999'

            resp = requests.post(
                url,
                data=data,
                headers={**HEADERS, "Content-Type": "application/x-www-form-urlencoded"},
                timeout=30
            )
            resp.raise_for_status()

            csv_content = resp.text
            if csv_content.count('\n') >= 2:
                print(f"  ✓ Legacy getRunInfo successful", file=sys.stderr)
                platform = _parse_cncb_platform_response(csv_content, crr_id)
                if platform:
                    return platform
        except Exception as e:
            print(f"  Legacy getRunInfo failed: {str(e)}", file=sys.stderr)

    print(f"Warning: All CNCB query strategies failed for {crr_id}", file=sys.stderr)
    return None


def _parse_cncb_platform_response(csv_content, target_run_id=None):
    """
    Parse CNCB API CSV response and extract platform information

    Args:
        csv_content: CSV string from CNCB API
        target_run_id: Optional specific run ID to filter for

    Returns:
        Platform name or None
    """
    from io import StringIO

    try:
        df = pd.read_csv(StringIO(csv_content))

        if df.empty:
            return None

        print(f"  ✓ Retrieved {len(df)} runs from CNCB", file=sys.stderr)

        # If we have a specific run ID, filter for it
        if target_run_id and 'Run' in df.columns:
            run_df = df[df['Run'] == target_run_id]
            if not run_df.empty:
                df = run_df
                print(f"  ✓ Found metadata for run {target_run_id}", file=sys.stderr)
            else:
                print(f"Warning: Run {target_run_id} not found, using first run as fallback", file=sys.stderr)

        # Look for platform information in common column names
        platform_columns = ['Platform', 'Instrument', 'Model', 'Sequencing Platform', 'instrument']

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

        print(f"Warning: Platform column not found in CNCB response", file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        return None

    except Exception as e:
        print(f"Warning: Failed to parse CNCB response: {str(e)}", file=sys.stderr)
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


def adaptive_tail_trim(input_dir, output_dir, max_sample_reads=10000):
    """
    Adaptive tail trimming for 454 reads: analyse → trim → compute max-ambiguous.

    454 pyrosequencing reads accumulate N bases toward the 3' end due to
    signal decay.  This function:
      1. Scans per-position N frequency from the 3' end to decide how many
         bases to trim (data-driven, not a fixed number).
      2. Trims that many bases from every read and writes the results to
         *output_dir*.
      3. Counts remaining N bases per read after trimming and returns the
         95th-percentile value for use as ``--p-max-ambiguous``.

    Args:
        input_dir:  Directory containing FASTQ files (gzipped or plain).
        output_dir: Directory for trimmed FASTQ output.
        max_sample_reads: Max reads to sample per file for statistics
                          (0 = unlimited).

    Returns:
        dict with keys ``trim_length`` (int) and ``max_ambiguous`` (int).

    Prints machine-readable lines to stdout:
        TRIM_LENGTH=<int>
        MAX_AMBIGUOUS=<int>
    """
    os.makedirs(output_dir, exist_ok=True)

    # ── Collect FASTQ paths ──────────────────────────────────────────────
    fq_files = sorted(
        glob.glob(os.path.join(input_dir, '*.fastq*'))
    )
    if not fq_files:
        raise FileNotFoundError(f"No FASTQ files found in {input_dir}")

    # ── Step 1: per-position N frequency from 3' end ─────────────────────
    SCAN_WINDOW = 50          # positions from 3' end to examine
    n_counts = np.zeros(SCAN_WINDOW, dtype=np.int64)   # N count at each pos
    total_counts = np.zeros(SCAN_WINDOW, dtype=np.int64)  # reads covering pos

    for fq in fq_files:
        open_fn = gzip.open if fq.endswith('.gz') else open
        sampled = 0
        with open_fn(fq, 'rt') as fh:
            while max_sample_reads == 0 or sampled < max_sample_reads:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().rstrip('\n')
                fh.readline()   # +
                fh.readline()   # qual
                sampled += 1
                seq_len = len(seq)
                limit = min(seq_len, SCAN_WINDOW)
                for pos in range(limit):
                    idx = seq_len - 1 - pos          # actual index in seq
                    total_counts[pos] += 1
                    if seq[idx] == 'N' or seq[idx] == 'n':
                        n_counts[pos] += 1

    # Compute N frequency at each position from 3' end
    with np.errstate(divide='ignore', invalid='ignore'):
        n_freq = np.where(total_counts > 0,
                          n_counts / total_counts, 0.0)

    # Background N frequency: average over positions 25-45 from 3' end
    bg_start, bg_end = 25, 45
    bg_positions = n_freq[bg_start:bg_end]
    if len(bg_positions) > 0 and np.any(total_counts[bg_start:bg_end] > 0):
        background = float(np.mean(bg_positions[total_counts[bg_start:bg_end] > 0]))
    else:
        background = 0.0

    threshold = max(background * 3, 0.01)

    # Scan from position 0 (3' end) inward
    trim_length = 0
    for pos in range(SCAN_WINDOW):
        if n_freq[pos] > threshold:
            trim_length = pos + 1
        else:
            break

    # Cap at 30 to avoid over-trimming
    trim_length = min(trim_length, 30)

    print(f"  Tail N-frequency analysis (first 15 positions from 3' end):",
          file=sys.stderr)
    for p in range(min(15, SCAN_WINDOW)):
        bar = '#' * int(n_freq[p] * 50)
        print(f"    pos -{p+1:>2d}: {n_freq[p]*100:5.1f}%  {bar}",
              file=sys.stderr)
    print(f"  Background N freq (pos -26 to -45): {background*100:.2f}%",
          file=sys.stderr)
    print(f"  Threshold: {threshold*100:.2f}%", file=sys.stderr)
    print(f"  → Trim length: {trim_length} bp", file=sys.stderr)

    # ── Step 2: trim tails and write output ──────────────────────────────
    n_per_read = []   # collect N counts for step 3

    for fq in fq_files:
        open_fn = gzip.open if fq.endswith('.gz') else open
        out_name = os.path.basename(fq)
        # Always write gzipped output
        if not out_name.endswith('.gz'):
            out_name += '.gz'
        out_path = os.path.join(output_dir, out_name)

        sampled_for_stats = 0
        with open_fn(fq, 'rt') as fin, gzip.open(out_path, 'wt') as fout:
            while True:
                header = fin.readline()
                if not header:
                    break
                seq = fin.readline().rstrip('\n')
                plus = fin.readline()
                qual = fin.readline().rstrip('\n')

                if trim_length > 0 and len(seq) > trim_length:
                    seq = seq[:-trim_length]
                    qual = qual[:-trim_length]

                fout.write(header)
                fout.write(seq + '\n')
                fout.write(plus)
                fout.write(qual + '\n')

                # Sample N counts for step 3
                if max_sample_reads == 0 or sampled_for_stats < max_sample_reads:
                    n_per_read.append(seq.upper().count('N'))
                    sampled_for_stats += 1

    # ── Step 3: P95 of per-read N count → max_ambiguous ──────────────────
    if n_per_read:
        p95 = int(np.percentile(n_per_read, 95))
    else:
        p95 = 0

    max_ambiguous = max(p95, 0)

    print(f"  Post-trim N distribution (sampled {len(n_per_read)} reads):",
          file=sys.stderr)
    n_counter = Counter(n_per_read)
    for k in sorted(n_counter.keys())[:10]:
        pct = n_counter[k] / len(n_per_read) * 100
        print(f"    {k} Ns: {n_counter[k]:>7d} ({pct:.1f}%)", file=sys.stderr)
    print(f"  P95 = {p95}  →  --p-max-ambiguous {max_ambiguous}",
          file=sys.stderr)

    # Machine-readable output on stdout
    print(f"TRIM_LENGTH={trim_length}")
    print(f"MAX_AMBIGUOUS={max_ambiguous}")

    return {"trim_length": trim_length, "max_ambiguous": max_ambiguous}


def check_quality_diversity(fastq_dir, n_samples=3, n_reads=1000):
    """
    Check quality score diversity in FASTQ files to determine if DADA2 is viable.

    DADA2 requires diverse quality scores to learn its error model. Binned or
    dummy quality scores (e.g., from NCBI fasterq-dump) have too few distinct
    values, causing DADA2 to fail or produce unreliable results.

    Checks up to n_samples FASTQ files, reading n_reads from each. Counts the
    number of unique Phred quality ASCII characters across all sampled reads.

    Threshold: >= 10 unique quality values → "normal" (DADA2 viable)
               <  10 unique quality values → "degraded" (use VSEARCH instead)

    Args:
        fastq_dir: Directory containing FASTQ files
        n_samples: Number of files to sample (default: 3)
        n_reads: Number of reads per file to sample (default: 1000)

    Prints machine-readable lines to stdout:
        QUALITY_STATUS=normal|degraded
        UNIQUE_QUALS=<int>
    """
    fq_files = sorted(glob.glob(os.path.join(fastq_dir, '*.fastq*')))
    if not fq_files:
        print("WARNING: No FASTQ files found, assuming degraded", file=sys.stderr)
        print("QUALITY_STATUS=degraded")
        print("UNIQUE_QUALS=0")
        return "degraded"

    sample_files = fq_files[:n_samples]
    all_qual_chars = set()

    for fq in sample_files:
        open_fn = gzip.open if fq.endswith('.gz') else open
        count = 0
        with open_fn(fq, 'rt') as fh:
            while count < n_reads:
                header = fh.readline()
                if not header:
                    break
                fh.readline()          # sequence
                fh.readline()          # +
                qual = fh.readline().strip()
                count += 1
                all_qual_chars.update(qual)

    n_unique = len(all_qual_chars)
    threshold = 10
    status = "normal" if n_unique >= threshold else "degraded"

    print(f"  Quality diversity: {n_unique} unique Q values from "
          f"{len(sample_files)} files x {n_reads} reads", file=sys.stderr)
    print(f"  Unique Q chars: {sorted(all_qual_chars)}", file=sys.stderr)
    print(f"  Status: {status} (threshold >= {threshold})", file=sys.stderr)

    print(f"QUALITY_STATUS={status}")
    print(f"UNIQUE_QUALS={n_unique}")
    return status


def _process_single_fastq(args):
    """Worker function for parallel FASTQ preprocessing (must be top-level for pickling)."""
    fq, output_dir, trim_front, truncate_length, max_n, is_paired = args
    open_fn = gzip.open if fq.endswith('.gz') else open
    out_name = os.path.basename(fq)
    if is_paired:
        out_name = out_name.replace('_1.fastq', '.fastq').replace('_R1', '')
    if not out_name.endswith('.gz'):
        out_name += '.gz'
    out_path = os.path.join(output_dir, out_name)

    file_in = 0
    file_out = 0
    file_short = 0
    file_n_filtered = 0

    with open_fn(fq, 'rt') as fin, gzip.open(out_path, 'wt') as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline().rstrip('\n')
            plus = fin.readline()
            qual = fin.readline().rstrip('\n')
            file_in += 1

            seq_after_front = seq[trim_front:]
            qual_after_front = qual[trim_front:]
            if len(seq_after_front) < truncate_length:
                file_short += 1
                continue
            seq_trimmed = seq_after_front[:truncate_length]
            qual_trimmed = qual_after_front[:truncate_length]

            n_count = seq_trimmed.upper().count('N')
            if n_count > max_n:
                file_n_filtered += 1
                continue

            fout.write(header)
            fout.write(seq_trimmed + '\n')
            fout.write(plus)
            fout.write(qual_trimmed + '\n')
            file_out += 1

    return os.path.basename(fq), file_in, file_out, file_short, file_n_filtered


def degraded_quality_preprocess(input_dir, output_dir, trim_front=15, truncate_length=0,
                                max_n=1, sequence_type="single", threads=4):
    """
    Preprocess FASTQ files for the degraded quality score pipeline.

    When quality scores are binned/unreliable (DADA2 cannot learn error model),
    this function applies simple deterministic preprocessing:
      1. Trim first trim_front bp from each read
      2. Truncate to a fixed length (all output reads are identical length)
      3. Discard reads with N count > max_n
      4. For PE data: use only forward (R1) reads to avoid merge Q value issues

    Files are processed in parallel using multiple workers.

    Args:
        input_dir: Directory with FASTQ files (after adapter removal + primer trimming)
        output_dir: Output directory for preprocessed files
        trim_front: Bases to trim from 5' end (default: 15)
        truncate_length: Fixed length to keep after trim_front (0=auto-detect
                         using 10th percentile of read lengths)
        max_n: Maximum N bases allowed per read; reads exceeding this are discarded (default: 1)
        sequence_type: "paired" or "single"
        threads: Number of parallel workers (default: 4)

    Prints machine-readable lines to stdout:
        TRUNCATE_LENGTH=<int>
        TOTAL_IN=<int>
        TOTAL_OUT=<int>
    """
    os.makedirs(output_dir, exist_ok=True)

    fq_files = sorted(glob.glob(os.path.join(input_dir, '*.fastq*')))
    if not fq_files:
        raise FileNotFoundError(f"No FASTQ files found in {input_dir}")

    # For PE: only process R1 files (forward reads only)
    if sequence_type == "paired":
        r1_files = [f for f in fq_files if '_1.fastq' in f or '_R1' in f]
        if not r1_files:
            print("  WARNING: No R1 files found for PE, using all files", file=sys.stderr)
            r1_files = fq_files
        process_files = r1_files
        print(f"  PE mode: using forward reads only ({len(r1_files)} R1 files)", file=sys.stderr)
    else:
        process_files = fq_files
        print(f"  SE mode: processing {len(process_files)} files", file=sys.stderr)

    # Auto-detect truncation length if not specified
    if truncate_length <= 0:
        print("  Auto-detecting truncation length...", file=sys.stderr)
        sampled_lengths = []
        sample_files = process_files[:min(5, len(process_files))]
        for fq in sample_files:
            open_fn = gzip.open if fq.endswith('.gz') else open
            count = 0
            with open_fn(fq, 'rt') as fin:
                while count < 2000:
                    header = fin.readline()
                    if not header:
                        break
                    seq = fin.readline().rstrip('\n')
                    _ = fin.readline()  # +
                    _ = fin.readline()  # qual
                    remaining = len(seq) - trim_front
                    if remaining > 0:
                        sampled_lengths.append(remaining)
                    count += 1

        if not sampled_lengths:
            raise ValueError("No reads found for truncation length detection")

        truncate_length = int(np.percentile(sampled_lengths, 10))
        if truncate_length < 50:
            print(f"  WARNING: Auto-detected truncation length is very short ({truncate_length}bp)",
                  file=sys.stderr)
        print(f"  Auto-detected truncation length: {truncate_length}bp "
              f"(10th percentile of {len(sampled_lengths)} sampled reads)", file=sys.stderr)
    else:
        print(f"  Using specified truncation length: {truncate_length}bp", file=sys.stderr)

    # Process files in parallel
    is_paired = (sequence_type == "paired")
    worker_args = [
        (fq, output_dir, trim_front, truncate_length, max_n, is_paired)
        for fq in process_files
    ]

    n_workers = min(len(process_files), max(1, threads))
    total_in = 0
    total_out = 0
    total_short = 0
    total_n_filtered = 0

    if n_workers > 1:
        print(f"  Processing {len(process_files)} files with {n_workers} parallel workers",
              file=sys.stderr)
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            for fname, f_in, f_out, f_short, f_nfilt in executor.map(
                    _process_single_fastq, worker_args):
                total_in += f_in
                total_out += f_out
                total_short += f_short
                total_n_filtered += f_nfilt
                print(f"  {fname}: {f_in} -> {f_out} reads", file=sys.stderr)
    else:
        for args in worker_args:
            fname, f_in, f_out, f_short, f_nfilt = _process_single_fastq(args)
            total_in += f_in
            total_out += f_out
            total_short += f_short
            total_n_filtered += f_nfilt
            print(f"  {fname}: {f_in} -> {f_out} reads", file=sys.stderr)

    print(f"\n  Summary: {total_in} reads in, {total_out} out", file=sys.stderr)
    print(f"  Filtered: {total_short} too short (<{truncate_length}bp after trim), "
          f"{total_n_filtered} N>{max_n}", file=sys.stderr)

    print(f"TRUNCATE_LENGTH={truncate_length}")
    print(f"TOTAL_IN={total_in}")
    print(f"TOTAL_OUT={total_out}")


def derep_fastq_for_vsearch(input_dir, output_fasta, threads=4):
    """
    Dereplicate preprocessed FASTQ files directly using vsearch CLI.

    Bypasses QIIME2 entirely: reads all preprocessed FASTQs, converts to FASTA,
    runs vsearch --derep_fulllength, and outputs size-annotated FASTA sorted by
    abundance (ready for vsearch --cluster_size).

    This replaces the previous 4-step QIIME2 round-trip:
      ImportFastqToQiime2 → QualityControlForQZA → LS454_Deduplication → ExportForVsearch

    Args:
        input_dir: Directory with preprocessed FASTQ files (from degraded_quality_preprocess)
        output_fasta: Output path for size-annotated dereplicated FASTA
        threads: Number of threads for vsearch (default: 4)

    Prints machine-readable lines to stdout:
        DEREP_INPUT_READS=<int>
        DEREP_UNIQUE_SEQS=<int>
    """
    import tempfile

    fq_files = sorted(glob.glob(os.path.join(input_dir, '*.fastq*')))
    if not fq_files:
        raise FileNotFoundError(f"No FASTQ files found in {input_dir}")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Step 1: Convert all FASTQ → single FASTA (vsearch needs FASTA input)
        combined_fasta = os.path.join(tmpdir, 'all_reads.fasta')
        total_reads = 0
        with open(combined_fasta, 'w') as fout:
            for fq in fq_files:
                open_fn = gzip.open if fq.endswith('.gz') else open
                with open_fn(fq, 'rt') as fin:
                    while True:
                        header = fin.readline()
                        if not header:
                            break
                        seq = fin.readline().rstrip('\n')
                        _ = fin.readline()  # +
                        _ = fin.readline()  # qual
                        total_reads += 1
                        fout.write(f">{header[1:]}{seq}\n")

        print(f"  Combined {total_reads} reads from {len(fq_files)} files", file=sys.stderr)

        # Step 2: vsearch --derep_fulllength
        derep_fasta = os.path.join(tmpdir, 'derep.fasta')
        cmd = [
            'vsearch', '--derep_fulllength', combined_fasta,
            '--output', derep_fasta,
            '--sizein', '--sizeout',
            '--minuniquesize', '1',
            '--threads', str(threads)
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"vsearch derep_fulllength failed:\n{result.stderr}")

        # Step 3: Sort by abundance descending (required by vsearch --cluster_size)
        cmd_sort = [
            'vsearch', '--sortbysize', derep_fasta,
            '--output', output_fasta,
            '--sizein', '--sizeout'
        ]
        result_sort = subprocess.run(cmd_sort, capture_output=True, text=True)
        if result_sort.returncode != 0:
            raise RuntimeError(f"vsearch sortbysize failed:\n{result_sort.stderr}")

        # Count unique sequences in output
        n_unique = 0
        with open(output_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    n_unique += 1

    print(f"  Dereplicated: {total_reads} reads → {n_unique} unique sequences", file=sys.stderr)
    print(f"DEREP_INPUT_READS={total_reads}")
    print(f"DEREP_UNIQUE_SEQS={n_unique}")


def export_derep_for_vsearch(repseq_qza, table_qza, output_fasta):
    """
    Export QIIME2 dereplicated artifacts to a size-annotated FASTA for vsearch.

    Reads the rep-seqs.qza (FASTA) and table.qza (BIOM feature table),
    annotates each sequence header with its total abundance (;size=N),
    and sorts by abundance descending (required by vsearch --cluster_size).

    Args:
        repseq_qza: Path to dereplicated rep-seqs .qza
        table_qza: Path to dereplicated table .qza
        output_fasta: Output path for size-annotated FASTA

    Prints machine-readable lines to stdout:
        EXPORT_FEATURES=<int>
        EXPORT_TOTAL_ABUNDANCE=<int>
    """
    import tempfile
    from biom import load_table

    with tempfile.TemporaryDirectory() as tmpdir:
        seq_dir = os.path.join(tmpdir, 'seqs')
        tbl_dir = os.path.join(tmpdir, 'table')

        # Export both .qza artifacts
        subprocess.run([
            'qiime', 'tools', 'export',
            '--input-path', repseq_qza, '--output-path', seq_dir
        ], check=True, capture_output=True)

        subprocess.run([
            'qiime', 'tools', 'export',
            '--input-path', table_qza, '--output-path', tbl_dir
        ], check=True, capture_output=True)

        # Read BIOM table → per-feature total abundance
        biom_path = os.path.join(tbl_dir, 'feature-table.biom')
        table = load_table(biom_path)
        feature_abundances = {}
        for fid in table.ids(axis='observation'):
            total = int(table.data(fid, axis='observation', dense=True).sum())
            feature_abundances[fid] = total

        # Read FASTA sequences
        fasta_path = os.path.join(seq_dir, 'dna-sequences.fasta')
        sequences = []
        for record in SeqIO.parse(fasta_path, 'fasta'):
            sequences.append((str(record.id), str(record.seq)))

        # Validate: FASTA IDs == BIOM feature IDs
        fasta_ids = set(sid for sid, _ in sequences)
        biom_ids = set(feature_abundances.keys())
        if fasta_ids != biom_ids:
            only_fasta = fasta_ids - biom_ids
            only_biom = biom_ids - fasta_ids
            raise ValueError(
                f"Feature ID mismatch between rep-seqs and table.\n"
                f"  In FASTA only ({len(only_fasta)}): {list(only_fasta)[:5]}\n"
                f"  In BIOM only ({len(only_biom)}): {list(only_biom)[:5]}"
            )

        # Sort by abundance descending (vsearch --cluster_size expects this)
        sequences.sort(key=lambda x: feature_abundances[x[0]], reverse=True)

        # Write size-annotated FASTA
        total_abundance = 0
        with open(output_fasta, 'w') as fout:
            for sid, seq in sequences:
                size = feature_abundances[sid]
                total_abundance += size
                fout.write(f">{sid};size={size}\n{seq}\n")

    print(f"  Exported {len(sequences)} features, total abundance {total_abundance}",
          file=sys.stderr)
    print(f"EXPORT_FEATURES={len(sequences)}")
    print(f"EXPORT_TOTAL_ABUNDANCE={total_abundance}")


def _relabel_single_sample(args):
    """Worker function for parallel sample relabeling (must be top-level for pickling)."""
    sample_id, filepath, output_path = args
    open_fn = gzip.open if filepath.endswith('.gz') else open
    read_num = 0
    with open_fn(filepath, 'rt') as fin, open(output_path, 'w') as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline().rstrip('\n')
            _ = fin.readline()  # +
            _ = fin.readline()  # qual
            read_num += 1
            fout.write(f">{sample_id};read_{read_num}\n{seq}\n")
    return sample_id, read_num


def relabel_reads_for_mapping(manifest_path, output_fasta, threads=4):
    """
    Merge all samples' preprocessed reads into a single FASTA with sample labels.

    Reads the QIIME2 manifest file to get authoritative sample-id → filepath
    mappings, converts each sample's FASTQ to FASTA with headers:
        >sample_id;read_1
        >sample_id;read_2
        ...

    Samples are processed in parallel, then concatenated in order.

    This combined FASTA is used by vsearch --usearch_global to build the OTU table.

    Args:
        manifest_path: Path to QIIME2 manifest TSV (sample-id, absolute-filepath)
        output_fasta: Output path for combined labeled FASTA
        threads: Number of parallel workers (default: 4)

    Prints machine-readable lines to stdout:
        RELABEL_SAMPLES=<int>
        RELABEL_TOTAL_READS=<int>
    """
    import tempfile

    manifest = pd.read_csv(manifest_path, sep='\t')
    n_samples = len(manifest)
    n_workers = min(n_samples, max(1, threads))

    with tempfile.TemporaryDirectory() as tmpdir:
        # Build worker args: each sample writes to a temp file
        worker_args = []
        tmp_paths = []
        for idx, row in manifest.iterrows():
            sample_id = str(row['sample-id'])
            filepath = str(row['absolute-filepath'])
            tmp_path = os.path.join(tmpdir, f"sample_{idx}.fasta")
            tmp_paths.append(tmp_path)
            worker_args.append((sample_id, filepath, tmp_path))

        # Process samples in parallel
        total_reads = 0
        if n_workers > 1:
            print(f"  Relabeling {n_samples} samples with {n_workers} parallel workers",
                  file=sys.stderr)
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                for sample_id, read_num in executor.map(_relabel_single_sample, worker_args):
                    total_reads += read_num
                    print(f"  {sample_id}: {read_num} reads", file=sys.stderr)
        else:
            for args in worker_args:
                sample_id, read_num = _relabel_single_sample(args)
                total_reads += read_num
                print(f"  {sample_id}: {read_num} reads", file=sys.stderr)

        # Concatenate temp files in manifest order
        with open(output_fasta, 'wb') as fout:
            for tmp_path in tmp_paths:
                with open(tmp_path, 'rb') as fin:
                    while True:
                        chunk = fin.read(1024 * 1024)  # 1MB chunks
                        if not chunk:
                            break
                        fout.write(chunk)

    print(f"\n  Total: {n_samples} samples, {total_reads} reads", file=sys.stderr)
    print(f"RELABEL_SAMPLES={n_samples}")
    print(f"RELABEL_TOTAL_READS={total_reads}")


def import_vsearch_to_qiime2(zotu_fasta, otu_table_tsv, manifest_path,
                              output_table_qza, output_repseq_qza):
    """
    Import vsearch results (ZOTU FASTA + OTU table) back into QIIME2 artifacts.

    Steps:
      1. Strip ;size= annotations from ZOTU FASTA
      2. Validate feature ID consistency (OTU table IDs subset of FASTA IDs)
      3. Validate sample name consistency (OTU table samples subset of manifest)
      4. Convert OTU table TSV → BIOM V2.1 (HDF5)
      5. Import BIOM → FeatureTable[Frequency] .qza
      6. Import FASTA → FeatureData[Sequence] .qza

    Args:
        zotu_fasta: Path to ZOTU FASTA (may contain ;size= annotations)
        otu_table_tsv: Path to vsearch --otutabout TSV
        manifest_path: Path to QIIME2 manifest TSV
        output_table_qza: Output path for feature table .qza
        output_repseq_qza: Output path for rep-seqs .qza

    Prints machine-readable lines to stdout:
        IMPORT_FEATURES=<int>
        IMPORT_SAMPLES=<int>
        IMPORT_TOTAL_READS=<int>
    """
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        # 1. Strip ;size= from ZOTU FASTA
        clean_fasta = os.path.join(tmpdir, 'zotus_clean.fasta')
        n_features = 0
        fasta_ids = set()
        with open(zotu_fasta, 'r') as fin, open(clean_fasta, 'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    # Remove ;size=N and anything after it
                    header = line.rstrip('\n')
                    seq_id = header[1:].split(';')[0]
                    fasta_ids.add(seq_id)
                    fout.write(f">{seq_id}\n")
                    n_features += 1
                else:
                    fout.write(line)

        # Verify no ;size= remains
        with open(clean_fasta, 'r') as f:
            for line in f:
                if line.startswith('>') and ';size=' in line:
                    raise ValueError(f"Residual ;size= in cleaned FASTA: {line.rstrip()}")

        print(f"  Cleaned FASTA: {n_features} ZOTUs", file=sys.stderr)

        # 2. Read OTU table TSV to validate
        table_feature_ids = set()
        table_sample_ids = set()
        total_reads = 0
        with open(otu_table_tsv, 'r') as f:
            header_line = None
            for line in f:
                if line.startswith('#'):
                    continue
                if header_line is None:
                    header_line = line.rstrip('\n').split('\t')
                    table_sample_ids = set(header_line[1:])
                    continue
                parts = line.rstrip('\n').split('\t')
                # Feature ID may contain ;size= from vsearch output
                feature_id = parts[0].split(';')[0]
                table_feature_ids.add(feature_id)
                total_reads += sum(int(float(x)) for x in parts[1:] if x)

        # 3. Validate feature IDs: table features must be in FASTA
        missing_features = table_feature_ids - fasta_ids
        if missing_features:
            raise ValueError(
                f"OTU table has {len(missing_features)} features not in ZOTU FASTA: "
                f"{list(missing_features)[:5]}"
            )

        # 4. Validate sample names against manifest
        manifest = pd.read_csv(manifest_path, sep='\t')
        manifest_samples = set(manifest['sample-id'].astype(str))
        unknown_samples = table_sample_ids - manifest_samples
        if unknown_samples:
            print(f"  WARNING: {len(unknown_samples)} samples in OTU table not in manifest: "
                  f"{list(unknown_samples)[:5]}", file=sys.stderr)

        print(f"  OTU table: {len(table_feature_ids)} features, "
              f"{len(table_sample_ids)} samples, {total_reads} total reads", file=sys.stderr)

        # 5. Convert TSV → BIOM V2.1 (HDF5)
        biom_path = os.path.join(tmpdir, 'otu_table.biom')
        subprocess.run([
            'biom', 'convert',
            '-i', otu_table_tsv,
            '-o', biom_path,
            '--table-type', 'OTU table',
            '--to-hdf5'
        ], check=True, capture_output=True)

        # Validate BIOM
        result = subprocess.run(
            ['biom', 'validate-table', '-i', biom_path],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            raise ValueError(f"BIOM validation failed: {result.stderr}")
        print(f"  BIOM validation passed", file=sys.stderr)

        # 6. Import into QIIME2
        subprocess.run([
            'qiime', 'tools', 'import',
            '--input-path', biom_path,
            '--type', 'FeatureTable[Frequency]',
            '--input-format', 'BIOMV210Format',
            '--output-path', output_table_qza
        ], check=True, capture_output=True)

        subprocess.run([
            'qiime', 'tools', 'import',
            '--input-path', clean_fasta,
            '--type', 'FeatureData[Sequence]',
            '--output-path', output_repseq_qza
        ], check=True, capture_output=True)

    print(f"  Imported to QIIME2: {output_table_qza}", file=sys.stderr)
    print(f"  Imported to QIIME2: {output_repseq_qza}", file=sys.stderr)
    print(f"IMPORT_FEATURES={len(table_feature_ids)}")
    print(f"IMPORT_SAMPLES={len(table_sample_ids)}")
    print(f"IMPORT_TOTAL_READS={total_reads}")


def append_summary(dataset_id, sra_file, raw_counts_file, final_table, output_csv):
    """
    Append per-sample summary (raw reads + final reads) to a unified CSV.

    Args:
        dataset_id: BioProject ID
        sra_file: Path to _sra.txt (Run<tab>SampleName)
        raw_counts_file: Path to _raw_read_counts.tsv (Run<tab>SampleName<tab>RawReads)
        final_table: Path to final-table.qza
        output_csv: Path to the unified summary CSV
    """
    import tempfile
    import subprocess
    from biom import load_table

    # 1. Read SRA mapping: Run → SampleName
    sra_df = pd.read_csv(sra_file, sep='\t', header=None, names=['Run', 'SampleName'])

    # 2. Read raw read counts
    raw_df = pd.read_csv(raw_counts_file, sep='\t', header=None, names=['Run', 'SampleName', 'RawReads'])

    # 3. Export final-table.qza and get per-sample total reads
    final_reads = {}
    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.run([
            'qiime', 'tools', 'export',
            '--input-path', final_table,
            '--output-path', tmpdir
        ], check=True, capture_output=True)

        biom_path = os.path.join(tmpdir, 'feature-table.biom')
        table = load_table(biom_path)

        for sample_id in table.ids(axis='sample'):
            total = int(table.data(sample_id, axis='sample', dense=True).sum())
            final_reads[sample_id] = total
            # Handle .fastq suffix artifact from SE manifest
            clean_id = sample_id.replace('.fastq', '')
            if clean_id != sample_id:
                final_reads[clean_id] = total

    # 4. Build summary rows
    rows = []
    for _, row in sra_df.iterrows():
        run_id = row['Run']
        sample_name = row['SampleName']

        raw_match = raw_df[raw_df['Run'] == run_id]
        raw_read_count = int(raw_match['RawReads'].iloc[0]) if not raw_match.empty else 0

        final_read_count = final_reads.get(sample_name, 0)

        rows.append({
            'BioProject': dataset_id,
            'Run': run_id,
            'SampleName': sample_name,
            'RawReads': raw_read_count,
            'FinalReads': final_read_count
        })

    result_df = pd.DataFrame(rows)

    # 5. Append to unified CSV (write header only if file doesn't exist or is empty)
    # Use file lock to prevent interleaving when datasets run in parallel
    import fcntl
    lock_path = output_csv + '.lock'
    with open(lock_path, 'w') as lockf:
        fcntl.flock(lockf, fcntl.LOCK_EX)
        write_header = not os.path.exists(output_csv) or os.path.getsize(output_csv) == 0
        result_df.to_csv(output_csv, mode='a', header=write_header, index=False)
        fcntl.flock(lockf, fcntl.LOCK_UN)

    print(f"✓ Summary appended for {dataset_id}: {len(rows)} samples")


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
            "subset_meta_for_test",
            "mk_manifest_SE",
            "mk_manifest_PE",
            "trim_pos_deblur",
            "get_sequencing_platform",
            "adaptive_tail_trim",
            "check_quality_diversity",
            "degraded_quality_preprocess",
            "derep_fastq_for_vsearch",
            "export_derep_for_vsearch",
            "relabel_reads_for_mapping",
            "import_vsearch_to_qiime2",
            "append_summary"
        ],
        help="The function to execute."
    )
    parser.add_argument("--FilePath", help="Path to input file")
    parser.add_argument("--SequencingPlatform", help="Sequencing Platform column name")
    parser.add_argument("--Bioproject", help="Bioproject column name")
    parser.add_argument("--SRA_Number", help="SRA_Number column name")
    parser.add_argument("--Biosample", help="Biosample column name")
    parser.add_argument("--OutputDir", help="Output directory for generated files")
    parser.add_argument("--srr_id", help="SRA accession number")
    parser.add_argument("--bioproject_id", help="BioProject ID (required for CNCB/CRR accessions)")
    parser.add_argument("--input_dir", help="Input directory (for adaptive_tail_trim)")
    parser.add_argument("--output_dir", help="Output directory (for adaptive_tail_trim)")
    parser.add_argument("--max_sample_reads", type=int, default=10000,
                        help="Max reads to sample per file (default: 10000, 0=all)")
    parser.add_argument("--dataset_id", help="Dataset/BioProject ID for summary")
    parser.add_argument("--sra_file", help="Path to _sra.txt file")
    parser.add_argument("--raw_counts", help="Path to raw read counts TSV")
    parser.add_argument("--final_table", help="Path to final-table.qza")
    parser.add_argument("--output_csv", help="Path to output summary CSV")
    parser.add_argument("--n_samples", type=int, default=3,
                        help="Number of files to sample for quality check (default: 3)")
    parser.add_argument("--n_reads", type=int, default=1000,
                        help="Number of reads per file to sample (default: 1000)")
    parser.add_argument("--trim_front", type=int, default=15,
                        help="Bases to trim from 5' end (default: 15)")
    parser.add_argument("--trim_tail", type=int, default=30,
                        help="Bases to trim from 3' end (default: 30)")
    parser.add_argument("--truncate_length", type=int, default=0,
                        help="Fixed truncation length after trim_front (0=auto-detect)")
    parser.add_argument("--max_n", type=int, default=1,
                        help="Max N bases allowed per read (default: 1)")
    parser.add_argument("--sequence_type", default="single",
                        help="Sequence type: 'paired' or 'single' (default: single)")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads for parallel operations (default: 4)")
    # Arguments for vsearch pipeline functions
    parser.add_argument("--repseq_qza", help="Path to rep-seqs .qza (for export_derep_for_vsearch)")
    parser.add_argument("--table_qza", help="Path to table .qza (for export_derep_for_vsearch)")
    parser.add_argument("--output_fasta", help="Output FASTA path")
    parser.add_argument("--manifest_path", help="Path to QIIME2 manifest TSV")
    parser.add_argument("--zotu_fasta", help="Path to ZOTU FASTA (for import_vsearch_to_qiime2)")
    parser.add_argument("--otu_table_tsv", help="Path to vsearch OTU table TSV")
    parser.add_argument("--output_table_qza", help="Output table .qza path")
    parser.add_argument("--output_repseq_qza", help="Output rep-seqs .qza path")

    args = parser.parse_args()
    
    # Execute functions
    if args.function == "mk_manifest_SE":
        mk_manifest_SE(args.FilePath)

    elif args.function == "mk_manifest_PE":
        mk_manifest_PE(args.FilePath)

    elif args.function == "trim_pos_deblur":
        trim_pos_deblur(args.FilePath)

    elif args.function == "GenerateDatasetsIDsFile":
        GenerateDatasetsIDsFile(args.FilePath, args.Bioproject, args.SequencingPlatform, args.OutputDir)

    elif args.function == "GenerateSRAsFile":
        GenerateSRAsFile(args.FilePath, args.Bioproject, args.SRA_Number, args.Biosample, args.OutputDir)

    elif args.function == "subset_meta_for_test":
        result = subset_meta_for_test(args.FilePath, args.Bioproject, args.SRA_Number, args.OutputDir)
        print(result)

    elif args.function == "get_sequencing_platform":
        platform = get_sequencing_platform(args.srr_id, args.bioproject_id)
        if platform:
            print(platform)
        else:
            sys.exit(1)

    elif args.function == "adaptive_tail_trim":
        adaptive_tail_trim(args.input_dir, args.output_dir, args.max_sample_reads)

    elif args.function == "check_quality_diversity":
        check_quality_diversity(args.input_dir, args.n_samples, args.n_reads)

    elif args.function == "degraded_quality_preprocess":
        degraded_quality_preprocess(args.input_dir, args.output_dir,
                                    args.trim_front, args.truncate_length,
                                    args.max_n, args.sequence_type,
                                    args.threads)

    elif args.function == "derep_fastq_for_vsearch":
        derep_fastq_for_vsearch(args.input_dir, args.output_fasta, args.threads)

    elif args.function == "export_derep_for_vsearch":
        export_derep_for_vsearch(args.repseq_qza, args.table_qza, args.output_fasta)

    elif args.function == "relabel_reads_for_mapping":
        relabel_reads_for_mapping(args.manifest_path, args.output_fasta, args.threads)

    elif args.function == "import_vsearch_to_qiime2":
        import_vsearch_to_qiime2(args.zotu_fasta, args.otu_table_tsv,
                                  args.manifest_path,
                                  args.output_table_qza, args.output_repseq_qza)

    elif args.function == "append_summary":
        append_summary(args.dataset_id, args.sra_file, args.raw_counts, args.final_table, args.output_csv)