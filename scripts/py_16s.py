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


def _clean_id_series(series):
    """Strip spaces, newlines and tabs from an ID column (as strings).

    Centralizes the identical .astype(str)+strip chain previously inlined in
    GenerateDatasetsIDsFile and GenerateSRAsFile; the transformation is purely
    on the Series value, so callers assigning the result back are unchanged.
    """
    return (
        series
        .astype(str)
        .str.replace(' ', '', regex=True)
        .str.replace('\n', '', regex=True)
        .str.replace('\t', '', regex=True)
    )


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
    if output_dir:
        directory_path = output_dir
    else:
        directory_path = os.path.dirname(os.path.abspath(file_path))

    df = pd.read_csv(file_path)

    if Bioproject in df.columns:
        df[Bioproject] = _clean_id_series(df[Bioproject])

    # If platform column is provided and exists, include it (backward compatibility)
    if Data_SequencingPlatform and Data_SequencingPlatform in df.columns:
        df[Data_SequencingPlatform] = _clean_id_series(df[Data_SequencingPlatform])
        df_pair = (
            df[[Bioproject, Data_SequencingPlatform]]
            .dropna(subset=[Bioproject])
            .drop_duplicates()
        )
    else:
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
    if output_dir:
        directory_path = output_dir
    else:
        directory_path = os.path.dirname(os.path.abspath(file_path))
    df = pd.read_csv(file_path)

    columns_to_clean = [Bioproject, SRA_Number]
    if Biosample and Biosample in df.columns:
        columns_to_clean.append(Biosample)

    for col in columns_to_clean:
        if col in df.columns:
            df[col] = _clean_id_series(df[col])

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
        sample = basenames.split('.fastq')[0]
        # SE reads downloaded as "<bioproject>_<run>_1.fastq" carry a read-pair
        # suffix (e.g. ONT). Strip a trailing "_1" so the sample-id matches the
        # "<bioproject>_<run>" name used by Common_CountRawReads / append_summary;
        # otherwise the OTU-table sample (with _1) never matches the summary's
        # SampleName (without _1) and FinalReads is mis-reported as 0.
        if sample.endswith('_1'):
            sample = sample[:-2]
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
    """
    Calculate trim positions (start, end) for DADA2/Deblur from a
    QIIME2 seven-number-summaries TSV.

    Strategy:
      - Use count row to find the read-length boundary (outer1).
      - Use the 25th-percentile quality row as the primary decision maker
        (represents the quality that 75% of reads exceed).
      - Apply a sliding window (W=5) to smooth single-position noise.
      - Scan inward from both ends to find positions where the smoothed
        25th-percentile quality consistently meets Q_TRIM (25).
      - Require 5 consecutive passing positions for stability.
      - Enforce a minimum retained length of 50 bp.

    Returns:
      Prints "start,end" to stdout for shell capture.
      Returns (start, end) tuple, or (None, None) if no valid window.
    """
    Q_TRIM = 25          # quality threshold for 25th percentile
    W = 5                # sliding window size
    CONSEC = 5           # consecutive positions required
    MIN_RETAIN = 50      # minimum retained sequence length

    with open(file_path, "r") as f:
        tsv = [line.rstrip("\n") for line in f]

    pos = [int(x) for x in tsv[0].split("\t")[1:]]
    L = len(pos)
    counts = [float(x) for x in tsv[1].split("\t")[1:]]

    # ── Step 1: count-based read-length boundary ──
    count_thresh = min(9000, max(counts) * 0.7)
    under = [i for i, v in enumerate(counts) if v < count_thresh]
    outer1 = max(under[0] if under else (L - 1), 40)

    # ── Step 2: parse the 25th-percentile row (index 2 in data rows) ──
    data_rows = []
    for line in tsv[3:]:
        vals = [float(x) for x in line.split("\t")[1:][:L]]
        data_rows.append(vals)

    if len(data_rows) < 3:
        print("None,None")
        return (None, None)

    pct_25 = data_rows[1]  # row order: 9%, 25%, 50%, 75%, ... (2% skipped in tsv[3:])

    # ── Step 3: sliding-window smoothing ──
    sm = [None] * L
    for i in range(L - W + 1):
        sm[i] = sum(pct_25[i:i + W]) / W

    # ── Step 4: find start (5' trim) ──
    # Scan left-to-right within [0, outer1) for CONSEC consecutive
    # positions where smoothed quality >= Q_TRIM.
    final_start = None
    for i in range(outer1 - CONSEC + 1):
        if sm[i] is None:
            continue
        if all(sm[i + j] is not None and sm[i + j] >= Q_TRIM
               for j in range(CONSEC)):
            final_start = i
            break

    # ── Step 5: find end (3' truncation) ──
    # Scan right-to-left from outer1 for CONSEC consecutive positions
    # where smoothed quality >= Q_TRIM.
    final_end = None
    scan_limit = L - W  # last index with a valid smoothed value
    for i in range(min(outer1, scan_limit), CONSEC - 1, -1):
        if all(sm[i - j] is not None and sm[i - j] >= Q_TRIM
               for j in range(CONSEC)):
            final_end = i
            break

    # ── Step 6: validate ──
    if final_start is None or final_end is None:
        print("None,None")
        return (None, None)

    if final_end - final_start < MIN_RETAIN:
        print("None,None")
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
    if srr_id and srr_id.startswith('CRR'):
        return _get_platform_from_cncb(srr_id, bioproject_id)

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
        CNCB's API endpoint was changed from getRunInfo to getRunInfoByCra (2024).
    """
    import requests

    BASE_URL = "https://ngdc.cncb.ac.cn/gsa"
    HEADERS = {"User-Agent": "Mozilla/5.0"}

    def _try_cncb(endpoint, search_term, success_label, fail_label):
        """Run one CNCB query strategy; return platform name or None.

        Differs across strategies only by endpoint, search term, and log
        labels; the request/parse logic is identical, so each former inline
        block maps to one call of this closure with byte-identical effects.
        """
        try:
            url = f"{BASE_URL}/search/{endpoint}"
            data = f'searchTerm=%26quot%3B{search_term}%26quot%3BtotalDatas=9999%3BdownLoadCount=9999'

            resp = requests.post(
                url,
                data=data,
                headers={**HEADERS, "Content-Type": "application/x-www-form-urlencoded"},
                timeout=30
            )
            resp.raise_for_status()

            csv_content = resp.text
            if csv_content.count('\n') >= 2:
                print(f"  ✓ {success_label}", file=sys.stderr)
                platform = _parse_cncb_platform_response(csv_content, crr_id)
                if platform:
                    return platform
        except Exception as e:
            print(f"  {fail_label}: {str(e)}", file=sys.stderr)
        return None

    # Strategy 1: Try querying by run ID directly with getRunInfoByCra endpoint
    print(f"  Attempting CNCB query with run ID: {crr_id}", file=sys.stderr)

    platform = _try_cncb("getRunInfoByCra", crr_id,
                         "getRunInfoByCra with CRR ID successful",
                         "getRunInfoByCra with CRR failed")
    if platform:
        return platform

    # Strategy 2: If BioProject ID provided, try querying with it
    if bioproject_id:
        print(f"  Attempting CNCB query with BioProject: {bioproject_id}", file=sys.stderr)

        platform = _try_cncb("getRunInfoByCra", bioproject_id,
                             "getRunInfoByCra with BioProject successful",
                             "getRunInfoByCra with BioProject failed")
        if platform:
            return platform

        # Strategy 3: Fallback to old getRunInfo endpoint with BioProject
        print(f"  Trying legacy getRunInfo endpoint with BioProject", file=sys.stderr)

        platform = _try_cncb("getRunInfo", bioproject_id,
                             "Legacy getRunInfo successful",
                             "Legacy getRunInfo failed")
        if platform:
            return platform

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

        if target_run_id and 'Run' in df.columns:
            run_df = df[df['Run'] == target_run_id]
            if not run_df.empty:
                df = run_df
                print(f"  ✓ Found metadata for run {target_run_id}", file=sys.stderr)
            else:
                print(f"Warning: Run {target_run_id} not found, using first run as fallback", file=sys.stderr)

        platform_columns = ['Platform', 'Instrument', 'Model', 'Sequencing Platform', 'instrument']

        for col in platform_columns:
            if col in df.columns:
                platform_value = df[col].iloc[0] if not df[col].empty else None
                if platform_value and str(platform_value) != 'nan':
                    platform_str = str(platform_value).upper()

                    print(f"  Platform from CNCB: {platform_str}", file=sys.stderr)

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
                        return platform_str

        print(f"Warning: Platform column not found in CNCB response", file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        return None

    except Exception as e:
        print(f"Warning: Failed to parse CNCB response: {str(e)}", file=sys.stderr)
        return None


def _configure_entrez():
    """Configure Biopython Entrez for polite, rate-limit-tolerant NCBI access.

    - NCBI_API_KEY (env): raises the rate limit from 3 to 10 requests/s.
    - NCBI_EMAIL (env):   identifies the caller to NCBI (recommended).
    - max_tries/sleep_between_tries: Biopython's built-in retry on transient
      HTTP errors (complements the explicit 429 backoff in callers).
    """
    import os
    from Bio import Entrez
    Entrez.email = os.environ.get("NCBI_EMAIL", "your_email@example.com")
    _key = os.environ.get("NCBI_API_KEY")
    if _key:
        Entrez.api_key = _key
    Entrez.max_tries = 4
    Entrez.sleep_between_tries = 15


def _get_platform_from_ncbi(srr_id):
    """
    Get sequencing platform from NCBI for SRR/ERR/DRR accession

    Args:
        srr_id: NCBI SRA accession (SRR/ERR/DRR)

    Returns:
        Platform name or None
    """
    import os, time
    from Bio import Entrez
    import xml.etree.ElementTree as ET
    from urllib.error import HTTPError

    _configure_entrez()

    # NCBI rate-limits to 3 req/s without an API key (10 with one). Under
    # parallel datasets this fallback can hit HTTP 429; retry with exponential
    # backoff + jitter instead of failing the whole dataset.
    max_retries = 5
    for attempt in range(max_retries):
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

        except HTTPError as e:
            if e.code == 429 and attempt < max_retries - 1:
                wait = 3 * (2 ** attempt) + (hash(srr_id) % 5)   # 3-8, 6-11, 12-17, ... s (jittered)
                print(f"Warning: 429 for {srr_id}, retrying in {wait}s "
                      f"(attempt {attempt + 1}/{max_retries})", file=sys.stderr)
                time.sleep(wait)
                continue
            print(f"Warning: Failed to retrieve platform for {srr_id}: {str(e)}", file=sys.stderr)
            return None
        except Exception as e:
            print(f"Warning: Failed to retrieve platform for {srr_id}: {str(e)}", file=sys.stderr)
            return None
    return None


def batch_get_sequencing_platforms(pairs_file):
    """
    Batch-detect sequencing platforms for multiple datasets.

    Reads a TSV file with dataset_id<TAB>srr_id[<TAB>bioproject_id] per line.
    NCBI accessions (SRR/ERR/DRR) are queried in a single Entrez epost+efetch.
    CNCB accessions (CRR) are queried serially with a short delay.

    Prints results to stdout as: dataset_id<TAB>platform (one per line).
    """
    from Bio import Entrez
    import xml.etree.ElementTree as ET
    import time

    _configure_entrez()

    ncbi_pairs = []   # [(dataset_id, srr_id), ...]
    cncb_pairs = []   # [(dataset_id, srr_id, bioproject_id), ...]

    with open(pairs_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            ds_id, srr = parts[0], parts[1]
            bio_id = parts[2] if len(parts) > 2 else None
            if srr.startswith('CRR'):
                cncb_pairs.append((ds_id, srr, bio_id))
            else:
                ncbi_pairs.append((ds_id, srr))

    # --- NCBI batch query via epost + efetch ---
    if ncbi_pairs:
        srr_to_ds = {srr: ds_id for ds_id, srr in ncbi_pairs}
        srr_list = list(srr_to_ds.keys())

        try:
            search_term = " OR ".join(srr_list)
            search_handle = Entrez.esearch(db="sra", term=search_term, retmax=len(srr_list))
            search_results = Entrez.read(search_handle)
            search_handle.close()

            uid_list = search_results.get('IdList', [])
            if uid_list:
                fetch_handle = Entrez.efetch(db="sra", id=",".join(uid_list), retmode="xml")
                xml_data = fetch_handle.read()
                fetch_handle.close()

                root = ET.fromstring(xml_data)
                for pkg in root.findall('.//EXPERIMENT_PACKAGE'):
                    run_el = pkg.find('.//RUN')
                    platform_el = pkg.find('.//PLATFORM')
                    if run_el is not None and platform_el is not None and len(platform_el) > 0:
                        acc = run_el.get('accession', '')
                        plat = platform_el[0].tag
                        if acc in srr_to_ds:
                            print(f"{srr_to_ds[acc]}\t{plat}")
                            del srr_to_ds[acc]

            for srr, ds_id in srr_to_ds.items():
                print(f"Warning: No platform found for {srr} ({ds_id})", file=sys.stderr)

        except Exception as e:
            print(f"Warning: Batch NCBI query failed: {e}", file=sys.stderr)

    # --- CNCB serial query with delay ---
    for ds_id, srr, bio_id in cncb_pairs:
        plat = _get_platform_from_cncb(srr, bio_id)
        if plat:
            print(f"{ds_id}\t{plat}")
        else:
            print(f"Warning: No platform found for {srr} ({ds_id})", file=sys.stderr)
        if len(cncb_pairs) > 1:
            time.sleep(2)


def detect_length_window(input_dir, tolerance=0.15, floor=200, max_sample_reads=10000):
    """
    Auto-detect the amplicon length window for ONT reads via peak finding.

    Amplicon data (16S/ITS/18S, any platform) has a sharp read-length peak at
    the target amplicon size. This samples reads, builds a length histogram,
    locates the dominant peak (mode), and returns a window around it:
        [peak*(1-tolerance), peak*(1+tolerance)]   with a hard lower floor.

    This replaces ONT-AmpSeq's manually-set 1200-1600bp window so the pipeline
    stays fully automatic and adapts to whatever amplicon a dataset contains.

    Args:
        input_dir:        Directory with FASTQ files (gzipped or plain).
        tolerance:        Fractional half-width around the peak (0.15 = +/-15%).
        floor:            Hard minimum for the lower bound (bp); guards against
                          short junk reads dragging the window to zero.
        max_sample_reads: Max reads to sample per file (0 = unlimited).

    Prints machine-readable lines to stdout:
        LENGTH_LO=<int>
        LENGTH_HI=<int>
        PEAK=<int>
        PEAK_STATUS=<strong|weak>
    """
    fq_files = sorted(glob.glob(os.path.join(input_dir, '*.fastq*')))
    if not fq_files:
        raise FileNotFoundError(f"No FASTQ files found in {input_dir}")

    BIN = 50  # histogram bin width (bp)
    lengths = []
    for fq in fq_files:
        open_fn = gzip.open if fq.endswith('.gz') else open
        sampled = 0
        with open_fn(fq, 'rt') as fh:
            while max_sample_reads == 0 or sampled < max_sample_reads:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().rstrip('\n')
                fh.readline()  # +
                fh.readline()  # qual
                if not seq:
                    continue
                lengths.append(len(seq))
                sampled += 1

    if not lengths:
        raise ValueError("No reads found for length window detection")

    arr = np.asarray(lengths, dtype=np.int64)
    max_len = int(arr.max())
    hist = np.bincount(arr // BIN, minlength=max(1, max_len // BIN + 1))
    peak_bin = int(hist.argmax())
    peak = int(peak_bin * BIN + BIN // 2)  # bin midpoint

    # Peak-dominance check: fraction of reads within the peak bin +/- 1 bin.
    near = int(hist[max(0, peak_bin - 1):peak_bin + 2].sum())
    frac = near / len(arr)
    if frac >= 0.25:
        status = "strong"
        center = peak
    else:
        # Weak / multi-modal length distribution: fall back to the median for a
        # robust center rather than trusting a marginal histogram peak.
        status = "weak"
        center = int(np.median(arr))

    lo = int(round(center * (1.0 - tolerance)))
    hi = int(round(center * (1.0 + tolerance)))
    lo = max(lo, int(floor))
    if hi <= lo:
        hi = lo + BIN

    print(f"  Sampled {len(arr)} reads from {len(fq_files)} file(s)", file=sys.stderr)
    print(f"  Lengths: min={int(arr.min())} median={int(np.median(arr))} "
          f"max={max_len}", file=sys.stderr)
    print(f"  Dominant peak ~{peak}bp ({frac * 100:.0f}% of reads near peak, "
          f"{status})", file=sys.stderr)
    print(f"  Window: {lo}-{hi}bp (center {center}, tolerance +/-"
          f"{tolerance * 100:.0f}%, floor {floor})", file=sys.stderr)
    print(f"LENGTH_LO={lo}")
    print(f"LENGTH_HI={hi}")
    print(f"PEAK={peak}")
    print(f"PEAK_STATUS={status}")
    return {"length_lo": lo, "length_hi": hi, "peak": peak, "status": status}


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

    print(f"TRIM_LENGTH={trim_length}")
    print(f"MAX_AMBIGUOUS={max_ambiguous}")

    return {"trim_length": trim_length, "max_ambiguous": max_ambiguous}


def check_quality_diversity(fastq_dir, n_samples=3, n_reads=1000):
    """
    Check quality score diversity in FASTQ files to determine if DADA2 is viable.

    DADA2's learnErrors fits a regression of error rate against Q score and
    needs several distinct Q anchors to converge. Heavily binned data (e.g.,
    NovaSeq 2-bin: only Q3 and Q30) produces a degenerate error matrix and
    crashes with "Error matrix is NULL". Truly dummy/placeholder scores (1
    unique Q char) fail the same way.

    Threshold:
      < 5 unique Q chars across scanned reads → "degraded" (use VSEARCH)
      ≥ 5 unique Q chars                       → "normal" (DADA2 viable)

    Strategy: scan reads from all files, accumulating the union of Q chars.
    Early-exit once we cross the EARLY_EXIT threshold (well above the cutoff)
    so normal data finishes in milliseconds. Binned data reads everything but
    saturates quickly and avoids a later DADA2 crash.

    Args:
        fastq_dir: Directory containing FASTQ files
        n_samples: Not used (kept for CLI compatibility)
        n_reads: Not used (kept for CLI compatibility)

    Prints machine-readable lines to stdout:
        QUALITY_STATUS=normal|degraded_binned
        UNIQUE_QUALS=<int>
    """
    DEGRADED_CUTOFF = 5      # < 5 unique Q chars → degraded_binned (compressed/binned Q)
    EARLY_EXIT = 10          # ≥ 10 unique Q chars → certainly normal, stop scanning

    fq_files = sorted(glob.glob(os.path.join(fastq_dir, '*.fastq*')))
    if not fq_files:
        print("WARNING: No FASTQ files found, assuming degraded", file=sys.stderr)
        print("QUALITY_STATUS=degraded_binned")
        print("UNIQUE_QUALS=0")
        return "degraded_binned"

    qual_chars = set()
    files_scanned = 0
    reads_scanned = 0
    early_exit = False

    for fq in fq_files:
        files_scanned += 1
        open_fn = gzip.open if fq.endswith('.gz') else open
        with open_fn(fq, 'rt') as fh:
            for i, line in enumerate(fh):
                if i % 4 == 3:                    # qual line
                    qual_chars.update(line.rstrip())
                    reads_scanned += 1
                    if len(qual_chars) >= EARLY_EXIT:
                        early_exit = True
                        break
        if early_exit:
            break

    n_unique = len(qual_chars)

    if n_unique == 0:
        print("WARNING: No quality data found, assuming degraded", file=sys.stderr)
        print("QUALITY_STATUS=degraded_binned")
        print("UNIQUE_QUALS=0")
        return "degraded_binned"

    status = "degraded_binned" if n_unique < DEGRADED_CUTOFF else "normal"

    scope = "early-exit" if early_exit else "full scan"
    print(f"  Quality check: {n_unique} unique Q chars across "
          f"{files_scanned} file(s), {reads_scanned} reads ({scope})",
          file=sys.stderr)
    print(f"  Unique Q chars: {sorted(qual_chars)}", file=sys.stderr)
    print(f"  Status: {status} (cutoff: <{DEGRADED_CUTOFF})", file=sys.stderr)

    print(f"QUALITY_STATUS={status}")
    print(f"UNIQUE_QUALS={n_unique}")
    return status


def _process_single_fastq(args):
    """Worker function for parallel FASTQ preprocessing (must be top-level for pickling)."""
    fq, output_dir, trim_front, truncate_length, max_n, min_length, is_paired = args
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
            if len(seq_after_front) < max(truncate_length, min_length):
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


def sanitize_fastq(input_dir, min_length=50, sequence_type="single"):
    """
    Remove FASTQ records with actual sequence shorter than min_length.

    For PE data, reads R1/R2 in lockstep and drops both if either is too short.
    Files are overwritten in-place. Handles .fastq and .fastq.gz.

    Args:
        input_dir: Directory containing FASTQ files
        min_length: Minimum sequence length to keep (default: 50)
        sequence_type: "paired" or "single"
    """
    total_in = 0
    total_removed = 0

    if sequence_type == "paired":
        # Find R1/R2 pairs
        all_files = sorted(glob.glob(os.path.join(input_dir, '*.fastq*')))
        r1_files = [f for f in all_files if '_1.fastq' in f or '_R1' in f]

        for r1 in r1_files:
            r2 = r1.replace('_1.fastq', '_2.fastq').replace('_R1', '_R2')
            if not os.path.exists(r2):
                print(f"  WARNING: No R2 found for {os.path.basename(r1)}, skipping pair",
                      file=sys.stderr)
                continue

            open_r1 = gzip.open if r1.endswith('.gz') else open
            open_r2 = gzip.open if r2.endswith('.gz') else open

            kept_r1 = []
            kept_r2 = []
            file_in = 0
            file_removed = 0

            with open_r1(r1, 'rt') as f1, open_r2(r2, 'rt') as f2:
                while True:
                    h1, s1, p1, q1 = f1.readline(), f1.readline(), f1.readline(), f1.readline()
                    h2, s2, p2, q2 = f2.readline(), f2.readline(), f2.readline(), f2.readline()
                    if not h1 or not h2:
                        break
                    file_in += 1
                    seq1 = s1.rstrip('\n')
                    seq2 = s2.rstrip('\n')
                    if len(seq1) < min_length or len(seq2) < min_length:
                        file_removed += 1
                        continue
                    kept_r1.extend([h1, s1, p1, q1])
                    kept_r2.extend([h2, s2, p2, q2])

            # Write back
            write_r1 = gzip.open if r1.endswith('.gz') else open
            write_r2 = gzip.open if r2.endswith('.gz') else open
            with write_r1(r1, 'wt') as f1:
                f1.writelines(kept_r1)
            with write_r2(r2, 'wt') as f2:
                f2.writelines(kept_r2)

            total_in += file_in
            total_removed += file_removed
            if file_removed > 0:
                print(f"  {os.path.basename(r1)}: {file_in} pairs, removed {file_removed}",
                      file=sys.stderr)
    else:
        fq_files = sorted(glob.glob(os.path.join(input_dir, '*.fastq*')))
        for fq in fq_files:
            open_fn = gzip.open if fq.endswith('.gz') else open
            kept = []
            file_in = 0
            file_removed = 0

            with open_fn(fq, 'rt') as fin:
                while True:
                    h, s, p, q = fin.readline(), fin.readline(), fin.readline(), fin.readline()
                    if not h:
                        break
                    file_in += 1
                    if len(s.rstrip('\n')) < min_length:
                        file_removed += 1
                        continue
                    kept.extend([h, s, p, q])

            write_fn = gzip.open if fq.endswith('.gz') else open
            with write_fn(fq, 'wt') as fout:
                fout.writelines(kept)

            total_in += file_in
            total_removed += file_removed
            if file_removed > 0:
                print(f"  {os.path.basename(fq)}: {file_in} reads, removed {file_removed}",
                      file=sys.stderr)

    total_out = total_in - total_removed
    print(f"  Sanitize summary: {total_in} in, {total_out} kept, {total_removed} removed (<{min_length}bp)",
          file=sys.stderr)
    print(f"TOTAL_IN={total_in}")
    print(f"TOTAL_OUT={total_out}")
    print(f"REMOVED={total_removed}")


def degraded_quality_preprocess(input_dir, output_dir, trim_front=15, truncate_length=0,
                                max_n=1, min_length=50, sequence_type="single", threads=4):
    """
    Preprocess FASTQ files for the degraded quality score pipeline.

    When quality scores are binned/unreliable (DADA2 cannot learn error model),
    this function applies simple deterministic preprocessing:
      1. Trim first trim_front bp from each read
      2. Truncate to a fixed length (all output reads are identical length)
      3. Discard reads with N count > max_n
      4. Discard reads with sequence < min_length after trimming
      5. For PE data: use only forward (R1) reads to avoid merge Q value issues

    Files are processed in parallel using multiple workers.

    Args:
        input_dir: Directory with FASTQ files (after adapter removal + primer trimming)
        output_dir: Output directory for preprocessed files
        trim_front: Bases to trim from 5' end (default: 15)
        truncate_length: Fixed length to keep after trim_front (0=auto-detect
                         using 10th percentile of read lengths)
        max_n: Maximum N bases allowed per read; reads exceeding this are discarded (default: 1)
        min_length: Minimum sequence length after trimming; shorter reads are discarded (default: 50)
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

    is_paired = (sequence_type == "paired")
    worker_args = [
        (fq, output_dir, trim_front, truncate_length, max_n, min_length, is_paired)
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

        n_unique = 0
        with open(output_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    n_unique += 1

    print(f"  Dereplicated: {total_reads} reads → {n_unique} unique sequences", file=sys.stderr)
    print(f"DEREP_INPUT_READS={total_reads}")
    print(f"DEREP_UNIQUE_SEQS={n_unique}")


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
            fout.write(f">{sample_id}_{read_num};sample={sample_id};\n{seq}\n")
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
        worker_args = []
        tmp_paths = []
        for idx, row in manifest.iterrows():
            sample_id = str(row['sample-id'])
            filepath = str(row['absolute-filepath'])
            tmp_path = os.path.join(tmpdir, f"sample_{idx}.fasta")
            tmp_paths.append(tmp_path)
            worker_args.append((sample_id, filepath, tmp_path))

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
                    header = line.rstrip('\n')
                    seq_id = header[1:].split(';')[0]
                    fasta_ids.add(seq_id)
                    fout.write(f">{seq_id}\n")
                    n_features += 1
                else:
                    fout.write(line.upper())

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
                if header_line is None:
                    # vsearch --otutabout header starts with '#OTU ID'
                    if line.startswith('#OTU ID') or line.startswith('#OTU\t'):
                        header_line = line.rstrip('\n').split('\t')
                        table_sample_ids = set(header_line[1:])
                        continue
                    if line.startswith('#'):
                        continue
                    # Header without '#' prefix
                    header_line = line.rstrip('\n').split('\t')
                    table_sample_ids = set(header_line[1:])
                    continue
                if line.startswith('#'):
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
        ], check=True)

        subprocess.run([
            'qiime', 'tools', 'import',
            '--input-path', clean_fasta,
            '--type', 'FeatureData[Sequence]',
            '--input-format', 'DNAFASTAFormat',
            '--output-path', output_repseq_qza
        ], check=True)

    print(f"  Imported to QIIME2: {output_table_qza}", file=sys.stderr)
    print(f"  Imported to QIIME2: {output_repseq_qza}", file=sys.stderr)
    print(f"IMPORT_FEATURES={len(table_feature_ids)}")
    print(f"IMPORT_SAMPLES={len(table_sample_ids)}")
    print(f"IMPORT_TOTAL_READS={total_reads}")


def _upsert_csv(output_csv, new_rows, key_cols, fieldnames=None, sep=','):
    """Upsert rows into a CSV keyed by key_cols (idempotent across re-runs).

    Append rows whose key is new, skip rows identical to what is already stored,
    replace rows whose key exists but whose values changed. The whole file is
    rewritten atomically (temp file + os.replace) under an exclusive flock, so it
    is safe when several datasets run in parallel and safe against partial writes.
    Reusable for any per-sample / per-dataset summary table.
    """
    import csv
    import fcntl
    import tempfile

    if not new_rows:
        return
    if fieldnames is None:
        fieldnames = list(new_rows[0].keys())

    def _norm(r):
        def _cell(v):
            # None and NaN both render as empty, matching the old pandas.to_csv output
            if v is None or (isinstance(v, float) and v != v):
                return ''
            return str(v)
        return {fn: _cell(r.get(fn)) for fn in fieldnames}

    def _key(r):
        return tuple(r[k] for k in key_cols)

    lock_path = output_csv + '.lock'
    with open(lock_path, 'w') as lockf:
        fcntl.flock(lockf, fcntl.LOCK_EX)
        try:
            existing = []
            index = {}
            if os.path.exists(output_csv) and os.path.getsize(output_csv) > 0:
                with open(output_csv, newline='') as f:
                    for row in csv.DictReader(f, delimiter=sep):
                        nr = _norm(row)
                        index[_key(nr)] = len(existing)
                        existing.append(nr)

            for raw in new_rows:
                nr = _norm(raw)
                k = _key(nr)
                if k in index:
                    if existing[index[k]] != nr:
                        existing[index[k]] = nr
                else:
                    index[k] = len(existing)
                    existing.append(nr)

            out_dir = os.path.dirname(os.path.abspath(output_csv))
            fd, tmp = tempfile.mkstemp(dir=out_dir, suffix='.tmp')
            with os.fdopen(fd, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames,
                                        delimiter=sep, extrasaction='ignore',
                                        lineterminator='\n')
                writer.writeheader()
                writer.writerows(existing)
            os.replace(tmp, output_csv)
        finally:
            fcntl.flock(lockf, fcntl.LOCK_UN)


def append_summary(dataset_id, sra_file, raw_counts_file, final_table, output_csv, sequence_type="single"):
    """
    Append per-sample summary (raw reads + final reads) to a unified CSV.

    Args:
        dataset_id: BioProject ID
        sra_file: Path to _sra.txt (Run<tab>SampleName)
        raw_counts_file: Path to _raw_read_counts.tsv (Run<tab>SampleName<tab>RawReads)
        final_table: Path to final-table.qza
        output_csv: Path to the unified summary CSV
        sequence_type: Original sequence type ('paired' or 'single')
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
        # For PE data, raw_read_counts.tsv has R1+R2 total; divide by 2 for per-sample pair count
        if sequence_type == "paired" and raw_read_count > 0:
            raw_read_count = raw_read_count // 2

        final_read_count = final_reads.get(sample_name, 0)

        rows.append({
            'BioProject': dataset_id,
            'Run': run_id,
            'SampleName': sample_name,
            'RawReads': raw_read_count,
            'FinalReads': final_read_count
        })

    result_df = pd.DataFrame(rows)

    # 5. Upsert into the unified CSV, keyed by Run (idempotent across re-runs):
    #    append new rows, skip identical ones, replace changed ones. Parallel-safe
    #    and crash-safe (flock + atomic temp/rename) inside _upsert_csv. This keeps
    #    previously-succeeded datasets' rows intact when the pipeline is re-run.
    fieldnames = ['BioProject', 'Run', 'SampleName', 'RawReads', 'FinalReads']
    _upsert_csv(output_csv, result_df.to_dict('records'),
                key_cols=['Run'], fieldnames=fieldnames)

    print(f"✓ Summary upserted for {dataset_id}: {len(rows)} samples")


# ===========================================================================
# Amplified-region detection (align rep-seqs to E. coli 16S, map to V-regions)
# ===========================================================================
# V-region boundaries in E. coli 16S numbering (reference J01859), following
# Klindworth et al. 2013 (Nucleic Acids Res) and common usage. Tunable here.
VREGIONS = [("V1", 69, 99), ("V2", 137, 242), ("V3", 433, 497), ("V4", 576, 682),
            ("V5", 822, 879), ("V6", 986, 1043), ("V7", 1117, 1173),
            ("V8", 1243, 1294), ("V9", 1435, 1465)]


def _call_region(lo, hi, min_overlap=0.5):
    """Map an amplicon's E. coli span [lo, hi] to the V-region(s) it covers.
    A V-region counts as covered when the span overlaps >= min_overlap of its width."""
    covered = [name for (name, vs, ve) in VREGIONS
               if max(0, min(hi, ve) - max(lo, vs)) >= min_overlap * (ve - vs)]
    if not covered:
        return "unclassified"
    if len(covered) == len(VREGIONS):
        return "V1-V9"  # full-length
    return covered[0] if len(covered) == 1 else f"{covered[0]}-{covered[-1]}"


def detect_region(repseqs_qza, ecoli_ref, threads=4):
    """Infer the amplified 16S V-region by globally aligning representative
    sequences to the E. coli 16S reference (both strands) and taking the median
    target span. The variable middle aligns poorly but the conserved primer-flanking
    ends anchor the span, so the [start, end] is robust. Returns a dict with region,
    E. coli start/end, confidence (fraction of rep-seqs agreeing with the consensus
    call) and counts."""
    import glob
    import statistics
    import subprocess
    import tempfile

    blank = {"region": "NA", "ecoli_start": "", "ecoli_end": "",
             "confidence": "0.00", "n_repseqs": 0, "n_aligned": 0}
    if not (repseqs_qza and os.path.exists(repseqs_qza) and os.path.exists(ecoli_ref)):
        return blank
    with tempfile.TemporaryDirectory() as td:
        try:
            subprocess.run(["qiime", "tools", "export", "--input-path", repseqs_qza,
                            "--output-path", td], check=True, capture_output=True)
        except Exception:
            return blank
        fastas = glob.glob(os.path.join(td, "*.fasta"))
        if not fastas:
            return blank
        fa = fastas[0]
        n_rep = sum(1 for line in open(fa) if line.startswith(">"))
        aln = os.path.join(td, "region_aln.tsv")
        subprocess.run(["vsearch", "--usearch_global", fa, "--db", ecoli_ref,
                        "--id", "0.5", "--strand", "both", "--query_cov", "0.80",
                        "--userout", aln, "--userfields", "query+tilo+tihi+id",
                        "--maxaccepts", "1", "--top_hits_only",
                        "--threads", str(threads), "--quiet"],
                       check=False, capture_output=True)
        rows = []
        if os.path.exists(aln):
            for line in open(aln):
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3:
                    try:
                        rows.append((int(p[1]), int(p[2])))
                    except ValueError:
                        pass
        if not rows:
            return {**blank, "region": "unclassified", "n_repseqs": n_rep}
        mlo = int(statistics.median([r[0] for r in rows]))
        mhi = int(statistics.median([r[1] for r in rows]))
        region = _call_region(mlo, mhi)
        same = sum(1 for lo, hi in rows if _call_region(lo, hi) == region)
        return {"region": region, "ecoli_start": mlo, "ecoli_end": mhi,
                "confidence": f"{same / len(rows):.2f}",
                "n_repseqs": n_rep, "n_aligned": len(rows)}


def build_per_dataset_summary(output_dir, mode, ecoli_ref, summary_csv=None, threads=4):
    """After AmpliconPIP: write one row per SUCCESSFUL dataset (those that produced
    final rep-seqs) combining platform, quality status and the detected amplified
    region. Upserts into per_dataset_summary.tsv (keyed by Bioproject) so re-runs
    never drop previously-summarised datasets. Region results are cached per dataset
    and only recomputed when the rep-seqs are newer than the cache."""
    import csv
    import glob

    out_tsv = os.path.join(output_dir, "per_dataset_summary.tsv")

    platform_map = {}
    pc = os.path.join(output_dir, ".platform_cache.txt")
    if os.path.exists(pc):
        for line in open(pc):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                platform_map[p[0].strip()] = p[1].strip()

    reads = {}
    if summary_csv and os.path.exists(summary_csv):
        for r in csv.DictReader(open(summary_csv)):
            bp = r.get("BioProject", "")
            d = reads.setdefault(bp, {"n": 0, "raw": 0, "final": 0})
            d["n"] += 1
            for src, key in (("RawReads", "raw"), ("FinalReads", "final")):
                try:
                    d[key] += int(r.get(src) or 0)
                except ValueError:
                    pass

    rows = []
    for q in sorted(glob.glob(os.path.join(output_dir, "*", f"*-{mode}-final-rep-seqs.qza"))):
        ds = os.path.basename(os.path.dirname(q))
        rcache = os.path.join(output_dir, ds, f"{ds}_region.tsv")
        keys = ("region", "ecoli_start", "ecoli_end", "confidence", "n_repseqs", "n_aligned")
        reg = None
        if os.path.exists(rcache) and os.path.getmtime(rcache) >= os.path.getmtime(q):
            cached = open(rcache).read().rstrip("\n").split("\t")
            if len(cached) == len(keys):
                reg = dict(zip(keys, cached))          # reuse a valid cache
        if reg is None:
            reg = detect_region(q, ecoli_ref, threads=threads)
            # Cache real results only; never pin a transient failure (region == "NA"),
            # so a one-off hiccup is retried next run. A malformed cache also falls here
            # and gets rewritten -> self-healing.
            if reg.get("region") != "NA":
                with open(rcache, "w") as f:
                    f.write("\t".join(str(reg[k]) for k in keys) + "\n")
        qs = "NA"
        qf = os.path.join(output_dir, ds, f"{ds}_quality_status.txt")
        if os.path.exists(qf):
            qs = (open(qf).read().strip() or "NA")
        plat = "NA"
        pf = os.path.join(output_dir, ds, f"{ds}_platform.txt")
        if os.path.exists(pf):
            plat = (open(pf).read().strip() or "NA")
        elif ds in platform_map:
            plat = platform_map[ds]
        st = reads.get(ds, {})
        rows.append({
            "Bioproject": ds,
            "Platform": plat,
            "QualityStatus": qs,
            "Region": reg["region"],
            "RegionEcoliStart": reg["ecoli_start"],
            "RegionEcoliEnd": reg["ecoli_end"],
            "RegionConfidence": reg["confidence"],
            "NRepSeqs": reg["n_repseqs"],
            "NSamples": st.get("n", ""),
            "RawReads": st.get("raw", ""),
            "FinalReads": st.get("final", ""),
        })

    fields = ["Bioproject", "Platform", "QualityStatus", "Region",
              "RegionEcoliStart", "RegionEcoliEnd", "RegionConfidence",
              "NRepSeqs", "NSamples", "RawReads", "FinalReads"]
    _upsert_csv(out_tsv, rows, key_cols=["Bioproject"], fieldnames=fields, sep="\t")
    print(f"\u2713 per_dataset_summary.tsv: {len(rows)} dataset(s) -> {out_tsv}")


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
            "batch_get_sequencing_platforms",
            "adaptive_tail_trim",
            "detect_length_window",
            "check_quality_diversity",
            "sanitize_fastq",
            "degraded_quality_preprocess",
            "derep_fastq_for_vsearch",
            "relabel_reads_for_mapping",
            "import_vsearch_to_qiime2",
            "append_summary",
            "detect_region",
            "build_per_dataset_summary"
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
    parser.add_argument("--pairs_file", help="TSV file with dataset_id<TAB>srr_id[<TAB>bioproject_id] per line")
    parser.add_argument("--input_dir", help="Input directory (for adaptive_tail_trim)")
    parser.add_argument("--output_dir", help="Output directory (for adaptive_tail_trim)")
    parser.add_argument("--max_sample_reads", type=int, default=10000,
                        help="Max reads to sample per file (default: 10000, 0=all)")
    parser.add_argument("--tolerance", type=float, default=0.15,
                        help="Length-window half-width fraction around peak (default: 0.15)")
    parser.add_argument("--floor", type=int, default=200,
                        help="Hard minimum for the lower length bound in bp (default: 200)")
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
    parser.add_argument("--truncate_length", type=int, default=0,
                        help="Fixed truncation length after trim_front (0=auto-detect)")
    parser.add_argument("--max_n", type=int, default=1,
                        help="Max N bases allowed per read (default: 1)")
    parser.add_argument("--min_length", type=int, default=50,
                        help="Min sequence length to keep (default: 50)")
    parser.add_argument("--sequence_type", default="single",
                        help="Sequence type: 'paired' or 'single' (default: single)")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads for parallel operations (default: 4)")
    parser.add_argument("--output_fasta", help="Output FASTA path")
    parser.add_argument("--manifest_path", help="Path to QIIME2 manifest TSV")
    parser.add_argument("--zotu_fasta", help="Path to ZOTU FASTA (for import_vsearch_to_qiime2)")
    parser.add_argument("--otu_table_tsv", help="Path to vsearch OTU table TSV")
    parser.add_argument("--output_table_qza", help="Output table .qza path")
    parser.add_argument("--output_repseq_qza", help="Output rep-seqs .qza path")
    parser.add_argument("--repseqs", help="Path to rep-seqs .qza (for detect_region)")
    parser.add_argument("--ecoli_ref", help="Path to E. coli 16S reference FASTA (region detection)")
    parser.add_argument("--mode", help="Denoising mode token (asv|otu) used in output filenames")

    args = parser.parse_args()
    
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

    elif args.function == "batch_get_sequencing_platforms":
        batch_get_sequencing_platforms(args.pairs_file)

    elif args.function == "adaptive_tail_trim":
        adaptive_tail_trim(args.input_dir, args.output_dir, args.max_sample_reads)

    elif args.function == "detect_length_window":
        detect_length_window(args.input_dir, args.tolerance, args.floor,
                             args.max_sample_reads)

    elif args.function == "check_quality_diversity":
        check_quality_diversity(args.input_dir, args.n_samples, args.n_reads)

    elif args.function == "sanitize_fastq":
        sanitize_fastq(args.input_dir, args.min_length, args.sequence_type)

    elif args.function == "degraded_quality_preprocess":
        degraded_quality_preprocess(args.input_dir, args.output_dir,
                                    args.trim_front, args.truncate_length,
                                    args.max_n, args.min_length,
                                    args.sequence_type, args.threads)

    elif args.function == "derep_fastq_for_vsearch":
        derep_fastq_for_vsearch(args.input_dir, args.output_fasta, args.threads)

    elif args.function == "relabel_reads_for_mapping":
        relabel_reads_for_mapping(args.manifest_path, args.output_fasta, args.threads)

    elif args.function == "import_vsearch_to_qiime2":
        import_vsearch_to_qiime2(args.zotu_fasta, args.otu_table_tsv,
                                  args.manifest_path,
                                  args.output_table_qza, args.output_repseq_qza)

    elif args.function == "append_summary":
        append_summary(args.dataset_id, args.sra_file, args.raw_counts, args.final_table, args.output_csv, args.sequence_type)
    elif args.function == "detect_region":
        res = detect_region(args.repseqs, args.ecoli_ref, threads=args.threads)
        print(f"REGION={res['region']}")
        print(f"REGION_ECOLI_START={res['ecoli_start']}")
        print(f"REGION_ECOLI_END={res['ecoli_end']}")
        print(f"REGION_CONFIDENCE={res['confidence']}")
        print(f"REGION_N_ALIGNED={res['n_aligned']}/{res['n_repseqs']}")
    elif args.function == "build_per_dataset_summary":
        build_per_dataset_summary(args.output_dir, args.mode, args.ecoli_ref,
                                  summary_csv=args.output_csv, threads=args.threads)