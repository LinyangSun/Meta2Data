#!/usr/bin/env python3

import sys
import subprocess
import logging
import os
import pandas as pd
import numpy as np

# ----------------------------------------------------------------------
# ✅ Logging Configuration (consistent for entire script)
# ----------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stderr,
    force=True
)

# ----------------------------------------------------------------------
# ✅ Module Auto-Installer
# ----------------------------------------------------------------------
def check_and_install(module, module2):
    try:
        __import__(module)
    except ImportError:
        logging.warning(f"Module '{module}' not found. Installing '{module2}'...")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", module2],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        __import__(module)


# ----------------------------------------------------------------------
# ✅ Generate Datasets ID File (using Bioproject as dataset ID)
# ----------------------------------------------------------------------
def GenerateDatasetsIDsFile(file_path, Bioproject, Data_SequencingPlatform):
    directory_path = os.path.dirname(os.path.abspath(file_path))
    try:
        df = pd.read_csv(file_path)
        logging.info(f"CSV file loaded: {file_path}")
    except Exception as e:
        logging.error(f"Error reading the CSV file: {e}")
        sys.exit(1)

    columns_to_clean = [Bioproject, Data_SequencingPlatform]
    for col in columns_to_clean:
        if col in df.columns:
            df[col] = (
                df[col]
                .astype(str)
                .str.replace(' ', '', regex=True)
                .str.replace('\n', '', regex=True)
                .str.replace('\t', '', regex=True)
            )
        else:
            logging.warning(f"Column '{col}' not found in the file.")

    df_pair = (
        df[[Bioproject, Data_SequencingPlatform]]
        .dropna(subset=[Bioproject])
        .drop_duplicates()
    )

    datasets = df_pair.values.astype(object)
    out_path = f"{directory_path}/datasets_ID.txt"
    np.savetxt(out_path, datasets, fmt="%s")
    logging.info(f"datasets_ID.txt created at: {out_path}")
    return datasets


# ----------------------------------------------------------------------
# ✅ Generate SRA Files (using Bioproject as dataset ID)
# ----------------------------------------------------------------------
def GenerateSRAsFile(file_path, Bioproject, SRA_Number, Biosample, Forward, Reverse, region):
    directory_path = os.path.dirname(os.path.abspath(file_path))
    try:
        df = pd.read_csv(file_path)
        logging.info(f"CSV file loaded: {file_path}")
    except Exception as e:
        logging.error(f"Error reading the CSV file: {e}")
        sys.exit(1)

    columns_to_clean = [Bioproject, SRA_Number, Biosample, Forward, Reverse, region]
    for col in columns_to_clean:
        if col in df.columns:
            df[col] = (
                df[col]
                .astype(str)
                .str.replace(' ', '', regex=True)
                .str.replace('\n', '', regex=True)
                .str.replace('\t', '', regex=True)
            )
        else:
            logging.warning(f"Column '{col}' not found in the file.")

    df["rename"] = df[Bioproject] + '-' + df[Biosample]
    datasets = np.array(
        [str(x).strip().replace('\t', '') for x in df[Bioproject].dropna().unique()],
        dtype=object
    )

    for value in datasets:
        sub = df[df[Bioproject] == value].copy()
        out_df = pd.concat([
            sub[SRA_Number],
            sub["rename"],
            sub[region],
            sub[Forward],
            sub[Reverse]
        ], axis=1)

        folder_path = f"{directory_path}/{value}/"
        os.makedirs(folder_path, exist_ok=True)
        out_file = f'{folder_path}/{value}_sra.txt'
        out_df.to_csv(out_file, sep='\t', header=None, index=None)
        logging.info(f"SRA file created: {out_file}")

    logging.info("All SRA files finished.")


# ----------------------------------------------------------------------
# ✅ Manifest Generators
# ----------------------------------------------------------------------
def mk_manifest_SE(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t', header=None)
    except Exception as e:
        logging.error(f"Failed to read the file: {e}")
        sys.exit(1)

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
    logging.info(f"Manifest file created: {out_path}")


def mk_manifest_PE(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t', header=None)
    except Exception as e:
        logging.error(f"Failed to read the file: {e}")
        sys.exit(1)

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
    logging.info(f"Manifest file created: {out_path}")


# ----------------------------------------------------------------------
# ✅ Combine Stats File
# ----------------------------------------------------------------------
def combine_stats_file(filter_stats, deblur_stats, dataset_name):
    try:
        filter_stats_df = pd.read_csv(filter_stats)
        deblur_stats_df = pd.read_csv(deblur_stats)
    except Exception as e:
        logging.error(f"Error reading one of the stats files: {e}")
        sys.exit(1)

    merged_df = pd.merge(filter_stats_df, deblur_stats_df, on='sample-id', how='outer')
    out_path = os.path.join(os.path.dirname(filter_stats), f"{dataset_name}_process_stats.csv")
    merged_df.to_csv(out_path, index=False)
    logging.info(f"Combined stats file created: {out_path}")


# ----------------------------------------------------------------------
# ✅ Add Prefix to FASTA and BIOM
# ----------------------------------------------------------------------
def add_prefix_to_file(in_fasta, in_table, prefix):
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
    logging.info(f"Prefixed FASTA written to {out_fasta}")

    table = load_table(in_table)
    out_table = os.path.join(dir, f"{prefix}-feature-table.biom")
    observation_ids = table.ids(axis='observation')
    new_observation_ids = [prefix + obs_id for obs_id in observation_ids]
    new_data = table.matrix_data.tocsr()
    new_table = Table(new_data, new_observation_ids, table.ids())
    with h5py.File(out_table, 'w') as f:
        new_table.to_hdf5(f, 'representive sequence')
    logging.info(f"Prefixed BIOM written to {out_table}")


# ----------------------------------------------------------------------
# ✅ Trimming Functions
# ----------------------------------------------------------------------

def trim_pos_deblur(file_path):
    """
    Calculate trim positions for Deblur sequencing data.
    
    Args:
        file_path: Path to forward-seven-number-summaries.tsv file
    
    Returns:
        tuple: (start_position, end_position) or (None, None) if failed
    """
    logging.info("\n" + "="*60)
    logging.info("🔬 Deblur Sequencing Trim Position Calculator")
    logging.info("="*60)
    
    # ------------- read & parse -------------
    with open(file_path, "r") as f:
        tsv = [line.rstrip("\n") for line in f]
    pos = [int(x) for x in tsv[0].split("\t")[1:]]
    L = len(pos)
    counts = [float(x) for x in tsv[1].split("\t")[1:]]
    
    # 1) outer1 := first index where counts < threshold (fallback: last index)
    # Adjust threshold based on actual data range
    threshold = min(9000, max(counts) * 0.7)  # Use 70% of max count if max < 18000
    under = [i for i, v in enumerate(counts) if v < threshold]
    outer1 = under[0] if under else (L - 1)
    
    # Ensure outer1 is at least at position 40 to have valid inner/outer regions
    outer1 = max(outer1, 40)
    
    logging.info(f"\n📊 File diagnostics:")
    logging.info(f"   Total positions: {L}")
    logging.info(f"   Count threshold used: {threshold:.1f}")
    logging.info(f"   Outer1 threshold position: {outer1}")
    logging.info(f"   Count range: {min(counts):.1f} - {max(counts):.1f}")
    
    # 5) regions
    inner_region = range(0, 20)          # [0, 20), positions 0..19
    outer_region = range(max(0, outer1 - 20), outer1)  # [outer1-20, outer1)
    
    logging.info(f"   Inner region: positions 0-19")
    logging.info(f"   Outer region: positions {max(0, outer1 - 20)}-{outer1-1}")
    
    # rows after the counts: skip 2% by slicing from index 3
    rows = []
    for line in tsv[3:]:
        vals = [float(x) for x in line.split("\t")[1:][:L]]
        rows.append(vals)
    
    per_row = []
    for row_idx, vals in enumerate(rows):
        # 3–4) NA where value < 20
        na_idx = {i for i, v in enumerate(vals) if v < 20}
        
        logging.info(f"\n--- Row {row_idx} ---")
        logging.info(f"  Total positions with value < 20: {len(na_idx)}")
        
        # ---- LEFT: start (inner region) ----
        inner_na = sorted(i for i in na_idx if i in inner_region)
        n_inner = len(inner_na)
        
        logging.info(f"  Inner region (0-19): {n_inner} positions < 20")
        if n_inner > 0:
            logging.info(f"    Bad positions: {inner_na}")
        
        if n_inner == 0:
            start = 0
            logging.info(f"    ✓ start = {start} (no bad positions)")
        elif n_inner >= 20:
            # All positions bad - reject row
            logging.info(f"    ✗ REJECTED - all 20 inner region positions have values < 20")
            start = None
        elif n_inner >= 15:
            # Too many bad positions (≥75%) - reject row
            logging.info(f"    ✗ REJECTED - {n_inner}/20 inner region positions have values < 20 (≥75%)")
            start = None
        else:
            # Find contiguous run of NAs from position 0
            s = set(inner_na)
            r = 0
            while r in s and r < 20:
                r += 1
            if r > 0:
                start = r  # Position after the contiguous run
                logging.info(f"    ✓ start = {start} (after contiguous run of {r} bad positions from pos 0)")
            else:
                start = 0
                logging.info(f"    ✓ start = {start} (bad positions not contiguous from pos 0)")
        
        # ---- RIGHT: end (outer region) ----
        outer_na = sorted(i for i in na_idx if i in outer_region)
        n_outer = len(outer_na)
        
        logging.info(f"  Outer region ({max(0, outer1 - 20)}-{outer1-1}): {n_outer} positions < 20")
        if n_outer > 0:
            logging.info(f"    Bad positions: {outer_na}")
        
        if n_outer == 0:
            end = outer1
            logging.info(f"    ✓ end = {end} (no bad positions)")
        elif n_outer >= 20:
            # All positions bad - reject row
            logging.info(f"    ✗ REJECTED - all 20 outer region positions have values < 20")
            end = None
        elif n_outer >= 15:
            # Too many bad positions (≥75%) - reject row
            logging.info(f"    ✗ REJECTED - {n_outer}/20 outer region positions have values < 20 (≥75%)")
            end = None
        else:
            # Find contiguous run of NAs ending at outer1
            s = set(outer_na)
            r = 0
            while (outer1 - r) in s:
                r += 1
            if r > 0:
                # Contiguous NA run ending at outer1
                end = outer1 - r - 1
                logging.info(f"    ✓ end = {end} (before contiguous run of {r} bad positions ending at {outer1})")
            else:
                end = outer1
                logging.info(f"    ✓ end = {end} (bad positions not contiguous at end)")
        
        logging.info(f"  Result: start={start}, end={end}")
        per_row.append((start, end))
    
    # 7) combine rows
    starts = [s for s, e in per_row if s is not None]
    ends = [e for s, e in per_row if e is not None]
    
    logging.info(f"\n📈 Per-row results:")
    logging.info(f"   Total rows analyzed: {len(rows)}")
    logging.info(f"   Valid start positions: {starts if starts else 'None'}")
    logging.info(f"   Valid end positions: {ends if ends else 'None'}")
    
    # Check if we have enough valid rows
    if not starts or not ends:
        logging.warning("\n❌ FAILED: All rows were rejected - no valid trim points available")
        return (None, None)
    
    # Require at least 50% of rows to be valid
    valid_rows = sum(1 for s, e in per_row if s is not None and e is not None)
    if valid_rows < len(rows) * 0.5:
        logging.warning(f"\n❌ FAILED: Only {valid_rows}/{len(rows)} rows are valid (need ≥50%)")
        return (None, None)
    
    # Take most conservative trim points (no additional +1/-1 adjustment)
    final_start = max(starts)
    final_end = min(ends)
    
    logging.info(f"\n🎯 Trim calculation:")
    logging.info(f"   Most conservative start: {final_start} (max of valid starts)")
    logging.info(f"   Most conservative end: {final_end} (min of valid ends)")
    
    # Clamp to valid range
    final_start = max(0, min(final_start, L - 1))+1
    final_end = max(0, min(final_end, L - 1))-1
    
    logging.info(f"   After clamping: start={final_start}, end={final_end}")
    
    # Sanity check
    if final_start >= final_end:
        logging.warning(f"\n❌ FAILED: final_start ({final_start}) >= final_end ({final_end}) - no valid range")
        return (None, None)
    
    logging.info(f"\n✅ Everything good! Valid rows: {valid_rows}/{len(rows)}")
    logging.info(f"   Trim range: positions {final_start} to {final_end} (length: {final_end - final_start})")
    print(f"{final_start},{final_end}")
    return (final_start, final_end)


def trim_pos_454(file_path):
    """
    Calculate trim length for 454 sequencing data.
    Returns the position where read count drops below 9000.
    
    Args:
        file_path: Path to forward-seven-number-summaries.tsv file
    
    Returns:
        tuple: (trim_length,) - single value for 454 trimming
    """
    logging.info("\n" + "="*60)
    logging.info("🔬 454 Sequencing Trim length Calculator")
    logging.info("="*60)
    
    # ------------- Read & Parse TSV File -------------
    logging.info(f"\n📂 Reading file: {file_path}")
    
    try:
        with open(file_path, "r") as f:
            tsv = [line.rstrip("\n") for line in f]
    except FileNotFoundError:
        logging.error(f"❌ ERROR: File not found: {file_path}")
        return (None,)
    except Exception as e:
        logging.error(f"❌ ERROR: Failed to read file: {e}")
        return (None,)
    
    # Parse positions and counts from first two rows
    try:
        pos = [int(x) for x in tsv[0].split("\t")[1:]]
        counts = [float(x) for x in tsv[1].split("\t")[1:]]
    except (IndexError, ValueError) as e:
        logging.error(f"❌ ERROR: Failed to parse TSV data: {e}")
        return (None,)
    
    L = len(pos)
    
    logging.info(f"✓ Successfully parsed TSV file")
    logging.info(f"   Total positions: {L}")
    logging.info(f"   Position range: {pos[0]} to {pos[-1]}")
    
    # ------------- Calculate Statistics -------------
    logging.info(f"\n📊 Read Count Statistics:")
    logging.info(f"   Minimum count: {min(counts):,.1f}")
    logging.info(f"   Maximum count: {max(counts):,.1f}")
    logging.info(f"   Average count: {sum(counts)/len(counts):,.1f}")
    
    # ------------- Find Trim Position -------------
    threshold = 9000
    logging.info(f"\n🎯 Trim Position Calculation:")
    logging.info(f"   Threshold: {threshold:,} reads")
    
    # Find first position where count drops below threshold
    under_threshold = [i for i, v in enumerate(counts) if v < threshold]
    
    if not under_threshold:
        # All positions have >= 9000 reads
        trim_length = L
        logging.info(f"   ✓ All positions have ≥{threshold:,} reads")
        logging.info(f"   Using maximum length: {trim_length}")
    else:
        trim_length = under_threshold[0]
        logging.info(f"   ✓ Found {len(under_threshold)} positions with <{threshold:,} reads")
        logging.info(f"   First drop-off at position index: {trim_length}")
        logging.info(f"   Read count at position {trim_length}: {counts[trim_length]:,.1f}")
        if trim_length > 0:
            logging.info(f"   Read count at position {trim_length-1}: {counts[trim_length-1]:,.1f}")
    
    # ------------- Validation -------------
    logging.info(f"\n✅ Validation:")
    
    if trim_length == 0:
        logging.warning(f"   ⚠️  WARNING: Trim length is 0 - all positions below threshold!")
        logging.warning(f"   This indicates very low quality data")
        return (None,)
    
    if trim_length < 100:
        logging.warning(f"   ⚠️  WARNING: Trim length ({trim_length}) is very short (<100bp)")
        logging.warning(f"   Consider checking data quality")
    
    logging.info(f"   Trim length: {trim_length} bp")
    logging.info(f"   Retained: {trim_length}/{L} positions ({trim_length/L*100:.1f}%)")
    
    # ------------- Summary -------------
    logging.info(f"\n" + "="*60)
    logging.info(f"📋 FINAL RESULT:")
    logging.info(f"   Trim length for 454 data: {trim_length} bp")
    logging.info("="*60 + "\n")
    
    # Print result to stdout for shell capture
    print(f"{trim_length}")
    
    return (trim_length,)


def trim_pos_dada2(file_path):
    """
    Calculate trim positions for DADA2 sequencing data.
    Placeholder function - implement based on your DADA2 requirements.
    
    Args:
        file_path: Path to quality profile file
    
    Returns:
        tuple: (trim_left, trim_right) or appropriate values for DADA2
    """
    logging.info("\n" + "="*60)
    logging.info("🔬 DADA2 Trim Position Calculator")
    logging.info("="*60)
    logging.warning("⚠️  This function is not yet implemented!")
    logging.info("   Please implement DADA2-specific trimming logic")
    logging.info("="*60 + "\n")
    
    # TODO: Implement DADA2 trimming logic
    return (None, None)


# ----------------------------------------------------------------------
# ✅ Main Execution
# ----------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    
    # Check and install required modules
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
            "add_prefix_to_file"
        ],
        help="The function to execute."
    )
    parser.add_argument("--FilePath", help="Path to the input file (include file name).")
    parser.add_argument("--SequencingPlatform", help="Sequencing Platform column name (required for GenerateDatasetsIDsFile).")
    parser.add_argument("--Forward", help="Forward column name (required for GenerateSRAsFile).")
    parser.add_argument("--Reverse", help="Reverse column name (required for GenerateSRAsFile).")
    parser.add_argument("--Datasets_ID", help="Datasets_ID column name (required for some functions).")
    parser.add_argument("--Bioproject", help="Bioproject column name (required for GenerateSRAsFile).")
    parser.add_argument("--SRA_Number", help="SRA_Number column name (required for GenerateSRAsFile).")
    parser.add_argument("--Biosample", help="Biosample column name (required for GenerateSRAsFile).")
    parser.add_argument("--region", help="Region column name (required for GenerateSRAsFile).")
    parser.add_argument("--filter_stats", help="Path to the filter stats file (include file name).")
    parser.add_argument("--deblur_stats", help="Path to the deblur stats file (include file name).")
    parser.add_argument("--dataset_name", help="Datasets ID.")
    parser.add_argument("--In_fasta", help="Input fasta file.")
    parser.add_argument("--In_table", help="Input table file.")
    parser.add_argument("--Prefix", help="Prefix for output files.")

    args = parser.parse_args()
    
    # Execute the requested function
    if args.function == "mk_manifest_SE":
        if not args.FilePath:
            parser.error("mk_manifest_SE requires --FilePath")
        mk_manifest_SE(args.FilePath)
        
    elif args.function == "mk_manifest_PE":
        if not args.FilePath:
            parser.error("mk_manifest_PE requires --FilePath")
        mk_manifest_PE(args.FilePath)
        
    elif args.function == "combine_stats_file":
        if not all([args.filter_stats, args.deblur_stats, args.dataset_name]):
            parser.error("combine_stats_file requires --filter_stats, --deblur_stats, --dataset_name")
        combine_stats_file(args.filter_stats, args.deblur_stats, args.dataset_name)
        
    elif args.function == "trim_pos_dada2":
        if not args.FilePath:
            parser.error("trim_pos_dada2 requires --FilePath")
        trim_pos_dada2(args.FilePath)
        
    elif args.function == "trim_pos_454":
        if not args.FilePath:
            parser.error("trim_pos_454 requires --FilePath")
        trim_pos_454(args.FilePath)
        
    elif args.function == "trim_pos_deblur":
        if not args.FilePath:
            parser.error("trim_pos_deblur requires --FilePath")
        trim_pos_deblur(args.FilePath)
        
    elif args.function == "GenerateDatasetsIDsFile":
        if not all([args.FilePath, args.Datasets_ID, args.SequencingPlatform]):
            parser.error("GenerateDatasetsIDsFile requires --FilePath, --Datasets_ID, --SequencingPlatform")
        GenerateDatasetsIDsFile(args.FilePath, args.Datasets_ID, args.SequencingPlatform)
        
    elif args.function == "GenerateSRAsFile":
        if not all([args.FilePath, args.Datasets_ID, args.Bioproject, args.SRA_Number, 
                    args.Biosample, args.Forward, args.Reverse, args.region]):
            parser.error("GenerateSRAsFile requires --FilePath, --Datasets_ID, --Bioproject, "
                        "--SRA_Number, --Biosample, --Forward, --Reverse, --region")
        GenerateSRAsFile(args.FilePath, args.Datasets_ID, args.Bioproject, args.SRA_Number, 
                        args.Biosample, args.Forward, args.Reverse, args.region)
        
    elif args.function == "add_prefix_to_file":
        if not all([args.In_fasta, args.In_table, args.Prefix]):
            parser.error("add_prefix_to_file requires --In_fasta, --In_table, --Prefix")
        add_prefix_to_file(args.In_fasta, args.In_table, args.Prefix)
        
    else:
        logging.error("Unknown function specified.")
        sys.exit(1)