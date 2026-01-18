#!/usr/bin/env python3
"""
Unified Metadata Download Pipeline
Supports both NCBI (PRJNA/PRJEB/PRJDB) and CNCB/GSA (PRJCA) databases

Workflow:
1. Route BioProjects by ID prefix (PRJCA → CNCB, others → NCBI)
2. Download metadata:
   - NCBI: BioSample + SRA RunInfo (via Entrez API)
   - CNCB/GSA: Metadata files via iSeq tool (conda-installed)
3. Merge in stages:
   - Stage 1: Merge NCBI data (SRA + BioSample)
   - Stage 2: Parse CNCB metadata files (supports CSV, Excel, TSV)
   - Stage 3: Merge NCBI and CNCB by common columns

Dependencies:
   - NCBI Entrez API (Biopython)
   - iSeq tool (conda install -c bioconda iseq)

Usage:
    python3 unified_metadata_downloader.py \\
        --input /path/to/bioproject_ids/ \\
        --output /path/to/output/ \\
        --email your@email.com
"""

from Bio import Entrez
import time
import os
import pandas as pd
import re
from pathlib import Path
import glob
import argparse
import sys
import subprocess


# ============================================================================
# BioProject Routing
# ============================================================================

def classify_bioprojects(bioproject_list):
    """
    Classify BioProjects into NCBI and CNCB groups

    Args:
        bioproject_list: List of BioProject IDs

    Returns:
        dict: {'ncbi': [...], 'cncb': [...]}
    """
    ncbi_projects = []
    cncb_projects = []

    for bioproject in bioproject_list:
        bioproject = bioproject.strip()
        if bioproject.startswith(('PRJCA', 'CRA')):
            cncb_projects.append(bioproject)
        elif bioproject.startswith(('PRJNA', 'PRJEB', 'PRJDB')):
            ncbi_projects.append(bioproject)
        else:
            print(f"  Warning: Unknown BioProject format: {bioproject}")

    return {
        'ncbi': ncbi_projects,
        'cncb': cncb_projects
    }


def read_bioproject_ids(folder_path):
    """Read BioProject IDs from all txt files in folder"""
    txt_files = glob.glob(os.path.join(folder_path, "*.txt"))

    if not txt_files:
        print(f"ERROR: No txt files found in {folder_path}")
        return []

    all_ids = []
    for file in txt_files:
        try:
            with open(file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        all_ids.append(line)
        except Exception as e:
            print(f"  Warning: Failed to read {file}: {e}")
            continue

    # Remove duplicates while preserving order
    unique_ids = list(dict.fromkeys(all_ids))
    print(f"Found {len(unique_ids)} unique BioProject IDs from {len(txt_files)} files")

    return unique_ids


# ============================================================================
# NCBI Download Functions (from ncbi_metadata_downloader.py)
# ============================================================================

def get_biosamples_from_bioproject(bioproject_id):
    """Fetch all BioSample IDs linked to a BioProject"""
    try:
        search_handle = Entrez.esearch(db="bioproject", term=bioproject_id)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if not search_results["IdList"]:
            return []

        bioproject_uid = search_results["IdList"][0]

        link_handle = Entrez.elink(dbfrom="bioproject", db="biosample", id=bioproject_uid)
        link_results = Entrez.read(link_handle)
        link_handle.close()

        biosample_ids = []
        if link_results and link_results[0]["LinkSetDb"]:
            for link in link_results[0]["LinkSetDb"][0]["Link"]:
                biosample_ids.append(link["Id"])

        return biosample_ids

    except Exception as e:
        print(f"  Error fetching BioSamples: {str(e)}")
        return []


def download_biosample_data(bioproject_id, biosample_ids, output_dir="."):
    """Download BioSample metadata"""
    try:
        if not biosample_ids:
            return None

        fetch_handle = Entrez.efetch(db="biosample", id=",".join(biosample_ids),
                                     rettype="full", retmode="text")
        data = fetch_handle.read()
        fetch_handle.close()

        output_file = os.path.join(output_dir, f"{bioproject_id}_biosample.txt")
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(data)

        return output_file
    except Exception as e:
        print(f"  Error downloading BioSample: {str(e)}")
        return None


def download_sra_runinfo(bioproject_id, output_dir="."):
    """Download SRA RunInfo"""
    try:
        search_handle = Entrez.esearch(db="sra", term=f"{bioproject_id}[BioProject]", retmax=10000)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if not search_results["IdList"]:
            print(f"  No SRA data found for {bioproject_id}")
            return None

        fetch_handle = Entrez.efetch(db="sra", id=",".join(search_results["IdList"]),
                                     rettype="runinfo", retmode="text")
        data = fetch_handle.read()
        fetch_handle.close()

        # Handle both bytes and str
        if isinstance(data, bytes):
            data = data.decode('utf-8')

        output_file = os.path.join(output_dir, f"{bioproject_id}_sra_runinfo.csv")
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(data)

        return output_file
    except Exception as e:
        print(f"  Error downloading SRA: {str(e)}")
        return None


class BioSampleParser:
    """Parser for BioSample text files"""

    def __init__(self):
        self.samples = []

    def parse_file(self, file_path):
        """Parse BioSample text file"""
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        sample_blocks = re.split(r'\n(?=\d+:\s+)', content)

        for block in sample_blocks:
            if not block.strip():
                continue
            sample_data = self._parse_sample_block(block)
            if sample_data:
                self.samples.append(sample_data)

        return pd.DataFrame(self.samples)

    def _parse_sample_block(self, block):
        """Parse single BioSample block"""
        data = {}

        # Extract BioSample accession
        biosample_match = re.search(r'Accession:\s+(\S+)', block)
        if biosample_match:
            data['BioSample'] = biosample_match.group(1)

        # Extract other fields
        patterns = {
            'Sample_Name': r'Identifiers:.*?Label:\s+(\S+)',
            'Organism': r'Organism:\s+(.+?)(?:\n|$)',
            'Accession': r'Accession:\s+(\S+)',
        }

        for key, pattern in patterns.items():
            match = re.search(pattern, block)
            if match:
                data[key] = match.group(1).strip()

        # Extract attributes section
        attr_section = re.search(r'Attributes:(.+?)(?=\n\n|\Z)', block, re.DOTALL)
        if attr_section:
            attr_text = attr_section.group(1)
            # Parse key-value pairs
            attr_matches = re.findall(r'/(\w+)="([^"]*)"', attr_text)
            for key, value in attr_matches:
                data[key] = value

        return data if data else None


def batch_download_ncbi(bioproject_list, output_dir, email):
    """
    Batch download NCBI metadata (BioSample + SRA)

    Returns:
        dict: {'biosample_files': [...], 'sra_files': [...]}
    """
    Entrez.email = email

    print(f"\n{'='*70}")
    print(f"NCBI Download: {len(bioproject_list)} BioProjects")
    print('='*70)

    biosample_files = []
    sra_files = []

    for idx, bioproject_id in enumerate(bioproject_list, 1):
        print(f"\n[{idx}/{len(bioproject_list)}] Processing {bioproject_id}...")

        # Download BioSample
        biosample_ids = get_biosamples_from_bioproject(bioproject_id)
        if biosample_ids:
            biosample_file = download_biosample_data(bioproject_id, biosample_ids, output_dir)
            if biosample_file:
                biosample_files.append(biosample_file)
                print(f"  ✓ BioSample: {len(biosample_ids)} samples")

        time.sleep(0.4)

        # Download SRA
        sra_file = download_sra_runinfo(bioproject_id, output_dir)
        if sra_file:
            sra_files.append(sra_file)
            print(f"  ✓ SRA RunInfo downloaded")

        time.sleep(0.4)

    print(f"\n{'='*70}")
    print(f"NCBI Download Complete:")
    print(f"  BioSample files: {len(biosample_files)}")
    print(f"  SRA files: {len(sra_files)}")
    print('='*70)

    return {
        'biosample_files': biosample_files,
        'sra_files': sra_files
    }


def parse_ncbi_data(biosample_files, sra_files, output_dir):
    """
    Parse NCBI BioSample and SRA files

    Returns:
        dict: {'biosample_df': DataFrame, 'sra_df': DataFrame}
    """
    print(f"\n{'='*70}")
    print("Parsing NCBI Data")
    print('='*70)

    # Parse BioSample files
    parser = BioSampleParser()
    biosample_dfs = []

    for file in biosample_files:
        try:
            df = parser.parse_file(file)
            biosample_dfs.append(df)
            print(f"  ✓ Parsed {file}: {len(df)} samples")
        except Exception as e:
            print(f"  ✗ Error parsing {file}: {e}")

    biosample_df = pd.concat(biosample_dfs, ignore_index=True) if biosample_dfs else pd.DataFrame()

    # Parse SRA files
    sra_dfs = []
    for file in sra_files:
        try:
            df = pd.read_csv(file)
            sra_dfs.append(df)
            print(f"  ✓ Parsed {file}: {len(df)} runs")
        except Exception as e:
            print(f"  ✗ Error parsing {file}: {e}")

    sra_df = pd.concat(sra_dfs, ignore_index=True) if sra_dfs else pd.DataFrame()

    print(f"\n{'='*70}")
    print(f"NCBI Parsing Complete:")
    print(f"  Total BioSamples: {len(biosample_df)}")
    print(f"  Total SRA runs: {len(sra_df)}")
    print('='*70)

    return {
        'biosample_df': biosample_df,
        'sra_df': sra_df
    }


def merge_ncbi_data(biosample_df, sra_df, output_dir):
    """
    Merge NCBI BioSample and SRA data

    Returns:
        DataFrame: Merged NCBI data
    """
    print(f"\n{'='*70}")
    print("Merging NCBI Data (SRA + BioSample)")
    print('='*70)

    if biosample_df.empty or sra_df.empty:
        print("  Warning: Empty dataframe, skipping merge")
        return pd.DataFrame()

    # Strategy 1: Merge by BioSample ID only
    df_merged_simple = sra_df.merge(
        biosample_df,
        left_on='BioSample',
        right_on='BioSample',
        how='outer',
        suffixes=('_sra', '_biosample'),
        indicator=True
    )

    matched_simple = len(df_merged_simple[df_merged_simple['_merge'] == 'both'])

    # Strategy 2: Merge by BioProject + BioSample
    if 'BioProject' in sra_df.columns:
        df_biosample_temp = biosample_df.copy()
        df_sra_temp = sra_df.copy()

        # Create composite keys (handle missing values)
        # Use .get() for Series/DataFrame columns, with default empty string
        bioproject_col = df_biosample_temp['BioProject_ID'] if 'BioProject_ID' in df_biosample_temp.columns else pd.Series([''] * len(df_biosample_temp))
        biosample_col = df_biosample_temp['BioSample'] if 'BioSample' in df_biosample_temp.columns else pd.Series([''] * len(df_biosample_temp))
        df_biosample_temp['merge_key'] = bioproject_col.fillna('').astype(str) + '_' + biosample_col.fillna('').astype(str)

        bioproject_sra_col = df_sra_temp['BioProject'] if 'BioProject' in df_sra_temp.columns else pd.Series([''] * len(df_sra_temp))
        biosample_sra_col = df_sra_temp['BioSample'] if 'BioSample' in df_sra_temp.columns else pd.Series([''] * len(df_sra_temp))
        df_sra_temp['merge_key'] = bioproject_sra_col.fillna('').astype(str) + '_' + biosample_sra_col.fillna('').astype(str)

        df_merged_combined = df_sra_temp.merge(
            df_biosample_temp,
            on='merge_key',
            how='outer',
            suffixes=('_sra', '_biosample'),
            indicator=True
        )

        matched_combined = len(df_merged_combined[df_merged_combined['_merge'] == 'both'])
    else:
        matched_combined = 0
        df_merged_combined = pd.DataFrame()

    # Select best strategy
    if matched_simple >= matched_combined:
        print(f"  Using Strategy 1: BioSample ID matching")
        print(f"  Matched: {matched_simple} records")
        df_final = df_merged_simple
        method = "biosample_id"
    else:
        print(f"  Using Strategy 2: BioProject + BioSample matching")
        print(f"  Matched: {matched_combined} records")
        df_final = df_merged_combined
        method = "bioproject_biosample"

    # Save matched data
    df_matched = df_final[df_final['_merge'] == 'both'].copy()
    df_matched = df_matched.drop('_merge', axis=1)

    matched_file = Path(output_dir) / f"ncbi_merged_{method}.csv"
    df_matched.to_csv(matched_file, index=False, encoding='utf-8-sig')
    print(f"  ✓ Saved: {matched_file} ({len(df_matched)} records)")

    print('='*70)

    return df_matched


# ============================================================================
# CNCB Download Functions (using iSeq)
# ============================================================================

def download_gsa_with_iseq(bioproject_id, output_dir):
    """
    Download GSA/CNCB metadata using iSeq (conda-installed)

    Args:
        bioproject_id: GSA/CNCB/PRJCA accession ID (e.g., PRJCA004523, CRA001234)
        output_dir: Directory to save downloads

    Returns:
        str: Path to metadata file if successful, None otherwise
    """
    try:
        print(f"  Downloading via iSeq: {bioproject_id}")

        cmd = ['iseq', '-i', bioproject_id, '-m', '-o', output_dir]

        result = subprocess.run(cmd,
                              capture_output=True,
                              text=True,
                              timeout=300)

        if result.returncode != 0:
            print(f"  ✗ iSeq failed for {bioproject_id}")
            if result.stderr:
                print(f"    stderr: {result.stderr[:200]}")
            return None

        metadata_files = glob.glob(os.path.join(output_dir, f"{bioproject_id}*"))
        metadata_files = [f for f in metadata_files if os.path.isfile(f) and f.endswith(('.xlsx', '.xls', '.csv', '.txt'))]

        if metadata_files:
            metadata_file = metadata_files[0]
            print(f"  ✓ Downloaded: {os.path.basename(metadata_file)}")
            return metadata_file
        else:
            print(f"  Warning: No metadata file found after iSeq download")
            return None

    except subprocess.TimeoutExpired:
        print(f"  ✗ iSeq download timeout for {bioproject_id}")
        return None
    except Exception as e:
        print(f"  ✗ Error downloading {bioproject_id} with iSeq: {e}")
        return None


def batch_download_cncb(bioproject_list, output_dir):
    """
    Batch download CNCB/GSA metadata using iSeq (conda-installed)

    Returns:
        list: List of downloaded metadata file paths
    """
    print(f"\n{'='*70}")
    print(f"CNCB/GSA Download (iSeq): {len(bioproject_list)} BioProjects")
    print('='*70)

    metadata_files = []

    for idx, bioproject_id in enumerate(bioproject_list, 1):
        print(f"\n[{idx}/{len(bioproject_list)}] Processing {bioproject_id}...")

        metadata_file = download_gsa_with_iseq(bioproject_id, output_dir)
        if metadata_file:
            metadata_files.append(metadata_file)

        if idx < len(bioproject_list):
            time.sleep(2)

    print(f"\n{'='*70}")
    print(f"CNCB Download Complete: {len(metadata_files)} files")
    print('='*70)

    return metadata_files


def parse_cncb_data(metadata_files, output_dir):
    """
    Parse CNCB metadata files (can be Excel, CSV, or other formats)
    Combine into single CSV

    Returns:
        DataFrame: Combined CNCB data
    """
    print(f"\n{'='*70}")
    print("Parsing CNCB Data (iSeq downloads)")
    print('='*70)

    cncb_dfs = []

    for metadata_file in metadata_files:
        try:
            filename = os.path.basename(metadata_file)
            file_ext = os.path.splitext(metadata_file)[1].lower()

            # Parse based on file extension
            if file_ext == '.xlsx' or file_ext == '.xls':
                df = pd.read_excel(metadata_file, sheet_name=0)
            elif file_ext == '.csv':
                df = pd.read_csv(metadata_file, encoding='utf-8-sig')
            elif file_ext == '.txt':
                try:
                    df = pd.read_csv(metadata_file, sep='\t', encoding='utf-8-sig')
                except Exception:
                    df = pd.read_csv(metadata_file, encoding='utf-8-sig')
            else:
                # Try to read as CSV by default
                df = pd.read_csv(metadata_file, encoding='utf-8-sig')

            # Add source file info
            df['Source_File'] = filename
            df['Source_Database'] = 'CNCB/GSA'

            cncb_dfs.append(df)
            print(f"  ✓ Parsed {filename}: {len(df)} samples, {len(df.columns)} columns")
        except Exception as e:
            print(f"  ✗ Error parsing {metadata_file}: {e}")

    if not cncb_dfs:
        print("  Warning: No CNCB data parsed")
        return pd.DataFrame()

    cncb_df = pd.concat(cncb_dfs, ignore_index=True)

    # Save combined CNCB data
    cncb_file = Path(output_dir) / "cncb_combined.csv"
    cncb_df.to_csv(cncb_file, index=False, encoding='utf-8-sig')
    print(f"\n  ✓ Saved: {cncb_file} ({len(cncb_df)} records)")

    print('='*70)

    return cncb_df


# ============================================================================
# Final Merge: NCBI + CNCB
# ============================================================================

def normalize_column_names(df, database):
    """
    Normalize column names for merging
    Standardize common fields across NCBI and CNCB
    """
    # CNCB column mapping
    if database == 'CNCB':
        rename_map = {
            'Accession': 'BioSample',
            'Sample name': 'Sample_Name',
            'Project accession': 'BioProject',
            'Collection date': 'collection_date',
            'Geographic location': 'geo_loc_name',
            'Isolation source': 'isolation_source',
        }

        # Apply case-insensitive matching
        for old_col in df.columns:
            for key, value in rename_map.items():
                if old_col.strip() == key:
                    df = df.rename(columns={old_col: value})
                    break

    return df


def merge_ncbi_cncb_data(ncbi_df, cncb_df, output_dir):
    """
    Final merge: Combine NCBI and CNCB data by matching column names

    Returns:
        DataFrame: Final merged data
    """
    print(f"\n{'='*70}")
    print("Final Merge: NCBI + CNCB")
    print('='*70)

    if ncbi_df.empty and cncb_df.empty:
        print("  Warning: Both dataframes are empty")
        return pd.DataFrame()

    if ncbi_df.empty:
        print("  Only CNCB data available")
        final_df = cncb_df
    elif cncb_df.empty:
        print("  Only NCBI data available")
        final_df = ncbi_df
    else:
        # Normalize column names
        print("\n  Normalizing column names...")
        ncbi_df = normalize_column_names(ncbi_df, 'NCBI')
        cncb_df = normalize_column_names(cncb_df, 'CNCB')

        # Add database source tags
        ncbi_df['Source_Database'] = 'NCBI'
        cncb_df['Source_Database'] = 'CNCB/GSA'

        # Find common columns
        ncbi_cols = set(ncbi_df.columns)
        cncb_cols = set(cncb_df.columns)
        common_cols = ncbi_cols.intersection(cncb_cols)

        print(f"  NCBI columns: {len(ncbi_cols)}")
        print(f"  CNCB columns: {len(cncb_cols)}")
        print(f"  Common columns: {len(common_cols)}")

        if common_cols:
            print(f"  Common fields: {', '.join(sorted(common_cols)[:10])}...")

        # Concatenate (vertical merge - combine rows)
        print("\n  Combining datasets (vertical concatenation)...")
        final_df = pd.concat([ncbi_df, cncb_df], ignore_index=True, sort=False)

    # Save final merged data
    final_file = Path(output_dir) / "all_metadata_merged.csv"
    final_df.to_csv(final_file, index=False, encoding='utf-8-sig')

    print(f"\n{'='*70}")
    print(f"Final Merge Complete:")
    print(f"  Total records: {len(final_df)}")
    print(f"  Total columns: {len(final_df.columns)}")
    print(f"  Output: {final_file}")
    print('='*70)

    return final_df


# ============================================================================
# Main Pipeline
# ============================================================================

def run_unified_pipeline(input_folder, output_folder, email):
    """
    Execute unified metadata download pipeline

    Workflow:
        1. Read and classify BioProject IDs (NCBI vs CNCB)
        2. Download NCBI data (BioSample + SRA)
        3. Download CNCB data (Excel files)
        4. Parse and merge NCBI data (Stage 1)
        5. Parse CNCB data and convert to CSV (Stage 2)
        6. Final merge: NCBI + CNCB (Stage 3)

    Args:
        input_folder: Path to folder containing BioProject ID txt files
        output_folder: Path to output folder
        email: Email for NCBI Entrez API

    Returns:
        dict: Results summary
    """
    print("\n" + "="*70)
    print("UNIFIED METADATA DOWNLOAD PIPELINE")
    print("="*70)
    print(f"Input:  {input_folder}")
    print(f"Output: {output_folder}")
    print(f"Email:  {email}")
    print("="*70)

    # Create output directory
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    # Step 1: Read and classify BioProject IDs
    print("\n[Step 1] Reading and classifying BioProject IDs...")
    bioproject_ids = read_bioproject_ids(input_folder)

    if not bioproject_ids:
        print("ERROR: No BioProject IDs found")
        return None

    classified = classify_bioprojects(bioproject_ids)

    print(f"\n  Classification:")
    print(f"    NCBI (PRJNA/PRJEB/PRJDB): {len(classified['ncbi'])} projects")
    print(f"    CNCB/GSA (PRJCA):         {len(classified['cncb'])} projects")

    # Initialize results
    ncbi_df = pd.DataFrame()
    cncb_df = pd.DataFrame()

    # Step 2: Process NCBI BioProjects
    if classified['ncbi']:
        print("\n[Step 2] Processing NCBI BioProjects...")

        # Download
        ncbi_downloads = batch_download_ncbi(
            classified['ncbi'],
            output_folder,
            email
        )

        # Parse
        ncbi_parsed = parse_ncbi_data(
            ncbi_downloads['biosample_files'],
            ncbi_downloads['sra_files'],
            output_folder
        )

        # Merge (Stage 1)
        ncbi_df = merge_ncbi_data(
            ncbi_parsed['biosample_df'],
            ncbi_parsed['sra_df'],
            output_folder
        )
    else:
        print("\n[Step 2] No NCBI BioProjects to process")

    # Step 3: Process CNCB BioProjects
    if classified['cncb']:
        print("\n[Step 3] Processing CNCB/GSA BioProjects...")

        # Download
        xlsx_files = batch_download_cncb(
            classified['cncb'],
            output_folder
        )

        # Parse and convert to CSV (Stage 2)
        cncb_df = parse_cncb_data(xlsx_files, output_folder)
    else:
        print("\n[Step 3] No CNCB/GSA BioProjects to process")

    # Step 4: Final merge (Stage 3)
    print("\n[Step 4] Final merge: Combining NCBI and CNCB data...")
    final_df = merge_ncbi_cncb_data(ncbi_df, cncb_df, output_folder)

    # Summary
    print("\n" + "="*70)
    print("PIPELINE COMPLETE")
    print("="*70)
    print(f"Total BioProjects processed: {len(bioproject_ids)}")
    print(f"  - NCBI: {len(classified['ncbi'])}")
    print(f"  - CNCB/GSA: {len(classified['cncb'])}")
    print(f"\nFinal output:")
    print(f"  - Total records: {len(final_df)}")
    print(f"  - Total columns: {len(final_df.columns)}")
    print(f"  - File: {output_folder}/all_metadata_merged.csv")
    print("="*70 + "\n")

    return {
        'bioproject_ids': bioproject_ids,
        'classified': classified,
        'ncbi_df': ncbi_df,
        'cncb_df': cncb_df,
        'final_df': final_df
    }


# ============================================================================
# Command Line Interface
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Unified Metadata Download Pipeline (NCBI + CNCB/GSA)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download metadata for mixed NCBI and CNCB BioProjects
  python3 unified_metadata_downloader.py \\
      --input bioproject_ids/ \\
      --output metadata_output/ \\
      --email your@email.com

Input Format:
  Place .txt files containing BioProject IDs (one per line) in input folder:
    PRJNA123456
    PRJEB123456
    PRJCA123456

Output Files:
  - ncbi_merged_*.csv         : Merged NCBI data (SRA + BioSample)
  - cncb_combined.csv         : Combined CNCB data (converted from Excel)
  - all_metadata_merged.csv   : Final merged dataset (NCBI + CNCB)
        """
    )

    parser.add_argument('-i', '--input', required=True,
                       help='Input directory containing BioProject ID txt files')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory for metadata files')
    parser.add_argument('-e', '--email', required=True,
                       help='Email address (required by NCBI)')

    args = parser.parse_args()

    # Validate input directory
    if not os.path.isdir(args.input):
        print(f"ERROR: Input directory not found: {args.input}")
        sys.exit(1)

    # Run pipeline
    try:
        results = run_unified_pipeline(
            args.input,
            args.output,
            args.email
        )

        if results is None:
            sys.exit(1)

        sys.exit(0)

    except KeyboardInterrupt:
        print("\n\nPipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
