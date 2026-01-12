#!/usr/bin/env python3
"""
NCBI Metadata Download and Processing Pipeline
Complete workflow for downloading and merging BioProject, BioSample, and SRA data

Usage:
    python3 ncbi_metadata_downloader.py --input /path/to/bioproject_ids/ --output /path/to/output/ --email your@email.com
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


# ============================================================================
# Core Functions
# ============================================================================

def read_bioproject_ids(folder_path):
    """Read BioProject IDs from all txt files in folder"""
    txt_files = glob.glob(os.path.join(folder_path, "*.txt"))

    if not txt_files:
        print(f"ERROR: No txt files found in {folder_path}")
        return pd.DataFrame()

    df_list = []
    for file in txt_files:
        try:
            temp_df = pd.read_csv(file, sep='\t', engine='python', header=None,
                                  names=['BioProject_ID'])
            temp_df['source_file'] = os.path.basename(file)
            df_list.append(temp_df)
        except Exception as e:
            print(f"  Warning: Failed to read {file}: {e}")
            continue

    if not df_list:
        print("ERROR: No valid BioProject IDs found")
        return pd.DataFrame()

    df = pd.concat(df_list, ignore_index=True)
    df = df.drop_duplicates(subset=['BioProject_ID'], keep='first')

    print(f"Found {df.shape[0]} unique BioProject IDs from {len(txt_files)} files")
    return df


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


def batch_download_biosamples(bioproject_list, output_dir="."):
    """Batch download BioSample data for multiple BioProjects"""
    print(f"\n[Step 2/6] Downloading BioSample data for {len(bioproject_list)} BioProjects")

    results = []
    for i, bioproject_id in enumerate(bioproject_list, 1):
        print(f"  [{i}/{len(bioproject_list)}] {bioproject_id}...", end=" ", flush=True)

        biosample_ids = get_biosamples_from_bioproject(bioproject_id)

        result = {
            'BioProject_ID': bioproject_id,
            'BioSample_Count': len(biosample_ids),
            'TXT_File': None
        }

        if biosample_ids:
            txt_file = download_biosample_data(bioproject_id, biosample_ids, output_dir)
            result['TXT_File'] = txt_file
            print(f"{len(biosample_ids)} BioSamples")
        else:
            print("No BioSamples found")

        results.append(result)
        time.sleep(0.4)

    return pd.DataFrame(results)


class BioSampleParser:
    """Parser for BioSample text files"""

    def parse_file(self, file_path):
        """Parse BioSample file and extract sample information"""
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        sample_blocks = re.split(r'\n(?=\d+:\s+)', content)
        samples = []

        for block in sample_blocks:
            if not block.strip():
                continue
            sample_data = self._parse_sample_block(block)
            if sample_data:
                samples.append(sample_data)

        return pd.DataFrame(samples)

    def _parse_sample_block(self, block):
        """Parse individual sample block"""
        data = {}
        lines = block.split('\n')

        first_line = lines[0]
        match = re.match(r'(\d+):\s*(.+)', first_line)
        if match:
            data['Sample_Number'] = int(match.group(1))
            data['Sample_Name'] = match.group(2).strip()

        current_section = None
        attributes_lines = []
        description_text = []

        for line in lines[1:]:
            line_stripped = line.strip()

            if line_stripped.startswith('Identifiers:'):
                current_section = 'identifiers'
                self._parse_identifiers(line_stripped.replace('Identifiers:', '').strip(), data)
            elif line_stripped.startswith('Organism:'):
                data['Organism'] = line_stripped.replace('Organism:', '').strip()
            elif line_stripped.startswith('Attributes:'):
                current_section = 'attributes'
            elif line_stripped.startswith('Description:'):
                current_section = 'description'
            elif line_stripped.startswith('Keywords:'):
                data['Keywords'] = line_stripped.replace('Keywords:', '').strip()
            elif line_stripped.startswith('Accession:'):
                parts = line_stripped.split('\t')
                for part in parts:
                    if 'Accession:' in part:
                        data['Accession'] = part.replace('Accession:', '').strip()
                    elif 'ID:' in part:
                        data['ID'] = part.replace('ID:', '').strip()
            elif current_section == 'identifiers' and line_stripped:
                self._parse_identifiers(line_stripped, data)
            elif current_section == 'attributes' and line_stripped.startswith('/'):
                attributes_lines.append(line_stripped)
            elif current_section == 'description' and line_stripped:
                description_text.append(line_stripped)

        if description_text:
            data['Description'] = ' '.join(description_text)

        for attr_line in attributes_lines:
            self._parse_attribute(attr_line, data)

        return data

    def _parse_identifiers(self, text, data):
        """Parse identifier line"""
        biosample_match = re.search(r'BioSample:\s*(\S+)', text)
        if biosample_match:
            data['BioSample'] = biosample_match.group(1).rstrip(';')

        sample_match = re.search(r'Sample name:\s*([^;]+)', text)
        if sample_match:
            data['Sample_Name_Full'] = sample_match.group(1).strip()

        sra_match = re.search(r'SRA:\s*(\S+)', text)
        if sra_match:
            data['SRA'] = sra_match.group(1)

    def _parse_attribute(self, line, data):
        """Parse single attribute line"""
        match = re.match(r'/([^=]+)="([^"]*)"', line)
        if match:
            key = match.group(1).strip()
            value = match.group(2).strip()
            clean_key = key.replace(' ', '_').replace('-', '_')
            data[clean_key] = value


def batch_parse_biosample_files(directory, output_file="all_biosamples_combined.csv"):
    """Parse all BioSample txt files in directory"""
    print(f"\n[Step 3/6] Parsing BioSample files")

    directory = Path(directory)
    txt_files = list(directory.glob("*_biosample.txt"))

    if not txt_files:
        print("  No BioSample files found")
        return None

    all_dataframes = []
    parser = BioSampleParser()

    for i, file_path in enumerate(txt_files, 1):
        print(f"  [{i}/{len(txt_files)}] {file_path.name}...", end=" ", flush=True)
        try:
            df = parser.parse_file(file_path)
            df['Source_File'] = file_path.name

            bioproject_match = re.search(r'(PRJ[A-Z]+\d+)', file_path.name)
            if bioproject_match:
                df['BioProject_ID'] = bioproject_match.group(1)

            all_dataframes.append(df)
            print(f"{len(df)} samples")
        except Exception as e:
            print(f"Error: {str(e)}")
            continue

    if not all_dataframes:
        return None

    df_combined = pd.concat(all_dataframes, ignore_index=True)

    output_path = directory / output_file
    df_combined.to_csv(output_path, index=False, encoding='utf-8-sig')

    print(f"  Combined: {len(df_combined)} samples, {len(df_combined.columns)} columns")
    print(f"  Saved to: {output_path}")

    return df_combined


def get_sra_runs_from_bioproject(bioproject_id):
    """Fetch all SRA Run IDs linked to a BioProject"""
    try:
        search_query = f"{bioproject_id}[BioProject]"
        search_handle = Entrez.esearch(db="sra", term=search_query, retmax=10000)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        sra_ids = search_results["IdList"]
        return sra_ids
    except Exception as e:
        print(f"  Error fetching SRA runs: {str(e)}")
        return []


def download_sra_metadata(bioproject_id, sra_ids, output_dir="."):
    """Download SRA RunInfo metadata"""
    if not sra_ids:
        return None

    try:
        fetch_handle = Entrez.efetch(db="sra", id=",".join(sra_ids),
                                     rettype="runinfo", retmode="text")
        txt_data = fetch_handle.read()
        fetch_handle.close()

        txt_file = os.path.join(output_dir, f"{bioproject_id}_sra_runinfo.csv")

        if isinstance(txt_data, bytes):
            with open(txt_file, "w", encoding="utf-8") as f:
                f.write(txt_data.decode('utf-8'))
        else:
            with open(txt_file, "w", encoding="utf-8") as f:
                f.write(txt_data)

        return txt_file
    except Exception as e:
        print(f"  Error downloading SRA: {str(e)}")
        return None


def batch_download_sra_data(bioproject_list, output_dir="."):
    """Batch download SRA Run data for multiple BioProjects"""
    print(f"\n[Step 4/6] Downloading SRA data for {len(bioproject_list)} BioProjects")

    results = []
    for i, bioproject_id in enumerate(bioproject_list, 1):
        print(f"  [{i}/{len(bioproject_list)}] {bioproject_id}...", end=" ", flush=True)

        sra_ids = get_sra_runs_from_bioproject(bioproject_id)

        result = {
            'BioProject_ID': bioproject_id,
            'SRA_Run_Count': len(sra_ids),
            'RunInfo_File': None
        }

        if sra_ids:
            txt_file = download_sra_metadata(bioproject_id, sra_ids, output_dir)
            result['RunInfo_File'] = txt_file
            print(f"{len(sra_ids)} runs")
        else:
            print("No SRA runs found")

        results.append(result)
        time.sleep(0.4)

    return pd.DataFrame(results)


def parse_sra_runinfo_files(directory, output_file="all_sra_runs_combined.csv"):
    """Parse all SRA RunInfo CSV files"""
    print(f"\n[Step 5/6] Parsing SRA RunInfo files")

    directory = Path(directory)
    csv_files = list(directory.glob("*_sra_runinfo.csv"))

    if not csv_files:
        print("  No SRA RunInfo files found")
        return None

    all_dfs = []

    for i, csv_file in enumerate(csv_files, 1):
        print(f"  [{i}/{len(csv_files)}] {csv_file.name}...", end=" ", flush=True)
        try:
            df = pd.read_csv(csv_file)
            df['Source_File'] = csv_file.name

            bioproject_match = re.search(r'(PRJ[A-Z]+\d+)', csv_file.name)
            if bioproject_match:
                df['BioProject_ID_from_file'] = bioproject_match.group(1)

            all_dfs.append(df)
            print(f"{len(df)} runs")
        except Exception as e:
            print(f"Error: {str(e)}")
            continue

    if not all_dfs:
        return None

    df_combined = pd.concat(all_dfs, ignore_index=True)

    output_path = directory / output_file
    df_combined.to_csv(output_path, index=False, encoding='utf-8-sig')

    print(f"  Combined: {len(df_combined)} runs, {len(df_combined.columns)} columns")
    print(f"  Saved to: {output_path}")

    return df_combined


def merge_sra_biosample_data(biosample_file, sra_file, output_dir):
    """Merge SRA and BioSample data using optimal strategy"""
    print(f"\n[Step 6/6] Merging SRA and BioSample data")

    # Load data
    df_biosample = pd.read_csv(biosample_file)
    df_sra = pd.read_csv(sra_file)

    print(f"  BioSample records: {len(df_biosample)}")
    print(f"  SRA records: {len(df_sra)}")

    # Check overlap
    biosample_ids_bio = set(df_biosample['BioSample'].dropna())
    biosample_ids_sra = set(df_sra['BioSample'].dropna())
    common_ids = biosample_ids_bio & biosample_ids_sra

    print(f"  Common BioSample IDs: {len(common_ids)}")

    # Strategy 1: Merge by BioSample ID only
    df_merged_simple = df_sra.merge(
        df_biosample,
        left_on='BioSample',
        right_on='BioSample',
        how='outer',
        suffixes=('_sra', '_biosample'),
        indicator=True
    )

    merge_stats = df_merged_simple['_merge'].value_counts()
    matched_simple = merge_stats.get('both', 0)

    print(f"  Strategy 1 (BioSample ID): {matched_simple} matches")

    # Strategy 2: Merge by BioProject + BioSample
    df_biosample_temp = df_biosample.copy()
    df_sra_temp = df_sra.copy()

    df_biosample_temp['merge_key'] = (
        df_biosample_temp['BioProject_ID'].astype(str) + '_' +
        df_biosample_temp['BioSample'].astype(str)
    )
    df_sra_temp['merge_key'] = (
        df_sra_temp['BioProject'].astype(str) + '_' +
        df_sra_temp['BioSample'].astype(str)
    )

    df_merged_combined = df_sra_temp.merge(
        df_biosample_temp,
        on='merge_key',
        how='outer',
        suffixes=('_sra', '_biosample'),
        indicator=True
    )

    merge_stats_combined = df_merged_combined['_merge'].value_counts()
    matched_combined = merge_stats_combined.get('both', 0)

    print(f"  Strategy 2 (BioProject+BioSample): {matched_combined} matches")

    # Select best strategy
    if matched_simple >= matched_combined:
        print(f"  Using Strategy 1")
        df_final = df_merged_simple
        method = "biosample_id"
    else:
        print(f"  Using Strategy 2")
        df_final = df_merged_combined
        method = "bioproject_biosample"

    # Save results
    output_dir = Path(output_dir)

    # Matched data
    df_matched = df_final[df_final['_merge'] == 'both'].copy()
    df_matched = df_matched.drop('_merge', axis=1)
    matched_file = output_dir / f"matched_sra_biosample_{method}.csv"
    df_matched.to_csv(matched_file, index=False, encoding='utf-8-sig')
    print(f"  Matched data saved: {matched_file} ({len(df_matched)} records)")

    # Unmatched data
    df_sra_only = df_final[df_final['_merge'] == 'left_only'].copy()
    df_biosample_only = df_final[df_final['_merge'] == 'right_only'].copy()

    if len(df_sra_only) > 0:
        sra_only_file = output_dir / f"unmatched_sra_only_{method}.csv"
        df_sra_only.to_csv(sra_only_file, index=False, encoding='utf-8-sig')
        print(f"  SRA-only data saved: {sra_only_file} ({len(df_sra_only)} records)")

    if len(df_biosample_only) > 0:
        biosample_only_file = output_dir / f"unmatched_biosample_only_{method}.csv"
        df_biosample_only.to_csv(biosample_only_file, index=False, encoding='utf-8-sig')
        print(f"  BioSample-only data saved: {biosample_only_file} ({len(df_biosample_only)} records)")

    return {
        'merged': df_final,
        'matched': df_matched,
        'sra_only': df_sra_only,
        'biosample_only': df_biosample_only,
        'method': method
    }


# ============================================================================
# Main Pipeline
# ============================================================================

def run_complete_pipeline(input_folder, output_folder, email):
    """
    Execute complete NCBI data download and processing pipeline

    Args:
        input_folder: Path to folder containing BioProject ID txt files
        output_folder: Path to output folder
        email: Email for NCBI Entrez

    Returns:
        Dictionary containing all results
    """

    # Configure Entrez
    Entrez.email = email

    os.makedirs(output_folder, exist_ok=True)

    print("="*70)
    print("NCBI Data Download and Processing Pipeline")
    print("="*70)

    # Step 1: Read BioProject IDs
    print("\n[Step 1/6] Reading BioProject IDs")
    df_bioprojects = read_bioproject_ids(input_folder)

    if df_bioprojects.empty:
        print("ERROR: No BioProject IDs found. Exiting.")
        return None

    bioproject_ids = df_bioprojects['BioProject_ID'].tolist()
    print(f"  BioProjects: {bioproject_ids}")

    # Step 2-3: Download and parse BioSample data
    biosample_download_results = batch_download_biosamples(bioproject_ids, output_folder)
    df_biosamples = batch_parse_biosample_files(output_folder, "all_biosamples_combined.csv")

    # Step 4-5: Download and parse SRA data
    sra_download_results = batch_download_sra_data(bioproject_ids, output_folder)
    df_sra = parse_sra_runinfo_files(output_folder, "all_sra_runs_combined.csv")

    # Step 6: Merge data
    biosample_file = os.path.join(output_folder, "all_biosamples_combined.csv")
    sra_file = os.path.join(output_folder, "all_sra_runs_combined.csv")

    if os.path.exists(biosample_file) and os.path.exists(sra_file):
        merge_results = merge_sra_biosample_data(biosample_file, sra_file, output_folder)

        # Final summary
        print("\n" + "="*70)
        print("Pipeline Complete - Summary")
        print("="*70)
        print(f"BioProjects: {len(bioproject_ids)}")
        print(f"BioSamples: {len(df_biosamples) if df_biosamples is not None else 0}")
        print(f"SRA Runs: {len(df_sra) if df_sra is not None else 0}")
        print(f"Matched records: {len(merge_results['matched'])}")
        print(f"SRA-only records: {len(merge_results['sra_only'])}")
        print(f"BioSample-only records: {len(merge_results['biosample_only'])}")
        print(f"\nAll files saved to: {output_folder}")
        print("="*70)

        return {
            'bioproject_ids': bioproject_ids,
            'biosamples': df_biosamples,
            'sra_runs': df_sra,
            'merge_results': merge_results
        }
    else:
        print("\nWarning: Missing BioSample or SRA files, skipping merge")
        return None


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='NCBI Metadata Download and Processing Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Download metadata for BioProjects in input folder
  python3 ncbi_metadata_downloader.py -i ./bioproject_ids/ -o ./output/ -e your@email.com

  # Use custom email
  python3 ncbi_metadata_downloader.py -i ./input -o ./output -e researcher@university.edu

Input folder should contain .txt files with BioProject IDs (one per line).
        '''
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input folder containing BioProject ID txt files')
    parser.add_argument('-o', '--output', required=True,
                        help='Output folder for downloaded metadata')
    parser.add_argument('-e', '--email', required=True,
                        help='Email address for NCBI Entrez (required by NCBI)')

    args = parser.parse_args()

    # Validate input folder
    if not os.path.isdir(args.input):
        print(f"ERROR: Input folder does not exist: {args.input}")
        sys.exit(1)

    # Run pipeline
    results = run_complete_pipeline(args.input, args.output, args.email)

    if results is None:
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
