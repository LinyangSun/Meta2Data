#!/usr/bin/env python3
"""
Unified Metadata Download Pipeline V4.3

Supports:
- BioProject IDs (PRJNA*, PRJEB*, PRJDB*, PRJCA*)
- BioSample IDs (SAMN*, SAMEA*, SAMC*)
- Keyword search mode (NCBI + CNCB)

Dependencies: biopython, pandas, openpyxl, requests
"""

from Bio import Entrez
import time
import os
import csv
import pandas as pd
import re
import requests
from io import StringIO
from pathlib import Path
import glob
import argparse
import sys
import json
import threading
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import quote
from itertools import product
import hashlib

# ============================================================================
# Configuration
# ============================================================================

DEFAULT_MAX_WORKERS = 3
DEFAULT_RETRY_ATTEMPTS = 5
DEFAULT_RETRY_DELAY = 1.5
DEFAULT_REQUEST_DELAY = 0.4
BATCH_SIZE = 500

UNIFIED_PATTERNS = {
    'Run': re.compile(r'^[CEDS]RR\d+$'),
    'BioProject': re.compile(r'^PRJ[CEDN][A-Z]\d+$'),
    'BioSample': re.compile(r'^SAM[CEDN][A-Z]?\d+$'),
    'Experiment': re.compile(r'^[CEDS]RX\d+$'),
}

ID_PATTERNS = {
    'biosample': re.compile(r'^SAM[CEDN][A-Z]?\d+$'),
    'cncb': re.compile(r'^PRJC[A-Z]\d+$'),
    'ncbi': re.compile(r'^PRJ[EDN][A-Z]\d+$'),
}

PRIORITY_COLUMNS = ['Run', 'BioProject', 'BioSample', 'Experiment']

STATUS_HAS_DATA = 'has_data'
STATUS_NO_DATA = 'no_data'
STATUS_NO_RUN = 'no_run_info'
STATUS_DOWNLOAD_ERROR = 'download_error'
STATUS_INVALID_FORMAT = 'invalid_format'

COLUMN_RENAME_DICT = None
COLUMN_RENAME_JSON_PATH = Path(__file__).parent.parent / "docs" / "NCBI_Biosample.json"


# ============================================================================
# Helpers
# ============================================================================

def generate_fake_email():
    timestamp = str(int(time.time() / 86400))
    hash_val = hashlib.md5(timestamp.encode()).hexdigest()[:8]
    return f"meta2data_{hash_val}@research.example.com"


def load_column_rename_dict():
    global COLUMN_RENAME_DICT
    if COLUMN_RENAME_DICT is None:
        try:
            with open(COLUMN_RENAME_JSON_PATH, 'r', encoding='utf-8') as f:
                COLUMN_RENAME_DICT = json.load(f)
            print(f"  Loaded {len(COLUMN_RENAME_DICT)} column rename rules")
        except Exception as e:
            print(f"  Warning: Failed to load column rename dict: {e}")
            COLUMN_RENAME_DICT = {}
    return COLUMN_RENAME_DICT


def setup_entrez(api_key=None):
    """Set up Entrez with auto-generated email. Call once at pipeline start."""
    Entrez.email = generate_fake_email()
    if api_key:
        Entrez.api_key = api_key


def retry_wrapper(func, max_retries=DEFAULT_RETRY_ATTEMPTS, retry_delay=DEFAULT_RETRY_DELAY):
    last_error = None
    for attempt in range(max_retries):
        try:
            return func()
        except Exception as e:
            last_error = e
            if attempt < max_retries - 1:
                time.sleep(retry_delay * (2 ** attempt))
    raise last_error


def classify_ids(id_list):
    """Classify input IDs into biosample, cncb, ncbi, and unknown groups."""
    result = {'biosample': [], 'cncb': [], 'ncbi': [], 'unknown': []}
    for item in id_list:
        matched = False
        for key, pattern in ID_PATTERNS.items():
            if pattern.match(item):
                result[key].append(item)
                matched = True
                break
        if not matched:
            result['unknown'].append(item)
    return result


def validate_run_data(df):
    """Check DataFrame has valid Run accessions."""
    if df.empty or 'Run' not in df.columns:
        return False
    return df['Run'].astype(str).apply(lambda x: bool(UNIFIED_PATTERNS['Run'].match(x))).any()


def read_input_ids(folder_path):
    """Read IDs from all txt files in folder."""
    txt_files = glob.glob(os.path.join(folder_path, "*.txt"))
    if not txt_files:
        print(f"ERROR: No txt files found in {folder_path}")
        return []

    all_ids = []
    for f in txt_files:
        try:
            with open(f, 'r') as fh:
                for line in fh:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        all_ids.append(line)
        except Exception as e:
            print(f"  Warning: Failed to read {f}: {e}")

    unique_ids = list(dict.fromkeys(all_ids))
    print(f"Found {len(unique_ids)} unique IDs from {len(txt_files)} files")
    return unique_ids


# ============================================================================
# Column Cleaning
# ============================================================================

def remove_empty_columns(df):
    """Remove columns where all values are empty/NaN."""
    df = df.dropna(axis=1, how='all')
    mask = df.astype(str).apply(lambda col: col.str.strip()).replace('', pd.NA).ne('nan').any()
    return df.loc[:, mask]


def remove_duplicate_columns(df):
    """Remove columns with identical content."""
    seen = {}
    keep = []
    for col in df.columns:
        key = tuple(df[col].astype(str).fillna('').tolist())
        if key not in seen:
            seen[key] = col
            keep.append(col)
    return df[keep]


def clean_and_standardize_columns(df, source_prefix=None):
    """Clean column names and auto-detect standard fields by content."""
    if df.empty:
        return df

    df = remove_empty_columns(df)
    df = remove_duplicate_columns(df)

    new_columns = []
    cols_to_drop = []
    renamed = {}

    for col in df.columns:
        new_name = col
        new_name = re.sub(r'^(Sample|Experiment|Run)_', '', new_name)
        new_name = re.sub(r'(_sra|_biosample|_x|_y)$', '', new_name)

        if new_name == 'ID' and len(new_columns) == 0:
            cols_to_drop.append(col)
            continue

        col_values = df[col].astype(str).str.strip().unique()
        col_values = [v for v in col_values if v and v.lower() not in ['nan', 'none', '']]

        for standard_name, pattern in UNIFIED_PATTERNS.items():
            if col_values and all(pattern.match(str(v)) for v in col_values):
                if standard_name not in new_columns and standard_name not in renamed.values():
                    renamed[col] = standard_name
                    new_name = standard_name
                    break

        new_columns.append(new_name)

    if cols_to_drop:
        df = df.drop(columns=cols_to_drop)

    df.columns = new_columns

    first_cols = [c for c in PRIORITY_COLUMNS if c in df.columns]
    other_cols = [c for c in df.columns if c not in PRIORITY_COLUMNS]
    return df[first_cols + other_cols]


def apply_column_rename_from_dict(df):
    """Apply column renaming from NCBI_Biosample.json dictionary."""
    rename_dict = load_column_rename_dict()
    if not rename_dict:
        return df

    variant_to_standard = {}
    for standard_name, variants in rename_dict.items():
        for variant in variants:
            variant_to_standard[variant.lower().strip()] = standard_name

    rename_map = {}
    for col in df.columns:
        standard = variant_to_standard.get(col.lower().strip())
        if standard and standard not in df.columns:
            rename_map[col] = standard

    if rename_map:
        df = df.rename(columns=rename_map)
        print(f"  Renamed {len(rename_map)} columns via dictionary")

    return df


def _to_camelcase(name):
    """Convert 'library_strategy' -> 'LibraryStrategy'."""
    parts = re.split(r'[\s_-]+', name.strip())
    return ''.join(part.capitalize() for part in parts if part)


def apply_camelcase_normalization(df):
    """Normalize column names to CamelCase. Merge duplicates intelligently."""
    rename_map = {}
    merge_needed = {}

    for col in df.columns:
        if ' ' in col or '_' in col or '-' in col:
            normalized = _to_camelcase(col)
            if normalized != col:
                if normalized in df.columns:
                    merge_needed[col] = normalized
                else:
                    rename_map[col] = normalized

    # Merge duplicate columns row-by-row
    for source_col, target_col in merge_needed.items():
        merged = []
        for idx in df.index:
            t = str(df.loc[idx, target_col]).strip() if pd.notna(df.loc[idx, target_col]) else ""
            s = str(df.loc[idx, source_col]).strip() if pd.notna(df.loc[idx, source_col]) else ""

            if not t and not s:
                merged.append(None)
            elif t and not s:
                merged.append(df.loc[idx, target_col])
            elif not t and s:
                merged.append(df.loc[idx, source_col])
            elif t.lower() == s.lower():
                merged.append(df.loc[idx, target_col])
            else:
                merged.append(f"{df.loc[idx, target_col]}_{df.loc[idx, source_col]}")

        df[target_col] = merged
        df = df.drop(columns=[source_col])

    if merge_needed:
        print(f"  Merged {len(merge_needed)} duplicate columns")

    if rename_map:
        df = df.rename(columns=rename_map)
        print(f"  Normalized {len(rename_map)} columns to CamelCase")

    return df


def standardize_columns(df):
    """Apply all three column standardization steps."""
    df = clean_and_standardize_columns(df)
    df = apply_column_rename_from_dict(df)
    df = apply_camelcase_normalization(df)
    return df


# ============================================================================
# State Management (Checkpoint System)
# ============================================================================

class StateManager:
    def __init__(self, output_dir):
        self.output_dir = Path(output_dir)
        self.checkpoint_dir = self.output_dir / "checkpoints"
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.state_file = self.checkpoint_dir / "download_state.json"
        self.lock = threading.Lock()
        self.state = self._load_state()

    def _load_state(self):
        if self.state_file.exists():
            try:
                with open(self.state_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except Exception:
                pass
        return {
            'completed_projects': {},
            'failed_downloads': [],
            'status_records': {},
            'last_update': None
        }

    def save_state(self):
        with self.lock:
            self.state['last_update'] = datetime.now().isoformat()
            temp_file = self.state_file.with_suffix('.tmp')
            try:
                with open(temp_file, 'w', encoding='utf-8') as f:
                    json.dump(self.state, f, indent=2)
                temp_file.replace(self.state_file)
            except Exception:
                pass

    def is_completed(self, project_id):
        csv_path = self.output_dir / f"{project_id}.processed.csv"
        with self.lock:
            return project_id in self.state['completed_projects'] and csv_path.exists()

    def mark_complete(self, project_id, file_path, row_count, source):
        with self.lock:
            self.state['completed_projects'][project_id] = {
                'file': str(file_path), 'rows': row_count,
                'source': source, 'timestamp': datetime.now().isoformat()
            }
            self.state['status_records'][project_id] = STATUS_HAS_DATA
        self.save_state()

    STAGE_TO_STATUS = {
        'validation': STATUS_INVALID_FORMAT,
        'no_run': STATUS_NO_RUN,
        'no_data': STATUS_NO_DATA,
    }

    def mark_failed(self, project_id, error, stage):
        with self.lock:
            self.state['failed_downloads'].append({
                'bioproject_id': project_id, 'stage': stage,
                'error': str(error), 'timestamp': datetime.now().isoformat()
            })
            self.state['status_records'][project_id] = self.STAGE_TO_STATUS.get(
                stage, STATUS_DOWNLOAD_ERROR)
        self.save_state()

    def mark_status(self, project_id, status):
        with self.lock:
            self.state['status_records'][project_id] = status
        self.save_state()

    def get_stats(self):
        with self.lock:
            return {
                'completed': len(self.state['completed_projects']),
                'failed': len(self.state['failed_downloads'])
            }

    def get_all_status(self):
        with self.lock:
            return dict(self.state['status_records'])


# ============================================================================
# BioProject Keyword Search (NCBI + CNCB)
# ============================================================================

class BioProjectDownloader:
    def __init__(self):
        Entrez.email = generate_fake_email()
        self.cncb_base_url = "https://ngdc.cncb.ac.cn"
        self.headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
            "Referer": f"{self.cncb_base_url}/search/"
        }

    def _setup_directories(self, output_dir):
        output_dir = Path(output_dir)
        dirs = {
            'root': output_dir,
            'ncbi_search': output_dir / 'tmp' / 'ncbi_keywords_search',
            'cncb_search': output_dir / 'tmp' / 'cncb_keywords_search',
            'results': output_dir / 'searched_keywords'
        }
        for path in dirs.values():
            path.mkdir(parents=True, exist_ok=True)
        return dirs

    def generate_search_queries(self, field, organism, opt=None):
        if opt:
            return [(f'("{f}") AND ("{o}") AND ("{op}")',
                     {'field': f, 'organism': o, 'opt': op})
                    for f, o, op in product(field, organism, opt)]
        return [(f'("{f}") AND ("{o}")',
                 {'field': f, 'organism': o, 'opt': None})
                for f, o in product(field, organism)]

    def search_and_download_batch(self, field, organism, opt=None,
                                   output_dir=Path("./downloads"),
                                   databases=None, file_format="tsv",
                                   delay=1.0):
        databases = databases or ["ncbi", "cncb"]
        dirs = self._setup_directories(output_dir)
        queries = self.generate_search_queries(field, organism, opt)

        print(f"Queries: {len(queries)} | Databases: {', '.join(databases)}\n")

        all_results = {"queries": [], "dirs": dirs}

        for idx, (query_string, keywords) in enumerate(queries, 1):
            print(f"[{idx}/{len(queries)}] {query_string}")
            query_dirname = self._sanitize_filename(query_string)
            results = {}

            if "ncbi" in databases:
                ncbi_dir = dirs['ncbi_search'] / query_dirname
                ncbi_dir.mkdir(parents=True, exist_ok=True)
                results["ncbi"] = self._download_from_ncbi(query_string, ncbi_dir)

            if "cncb" in databases:
                cncb_dir = dirs['cncb_search'] / query_dirname
                cncb_dir.mkdir(parents=True, exist_ok=True)
                results["cncb"] = self._download_from_cncb(query_string, cncb_dir, file_format)

            all_results["queries"].append({
                "query": query_string, "keywords": keywords, "files": results
            })

            if idx < len(queries):
                time.sleep(delay)

        self._save_summary(all_results, dirs['results'])
        return all_results

    def _download_from_ncbi(self, query, output_dir):
        try:
            search_handle = Entrez.esearch(db="bioproject", term=query, retmax=10000)
            search_results = Entrez.read(search_handle)
            search_handle.close()

            id_list = search_results["IdList"]
            if not id_list:
                print(f"  NCBI: 0")
                return None

            print(f"  NCBI: {len(id_list)}")
            output_file = output_dir / f"ncbi_{int(time.time())}.xml"
            fetch_handle = Entrez.efetch(db="bioproject", id=id_list,
                                         rettype="xml", retmode="xml")
            output_file.write_bytes(fetch_handle.read())
            fetch_handle.close()
            time.sleep(0.5)
            return output_file
        except Exception as e:
            print(f"  NCBI error: {e}")
            return None

    def _download_from_cncb(self, query, output_dir, file_format="tsv"):
        try:
            encoded_query = quote(query)
            api_url = (
                f"{self.cncb_base_url}/search/api/download/specific"
                f"?db=bioproject&q={encoded_query}&type={file_format}"
            )

            response = requests.get(api_url, headers=self.headers, timeout=120)
            response.raise_for_status()

            if len(response.content) < 100:
                print(f"  CNCB: 0 results")
                return None

            output_file = output_dir / f"cncb_{int(time.time())}.{file_format}"
            output_file.write_bytes(response.content)

            try:
                df = pd.read_csv(output_file, sep='\t' if file_format == 'tsv' else ',')
                print(f"  CNCB: {len(df)} results")
            except Exception:
                pass

            return output_file
        except Exception as e:
            print(f"  CNCB error: {e}")
            return None

    def _sanitize_filename(self, text, max_length=80):
        for char in '<>:"/\\|?*':
            text = text.replace(char, '_')
        text = text.replace(' ', '_').replace('"', '').replace("'", '')
        while '__' in text:
            text = text.replace('__', '_')
        return text[:max_length].strip('_')

    def _save_summary(self, all_results, results_dir):
        summary_file = results_dir / "search_summary.txt"
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(f"Total queries: {len(all_results['queries'])}\n")
            f.write(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            for idx, q in enumerate(all_results['queries'], 1):
                f.write(f"[{idx}] {q['query']}\n")
                for db, filepath in q['files'].items():
                    f.write(f"  {db.upper()}: {filepath.name if filepath else 'None'}\n")
        print(f"\n  Summary: {summary_file}\n")

    def parse_ncbi_xml(self, xml_file):
        from xml.etree import ElementTree as ET
        tree = ET.parse(xml_file)
        records = []
        for project in tree.getroot().findall('.//Project'):
            record = {}
            pid = project.find('.//ProjectID/ArchiveID')
            if pid is not None:
                record['accession'] = pid.get('accession')
            for tag, key in [('.//ProjectDescr/Title', 'title'),
                             ('.//ProjectDescr/Description', 'description'),
                             ('.//Organism/OrganismName', 'organism')]:
                el = project.find(tag)
                if el is not None:
                    record[key] = el.text
            records.append(record)
        return pd.DataFrame(records)

    def combine_batch_results(self, all_results, output_filename="combined_results.csv"):
        all_dfs = []
        print("Combining results...")

        for query_info in all_results['queries']:
            files = query_info.get('files', {})

            if files.get('ncbi') and files['ncbi'].exists():
                try:
                    df = self.parse_ncbi_xml(files['ncbi'])
                    if not df.empty:
                        df['source'] = 'NCBI'
                        all_dfs.append(df)
                except Exception:
                    pass

            if files.get('cncb') and files['cncb'].exists():
                try:
                    cncb_file = files['cncb']
                    df = pd.read_csv(cncb_file,
                                     sep='\t' if cncb_file.suffix == '.tsv' else ',',
                                     encoding='utf-8', on_bad_lines='skip', engine='python')
                    if not df.empty:
                        col_map = {'Accession': 'accession', 'Title': 'title',
                                   'Description': 'description', 'Species': 'organism'}
                        available = {k: v for k, v in col_map.items() if k in df.columns}
                        df = df[list(available.keys())].rename(columns=available)
                        df['source'] = 'CNCB'
                        all_dfs.append(df)
                        print(f"  Parsed CNCB: {len(df)} records")
                except Exception as e:
                    print(f"  CNCB parse error: {e}")

        if not all_dfs:
            print("No data found")
            return pd.DataFrame()

        combined_df = pd.concat(all_dfs, ignore_index=True)
        combined_df = combined_df.dropna(subset=['accession'])
        combined_df = combined_df[combined_df['accession'].str.startswith('PRJ', na=False)]
        combined_df = combined_df.drop_duplicates(subset=['accession'], keep='first')
        combined_df = combined_df[['accession', 'title', 'description', 'organism', 'source']]

        print(f"Final: {len(combined_df)} unique records")

        output_file = all_results['dirs']['results'] / output_filename
        combined_df.to_csv(output_file, index=False)

        ids_file = all_results['dirs']['results'] / "bioproject_ids.txt"
        combined_df['accession'].to_csv(ids_file, index=False, header=False)

        return combined_df


# ============================================================================
# CNCB/GSA Download
# ============================================================================

def download_cncb_metadata(accession, output_dir):
    """Download and process CNCB/GSA metadata for a BioProject."""
    BASE_URL = "https://ngdc.cncb.ac.cn/gsa"
    HEADERS = {"User-Agent": "Mozilla/5.0"}

    url = f"{BASE_URL}/search/getRunInfo"
    data = f'searchTerm=%26quot%3B{accession}%26quot%3BtotalDatas=9999%3BdownLoadCount=9999'

    try:
        resp = requests.post(url, data=data,
                             headers={**HEADERS, "Content-Type": "application/x-www-form-urlencoded"},
                             timeout=60)
        resp.raise_for_status()
        csv_content = resp.text
    except requests.RequestException:
        return None

    if csv_content.count('\n') < 2:
        return None

    (output_dir / f"{accession}.temp.csv").write_text(csv_content, encoding='utf-8')

    try:
        csv_df = pd.read_csv(StringIO(csv_content))
        if csv_df.empty:
            return None
        print(f"  CNCB CSV: {len(csv_df)} rows, {len(csv_df.columns)} columns")
    except Exception as e:
        print(f"  Failed to parse CSV: {e}")
        return None

    # Fetch Excel data for each CRA submission
    cra_ids = sorted({
        row.get('Submission', '').strip()
        for row in csv.DictReader(StringIO(csv_content))
        if row.get('Submission', '').strip().startswith('CRA')
    })

    all_excel_dfs = []
    if cra_ids:
        print(f"  Found {len(cra_ids)} CRA ID(s): {', '.join(cra_ids)}")
        for cra in cra_ids:
            try:
                resp = requests.post(f"{BASE_URL}/file/exportExcelFile",
                                     data={"type": "3", "dlAcession": cra},
                                     headers=HEADERS, timeout=60)
                resp.raise_for_status()

                temp_xlsx = output_dir / f"{cra}.temp.xlsx"
                temp_xlsx.write_bytes(resp.content)
                xlsx = pd.ExcelFile(temp_xlsx)

                sheet_dfs = []
                for sheet in xlsx.sheet_names:
                    sdf = pd.read_excel(xlsx, sheet_name=sheet)
                    sdf.columns = [f"{sheet}_{col}" for col in sdf.columns]
                    sheet_dfs.append(sdf)

                if sheet_dfs:
                    cra_merged = pd.concat(sheet_dfs, axis=1)
                    all_excel_dfs.append(cra_merged)
                    print(f"    {cra}: {len(cra_merged)} rows, {len(cra_merged.columns)} columns")

                temp_xlsx.unlink()
            except Exception as e:
                print(f"    {cra} failed: {e}")

    # Merge CSV + Excel
    if all_excel_dfs:
        excel_combined = pd.concat(all_excel_dfs, axis=0, ignore_index=True) \
            if len(all_excel_dfs) > 1 else all_excel_dfs[0]

        if len(csv_df) == len(excel_combined):
            final_df = pd.concat([csv_df, excel_combined], axis=1)
        else:
            print(f"  Row count mismatch (CSV: {len(csv_df)}, Excel: {len(excel_combined)}), using CSV only")
            final_df = csv_df
    else:
        final_df = csv_df

    return standardize_columns(final_df)


# ============================================================================
# NCBI Download
# ============================================================================

def get_biosamples_from_bioproject(bioproject_id):
    """Fetch BioSample UIDs linked to a BioProject."""

    def _fetch():
        search_handle = Entrez.esearch(db="bioproject", term=f"{bioproject_id}[Accession]", retmax=1)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if not search_results["IdList"]:
            return []

        link_handle = Entrez.elink(
            dbfrom="bioproject", db="biosample",
            id=search_results["IdList"][0], retmax=10000
        )
        link_results = Entrez.read(link_handle)
        link_handle.close()

        if link_results and link_results[0].get("LinkSetDb"):
            return [link["Id"] for link in link_results[0]["LinkSetDb"][0]["Link"]]
        return []

    try:
        return retry_wrapper(_fetch)
    except Exception:
        return []


def download_biosample_data(group_id, biosample_ids, output_dir):
    """Download BioSample metadata text."""
    if not biosample_ids:
        return None

    def _fetch():
        all_data = []
        for i in range(0, len(biosample_ids), BATCH_SIZE):
            batch = biosample_ids[i:i + BATCH_SIZE]
            handle = Entrez.efetch(db="biosample", id=",".join(batch),
                                   rettype="full", retmode="text")
            all_data.append(handle.read())
            handle.close()
            if len(biosample_ids) > BATCH_SIZE:
                time.sleep(DEFAULT_REQUEST_DELAY)

        output_file = output_dir / f"{group_id}_biosample.txt"
        output_file.write_text("\n".join(all_data), encoding="utf-8")
        return output_file

    try:
        return retry_wrapper(_fetch)
    except Exception:
        return None


def _fetch_sra_runinfo(sra_ids, group_id, output_dir):
    """Fetch SRA RunInfo CSV from a list of SRA UIDs. Shared by both flows."""
    if not sra_ids:
        return None

    def _fetch():
        all_data = []
        for i in range(0, len(sra_ids), BATCH_SIZE):
            batch = sra_ids[i:i + BATCH_SIZE]
            handle = Entrez.efetch(db="sra", id=",".join(batch),
                                   rettype="runinfo", retmode="text")
            data = handle.read()
            handle.close()
            if isinstance(data, bytes):
                data = data.decode('utf-8')
            # Skip header row for subsequent batches
            if i == 0:
                all_data.append(data)
            else:
                lines = data.split('\n')
                all_data.append('\n'.join(lines[1:]))
            if len(sra_ids) > BATCH_SIZE:
                time.sleep(DEFAULT_REQUEST_DELAY)

        output_file = output_dir / f"{group_id}_sra_runinfo.csv"
        output_file.write_text('\n'.join(all_data), encoding="utf-8")
        return output_file

    try:
        return retry_wrapper(_fetch)
    except Exception:
        return None


def download_sra_runinfo(bioproject_id, output_dir):
    """Search SRA by BioProject, then fetch RunInfo."""

    def _search():
        handle = Entrez.esearch(db="sra", term=f"{bioproject_id}[BioProject]", retmax=10000)
        results = Entrez.read(handle)
        handle.close()
        return results["IdList"]

    try:
        sra_ids = retry_wrapper(_search)
    except Exception:
        return None

    return _fetch_sra_runinfo(sra_ids, bioproject_id, output_dir)


def get_sra_ids_from_biosample_ids(biosample_accessions):
    """Link BioSample accessions -> SRA UIDs via esearch + elink."""

    sra_uids = []
    batch_size = 200

    for i in range(0, len(biosample_accessions), batch_size):
        batch = biosample_accessions[i:i + batch_size]

        def _fetch(b=batch):
            query = " OR ".join(f"{acc}[Accession]" for acc in b)
            handle = Entrez.esearch(db="biosample", term=query, retmax=10000)
            results = Entrez.read(handle)
            handle.close()
            bs_uids = results["IdList"]
            if not bs_uids:
                return []
            link_handle = Entrez.elink(dbfrom="biosample", db="sra",
                                       id=bs_uids, retmax=10000)
            link_results = Entrez.read(link_handle)
            link_handle.close()
            return [link["Id"]
                    for linkset in link_results
                    for linksetdb in linkset.get("LinkSetDb", [])
                    for link in linksetdb.get("Link", [])]

        try:
            sra_uids.extend(retry_wrapper(_fetch))
        except Exception:
            pass

        if len(biosample_accessions) > batch_size:
            time.sleep(DEFAULT_REQUEST_DELAY)

    return list(dict.fromkeys(sra_uids))


# ============================================================================
# BioSample Parser
# ============================================================================

def parse_biosample_file(file_path):
    """Parse BioSample text file into DataFrame."""
    if file_path is None or not Path(file_path).exists():
        return pd.DataFrame()

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception:
        return pd.DataFrame()

    samples = []
    for block in re.split(r'\n(?=\d+:\s+)', content):
        if not block.strip():
            continue

        data = {}
        match = re.search(r'Accession:\s+(\S+)', block)
        if match:
            data['BioSample'] = match.group(1)

        for key, pattern in [('Sample_Name', r'Identifiers:.*?Label:\s+(\S+)'),
                             ('Organism', r'Organism:\s+(.+?)(?:\n|$)')]:
            m = re.search(pattern, block)
            if m:
                data[key] = m.group(1).strip()

        attr_section = re.search(r'Attributes:(.+?)(?=\n\n|\Z)', block, re.DOTALL)
        if attr_section:
            for key, value in re.findall(r'/([^=]+?)="([^"]*)"', attr_section.group(1)):
                data[key.strip().replace(' ', '_').replace('-', '_')] = value

        if data:
            samples.append(data)

    return pd.DataFrame(samples)


# ============================================================================
# Data Merge
# ============================================================================

def merge_ncbi_data_single(biosample_df, sra_df):
    """Merge BioSample and SRA data."""
    if biosample_df.empty and sra_df.empty:
        return pd.DataFrame()
    if biosample_df.empty:
        return sra_df
    if sra_df.empty:
        return biosample_df

    if 'BioSample' in sra_df.columns and 'BioSample' in biosample_df.columns:
        return sra_df.merge(biosample_df, on='BioSample', how='outer',
                            suffixes=('', '_biosample'))

    return pd.concat([sra_df, biosample_df], axis=1)


def _download_and_merge_ncbi(group_id, biosample_ids, output_dir, sra_fetch_fn):
    """Shared logic: download BioSample + SRA, merge, standardize.

    sra_fetch_fn: callable that returns SRA RunInfo file path.
    """
    biosample_df = pd.DataFrame()
    sra_df = pd.DataFrame()

    try:
        biosample_file = download_biosample_data(group_id, biosample_ids,
                                                  output_dir)
        if biosample_file:
            biosample_df = parse_biosample_file(biosample_file)
    except Exception:
        pass

    time.sleep(DEFAULT_REQUEST_DELAY)

    try:
        sra_file = sra_fetch_fn()
        if sra_file:
            sra_df = pd.read_csv(sra_file)
    except Exception:
        pass

    if biosample_df.empty and sra_df.empty:
        return None

    return standardize_columns(merge_ncbi_data_single(biosample_df, sra_df))


def download_ncbi_metadata(bioproject_id, output_dir):
    """Download NCBI metadata for a BioProject."""
    biosample_ids = get_biosamples_from_bioproject(bioproject_id)
    return _download_and_merge_ncbi(
        bioproject_id, biosample_ids, output_dir,
        sra_fetch_fn=lambda: download_sra_runinfo(bioproject_id, output_dir)
    )


def download_ncbi_metadata_from_biosamples(biosample_accessions, output_dir):
    """Download NCBI metadata starting from BioSample IDs."""
    GROUP_ID = "BIOSAMPLE_INPUT"

    def _sra_fetch():
        sra_uids = get_sra_ids_from_biosample_ids(biosample_accessions)
        return _fetch_sra_runinfo(sra_uids, GROUP_ID, output_dir)

    result = _download_and_merge_ncbi(
        GROUP_ID, biosample_accessions, output_dir,
        sra_fetch_fn=_sra_fetch
    )

    # Restore original input BioSample IDs if NCBI cross-referenced them
    # (e.g., SAMD00518144 -> SAMN00518144)
    if result is not None and not result.empty and 'BioSample' in result.columns:
        original_set = set(biosample_accessions)
        # Build mapping: numeric suffix -> original accession
        suffix_to_original = {}
        for acc in biosample_accessions:
            m = re.match(r'SAM[A-Z]*(\d+)$', acc)
            if m:
                suffix_to_original[m.group(1)] = acc

        def _restore(val):
            if val in original_set:
                return val
            m = re.match(r'SAM[A-Z]*(\d+)$', str(val))
            if m and m.group(1) in suffix_to_original:
                return suffix_to_original[m.group(1)]
            return val

        result['BioSample'] = result['BioSample'].map(_restore)

    return result


# ============================================================================
# Single Project Processing
# ============================================================================

def process_single_bioproject(bioproject_id, output_dir, state_manager):
    """Process a single BioProject (NCBI or CNCB)."""
    print(f"\n{'─'*50}")
    print(f"Processing: {bioproject_id}")
    print('─'*50)

    if state_manager.is_completed(bioproject_id):
        print(f"  Already processed (checkpoint)")
        csv_path = output_dir / f"{bioproject_id}.processed.csv"
        try:
            return {'bioproject_id': bioproject_id,
                    'df': pd.read_csv(csv_path), 'source': 'checkpoint'}
        except Exception:
            pass

    if ID_PATTERNS['cncb'].match(bioproject_id):
        source = 'CNCB'
        df = download_cncb_metadata(bioproject_id, output_dir)
    elif ID_PATTERNS['ncbi'].match(bioproject_id):
        source = 'NCBI'
        df = download_ncbi_metadata(bioproject_id, output_dir)
    else:
        print(f"  Unknown format: {bioproject_id}")
        state_manager.mark_failed(bioproject_id, "Unknown format", 'validation')
        return None

    if df is None or df.empty:
        print(f"  No data retrieved")
        state_manager.mark_failed(bioproject_id, "No data", 'no_data')
        return None

    if not validate_run_data(df):
        print(f"  No valid Run information")
        state_manager.mark_failed(bioproject_id, "No Run info", 'no_run')
        return None

    # Preserve the user's original BioProject ID (NCBI may cross-reference
    # e.g., PRJEB -> PRJNA, PRJDB -> PRJNA)
    if 'BioProject' in df.columns:
        df['BioProject'] = bioproject_id

    df['Source_Database'] = source
    csv_path = output_dir / f"{bioproject_id}.processed.csv"
    df.to_csv(csv_path, index=False, encoding='utf-8')
    state_manager.mark_complete(bioproject_id, csv_path, len(df), source)
    print(f"  Saved: {csv_path.name}")

    return {'bioproject_id': bioproject_id, 'df': df, 'source': source}


def parallel_process_bioprojects(bioproject_list, output_dir, max_workers,
                                  state_manager):
    """Process multiple BioProjects in parallel."""
    print(f"\n{'='*70}")
    print(f"Processing {len(bioproject_list)} BioProjects (workers: {max_workers})")
    print('='*70)

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_single_bioproject, bid, output_dir,
                            state_manager): bid
            for bid in bioproject_list
        }
        for future in as_completed(futures):
            bid = futures[future]
            try:
                result = future.result()
                if result and result['df'] is not None:
                    results.append(result)
            except Exception as e:
                print(f"Task failed for {bid}: {e}")
                state_manager.mark_failed(bid, str(e), 'processing')

    print(f"\n{'='*70}")
    print(f"Processing Complete: {len(results)} / {len(bioproject_list)} successful")
    print('='*70)
    return results


# ============================================================================
# Status Report & Final Merge
# ============================================================================

def generate_status_report(all_ids, state_manager, output_dir):
    """Generate status.tsv."""
    status_records = state_manager.get_all_status()
    rows = [{'BioProject': bid, 'Status': status_records.get(bid, 'unknown')}
            for bid in all_ids]

    status_df = pd.DataFrame(rows)
    status_file = output_dir / "status.tsv"
    status_df.to_csv(status_file, sep='\t', index=False, encoding='utf-8')

    counts = status_df['Status'].value_counts().to_dict()
    print(f"\nStatus Report: {status_file}")
    for status, count in counts.items():
        print(f"  {status}: {count}")

    return status_df


def merge_all_results(results, output_dir, tmp_dir=None):
    """Merge all .processed.csv files into final output."""
    print(f"\n{'='*70}")
    print("Final Merge")
    print('='*70)

    scan_path = Path(tmp_dir) if tmp_dir else Path(output_dir)
    processed_files = sorted(scan_path.glob('*.processed.csv'))

    if not processed_files:
        print("  No .processed.csv files found")
        return pd.DataFrame()

    print(f"  Found {len(processed_files)} .processed.csv files")

    all_dfs = []
    for csv_file in processed_files:
        try:
            df = pd.read_csv(csv_file)
            if not df.empty and 'Run' in df.columns:
                all_dfs.append(df)
        except Exception as e:
            print(f"  Failed to read {csv_file.name}: {e}")

    if not all_dfs:
        print("  No valid DataFrames to merge")
        return pd.DataFrame()

    final_df = pd.concat(all_dfs, axis=0, ignore_index=True, sort=False)

    # Reorder: priority columns first
    first_cols = [c for c in PRIORITY_COLUMNS if c in final_df.columns]
    other_cols = [c for c in final_df.columns if c not in PRIORITY_COLUMNS]
    final_df = final_df[first_cols + other_cols]

    final_file = Path(output_dir) / "all_metadata_merged.csv"
    final_df.to_csv(final_file, index=False, encoding='utf-8-sig')

    print(f"  Total records: {len(final_df)}")
    print(f"  Output: {final_file}")
    print('='*70)

    return final_df


# ============================================================================
# Main Pipeline
# ============================================================================

def run_unified_pipeline(input_folder, output_folder, api_key=None, max_workers=None):
    """Execute unified metadata download pipeline."""
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)

    tmp_path = output_path / "tmp"
    tmp_path.mkdir(parents=True, exist_ok=True)

    state_manager = StateManager(tmp_path)

    setup_entrez(api_key)

    if max_workers is None:
        max_workers = 8 if api_key else DEFAULT_MAX_WORKERS

    print("\n" + "="*70)
    print("UNIFIED METADATA DOWNLOAD PIPELINE")
    print("="*70)
    print(f"Input:       {input_folder}")
    print(f"Output:      {output_folder}")
    print(f"API Key:     {'Yes' if api_key else 'No'}")
    print(f"Max Workers: {max_workers}")
    print("="*70)

    load_column_rename_dict()

    # Step 1: Read and classify input IDs
    print("\n[Step 1] Reading input IDs...")
    all_input_ids = read_input_ids(input_folder)
    if not all_input_ids:
        print("ERROR: No IDs found")
        return None

    groups = classify_ids(all_input_ids)

    print(f"\nTotal input IDs: {len(all_input_ids)}")
    print(f"  BioSample (SAM*):   {len(groups['biosample'])}")
    print(f"  CNCB (PRJC*):       {len(groups['cncb'])}")
    print(f"  NCBI (PRJ[EDN]*):   {len(groups['ncbi'])}")
    if groups['unknown']:
        print(f"  Unknown format:     {len(groups['unknown'])}")
        for oid in groups['unknown']:
            state_manager.mark_status(oid, STATUS_INVALID_FORMAT)

    stats = state_manager.get_stats()
    print(f"\nCheckpoint: {stats['completed']} completed, {stats['failed']} failed")

    # Step 1b: Process BioSample IDs
    biosample_results = []
    if groups['biosample']:
        print(f"\n[Step 1b] Processing {len(groups['biosample'])} BioSample IDs...")
        bs_df = download_ncbi_metadata_from_biosamples(
            groups['biosample'], tmp_path
        )
        if bs_df is not None and not bs_df.empty and validate_run_data(bs_df):
            bs_df['Source_Database'] = 'NCBI'
            csv_path = tmp_path / "BIOSAMPLE_INPUT.processed.csv"
            bs_df.to_csv(csv_path, index=False, encoding='utf-8')
            biosample_results.append({'bioproject_id': 'BIOSAMPLE_INPUT',
                                       'df': bs_df, 'source': 'NCBI'})
            print(f"  BioSample metadata saved: {len(bs_df)} records")
        else:
            print("  No valid data retrieved from BioSample IDs")

    # Step 2: Process BioProjects
    print("\n[Step 2] Processing BioProjects...")
    bioproject_ids = groups['cncb'] + groups['ncbi']

    results = parallel_process_bioprojects(
        bioproject_ids, tmp_path, max_workers, state_manager
    ) + biosample_results

    # Step 3: Status report
    print("\n[Step 3] Generating status report...")
    status_df = generate_status_report(all_input_ids, state_manager, output_path)

    # Step 4: Final merge
    print("\n[Step 4] Final merge...")
    final_df = merge_all_results(results, output_path, tmp_dir=tmp_path)

    # Summary
    print("\n" + "="*70)
    print("PIPELINE COMPLETE")
    print("="*70)
    print(f"Total input IDs: {len(all_input_ids)}")
    print(f"  With valid Run data: {len(results)}")
    print(f"  No data/No Run info: {len(all_input_ids) - len(results)}")
    print(f"\nOutput:")
    print(f"  status.tsv:             {len(status_df)} records")
    print(f"  all_metadata_merged.csv: {len(final_df)} records")
    print("="*70 + "\n")

    return {
        'input_ids': all_input_ids,
        'results': results,
        'final_df': final_df,
        'status_df': status_df
    }


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Unified Metadata Download Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-o', '--output', required=True,
                        help='Output directory')
    parser.add_argument('--keywords', action='store_true',
                        help='Enable keyword search mode')
    parser.add_argument('--field', nargs='+', default=None,
                        help='Search fields (required if --keywords)')
    parser.add_argument('--organism', nargs='+', default=None,
                        help='Organism terms (required if --keywords)')
    parser.add_argument('--opt', nargs='+', default=None,
                        help='Optional search terms')
    parser.add_argument('-i', '--input', default=None,
                        help='Input directory with ID txt files')
    parser.add_argument('-k', '--api-key', default=None,
                        help='NCBI API key')
    parser.add_argument('-w', '--max-workers', type=int, default=None,
                        help='Max parallel workers')

    args = parser.parse_args()

    if args.keywords:
        if not args.field or not args.organism:
            print("ERROR: --keywords requires both --field and --organism")
            parser.print_help()
            sys.exit(1)

        print("\n=== KEYWORD SEARCH MODE ===")
        print(f"Fields: {args.field}")
        print(f"Organisms: {args.organism}")
        if args.opt:
            print(f"Optional: {args.opt}")
        print("===========================\n")

        try:
            downloader = BioProjectDownloader()
            search_results = downloader.search_and_download_batch(
                field=args.field, organism=args.organism, opt=args.opt,
                output_dir=Path(args.output), databases=["ncbi", "cncb"], delay=1.0
            )

            combined_df = downloader.combine_batch_results(search_results)
            if combined_df.empty:
                print("ERROR: No BioProjects found")
                sys.exit(1)

            bioproject_input_folder = search_results['dirs']['results']
            result = run_unified_pipeline(
                str(bioproject_input_folder), args.output,
                args.api_key, args.max_workers
            )
            sys.exit(0 if result else 1)

        except KeyboardInterrupt:
            print("\n\nInterrupted by user")
            sys.exit(1)
        except Exception as e:
            print(f"\n\nERROR: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

    else:
        if not args.input:
            print("ERROR: --input required when not using --keywords mode")
            parser.print_help()
            sys.exit(1)

        if not os.path.isdir(args.input):
            print(f"ERROR: Input directory not found: {args.input}")
            sys.exit(1)

        print("\n=== INPUT MODE ===")
        print(f"Input: {args.input}")
        print("==================\n")

        try:
            result = run_unified_pipeline(
                args.input, args.output, args.api_key, args.max_workers
            )
            sys.exit(0 if result else 1)

        except KeyboardInterrupt:
            print("\n\nInterrupted by user")
            sys.exit(1)
        except Exception as e:
            print(f"\n\nERROR: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)


if __name__ == '__main__':
    main()
