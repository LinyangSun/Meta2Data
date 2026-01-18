from Bio import Entrez
import requests
import pandas as pd
from pathlib import Path
from urllib.parse import quote
import time
from typing import Optional, Dict, List, Tuple
from itertools import product

class BioProjectDownloader:
    """
    统一的 BioProject 下载器
    支持从 NCBI 和 CNCB 两个数据库搜索和下载数据
    """
    
    def __init__(self, email: str, custom_user_agent: Optional[str] = None):
        """初始化下载器"""
        Entrez.email = email
        self.cncb_base_url = "https://ngdc.cncb.ac.cn"
        
        user_agent = custom_user_agent or (
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/120.0.0.0 Safari/537.36"
        )
        
        self.headers = {
            "User-Agent": user_agent,
            "Referer": f"{self.cncb_base_url}/search/",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8"
        }
    
    def _setup_directories(self, output_dir: Path) -> Dict[str, Path]:
        """设置目录结构"""
        output_dir = Path(output_dir)
        
        dirs = {
            'root': output_dir,
            'tmp': output_dir / 'tmp',
            'ncbi_search': output_dir / 'tmp' / 'ncbi_keywords_search',
            'cncb_search': output_dir / 'tmp' / 'cncb_keywords_search',
            'results': output_dir / 'searched_keywords'
        }
        
        for dir_path in dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)
        
        return dirs
    
    def generate_search_queries(self,
                               field: List[str],
                               organism: List[str],
                               opt: Optional[List[str]] = None) -> List[Tuple[str, Dict[str, str]]]:
        """生成搜索关键词的所有组合"""
        queries = []
        
        if opt is None or len(opt) == 0:
            for f, o in product(field, organism):
                query_string = f'("{f}") AND ("{o}")'
                keywords = {'field': f, 'organism': o, 'opt': None}
                queries.append((query_string, keywords))
        else:
            for f, o, op in product(field, organism, opt):
                query_string = f'("{f}") AND ("{o}") AND ("{op}")'
                keywords = {'field': f, 'organism': o, 'opt': op}
                queries.append((query_string, keywords))
        
        return queries
    
    def search_and_download_batch(self,
                                  field: List[str],
                                  organism: List[str],
                                  opt: Optional[List[str]] = None,
                                  output_dir: Path = Path("./downloads"),
                                  databases: List[str] = ["ncbi", "cncb"],
                                  file_format: str = "tsv",
                                  delay: float = 1.0) -> Dict[str, List[Dict]]:
        """批量搜索并下载 BioProject 数据"""
        dirs = self._setup_directories(output_dir)
        queries = self.generate_search_queries(field, organism, opt)
        
        print(f"Output directory: {dirs['root']}")
        print(f"Total queries to execute: {len(queries)}")
        print(f"Databases: {', '.join(databases)}\n")
        
        all_results = {"queries": [], "results": [], "dirs": dirs}
        
        for idx, (query_string, keywords) in enumerate(queries, 1):
            print(f"[{idx}/{len(queries)}] Query: {query_string}")
            
            query_dirname = self._sanitize_filename(query_string)
            
            try:
                results = {}
                
                if "ncbi" in databases:
                    ncbi_query_dir = dirs['ncbi_search'] / query_dirname
                    ncbi_query_dir.mkdir(parents=True, exist_ok=True)
                    try:
                        results["ncbi"] = self._download_from_ncbi(
                            query_string, ncbi_query_dir
                        )
                    except Exception as e:
                        print(f"  NCBI failed: {e}")
                        results["ncbi"] = None
                
                if "cncb" in databases:
                    cncb_query_dir = dirs['cncb_search'] / query_dirname
                    cncb_query_dir.mkdir(parents=True, exist_ok=True)
                    try:
                        results["cncb"] = self._download_from_cncb(
                            query_string, cncb_query_dir, file_format
                        )
                    except Exception as e:
                        print(f"  CNCB failed: {e}")
                        results["cncb"] = None
                
                query_info = {
                    "query": query_string,
                    "keywords": keywords,
                    "files": results,
                    "ncbi_dir": dirs['ncbi_search'] / query_dirname if "ncbi" in databases else None,
                    "cncb_dir": dirs['cncb_search'] / query_dirname if "cncb" in databases else None
                }
                all_results["queries"].append(query_info)
                all_results["results"].append(results)
                
            except Exception as e:
                print(f"  ✗ Error: {e}")
                all_results["queries"].append({
                    "query": query_string,
                    "keywords": keywords,
                    "files": None,
                    "error": str(e)
                })
                all_results["results"].append(None)
            
            if idx < len(queries):
                time.sleep(delay)
        
        self._save_summary(all_results, dirs['results'])
        return all_results
    
    def _download_from_ncbi(self, query: str, output_dir: Path) -> Path:
        """从 NCBI 下载数据"""
        search_handle = Entrez.esearch(db="bioproject", term=query, retmax=10000)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_results["IdList"]
        count = int(search_results["Count"])
        
        if not id_list:
            print(f"  NCBI: 0 results")
            return None
        
        print(f"  NCBI: {count} results")
        
        timestamp = int(time.time())
        output_file = output_dir / f"ncbi_{timestamp}.xml"
        
        fetch_handle = Entrez.efetch(db="bioproject", id=id_list, 
                                     rettype="xml", retmode="xml")
        output_file.write_bytes(fetch_handle.read())
        fetch_handle.close()
        
        time.sleep(0.5)
        return output_file
    
    def _download_from_cncb(self, query: str, output_dir: Path, 
                           file_format: str = "tsv") -> Path:
        """从 CNCB 下载数据"""
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
        
        timestamp = int(time.time())
        output_file = output_dir / f"cncb_{timestamp}.{file_format}"
        output_file.write_bytes(response.content)
        
        try:
            df = pd.read_csv(output_file, sep='\t' if file_format == 'tsv' else ',')
            print(f"  CNCB: {len(df)} results")
        except:
            pass
        
        return output_file
    
    def _sanitize_filename(self, text: str, max_length: int = 80) -> str:
        """清理文件名"""
        for char in '<>:"/\\|?*':
            text = text.replace(char, '_')
        text = text.replace(' ', '_').replace('"', '').replace("'", '')
        while '__' in text:
            text = text.replace('__', '_')
        return text[:max_length].strip('_')
    
    def _save_summary(self, all_results: Dict, results_dir: Path):
        """保存搜索摘要"""
        summary_file = results_dir / "search_summary.txt"
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("BioProject Search Summary\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Total queries: {len(all_results['queries'])}\n")
            f.write(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            dirs = all_results.get('dirs', {})
            if dirs:
                f.write("Directory Structure:\n")
                f.write(f"  Root: {dirs['root']}\n")
                f.write(f"  NCBI searches: {dirs['ncbi_search']}\n")
                f.write(f"  CNCB searches: {dirs['cncb_search']}\n")
                f.write(f"  Results: {dirs['results']}\n\n")
            
            f.write("=" * 80 + "\n")
            f.write("Query Results:\n")
            f.write("=" * 80 + "\n")
            
            for idx, query_info in enumerate(all_results['queries'], 1):
                f.write(f"\n[{idx}] {query_info['query']}\n")
                f.write(f"    Keywords: {query_info['keywords']}\n")
                
                if 'error' in query_info:
                    f.write(f"    Status: FAILED - {query_info['error']}\n")
                else:
                    files = query_info.get('files', {})
                    if query_info.get('ncbi_dir'):
                        f.write(f"    NCBI dir: {query_info['ncbi_dir']}\n")
                    if query_info.get('cncb_dir'):
                        f.write(f"    CNCB dir: {query_info['cncb_dir']}\n")
                    
                    for db, filepath in files.items():
                        status = filepath.name if filepath else "No results"
                        f.write(f"    {db.upper()}: {status}\n")
        
        print(f"\n✓ Summary saved to: {summary_file}")
    
    def parse_ncbi_xml(self, xml_file: Path) -> pd.DataFrame:
        """解析 NCBI XML 为 DataFrame"""
        from xml.etree import ElementTree as ET
        
        tree = ET.parse(xml_file)
        root = tree.getroot()
        records = []
        
        for project in root.findall('.//Project'):
            record = {}
            
            project_id = project.find('.//ProjectID/ArchiveID')
            if project_id is not None:
                record['accession'] = project_id.get('accession')
            
            title = project.find('.//ProjectDescr/Title')
            if title is not None:
                record['title'] = title.text
            
            description = project.find('.//ProjectDescr/Description')
            if description is not None:
                record['description'] = description.text
            
            organism = project.find('.//Organism/OrganismName')
            if organism is not None:
                record['organism'] = organism.text
            
            records.append(record)
        
        return pd.DataFrame(records)
    
    def combine_batch_results(self, all_results: Dict, 
                             output_filename: str = "combined_results.csv") -> pd.DataFrame:
        """
        合并批量搜索的所有结果
        
        处理步骤:
        1. 合并所有 DataFrame
        2. 删除 accession 为空的行
        3. 删除 accession 不以 'PRJ' 开头的行
        4. 按 accession 去重
        """
        all_dfs = []
        
        print("\n" + "=" * 80)
        print("Combining results...")
        print("=" * 80)
        
        for query_idx, query_info in enumerate(all_results['queries'], 1):
            print(f"\n[Query {query_idx}] {query_info['query']}")
            
            if 'error' in query_info:
                continue
            
            files = query_info.get('files', {})
            
            # NCBI
            if files.get('ncbi') and files['ncbi'].exists():
                try:
                    df = self.parse_ncbi_xml(files['ncbi'])
                    if not df.empty:
                        df['source'] = 'NCBI'
                        all_dfs.append(df)
                        print(f"  ✓ NCBI: {len(df)} records")
                except Exception as e:
                    print(f"  ✗ NCBI failed: {e}")
            
            # CNCB
            if files.get('cncb') and files['cncb'].exists():
                try:
                    cncb_file = files['cncb']
                    df = pd.read_csv(cncb_file, 
                                   sep='\t' if cncb_file.suffix == '.tsv' else ',',
                                   encoding='utf-8')
                    
                    if not df.empty:
                        columns_to_keep = {
                            'Accession': 'accession',
                            'Title': 'title',
                            'Description': 'description',
                            'Species': 'organism'
                        }
                        
                        available_columns = {k: v for k, v in columns_to_keep.items() 
                                           if k in df.columns}
                        
                        df = df[list(available_columns.keys())].rename(columns=available_columns)
                        df['source'] = 'CNCB'
                        
                        all_dfs.append(df)
                        print(f"  ✓ CNCB: {len(df)} records")
                except Exception as e:
                    print(f"  ✗ CNCB failed: {e}")
        
        if not all_dfs:
            print("\n✗ No data to combine")
            return pd.DataFrame()
        
        # 合并
        print(f"\nCombining {len(all_dfs)} DataFrames...")
        combined_df = pd.concat(all_dfs, axis=0, ignore_index=True)
        
        print(f"Initial combined: {len(combined_df)} records")
        
        # ========== 数据清洗 ==========
        
        # 1. 删除 accession 为空的行
        before_empty = len(combined_df)
        combined_df = combined_df.dropna(subset=['accession'])
        removed_empty = before_empty - len(combined_df)
        if removed_empty > 0:
            print(f"  ✓ Removed {removed_empty} rows with empty accession")
        
        # 2. 删除 accession 不以 'PRJ' 开头的行
        before_prj = len(combined_df)
        combined_df = combined_df[combined_df['accession'].str.startswith('PRJ', na=False)]
        removed_non_prj = before_prj - len(combined_df)
        if removed_non_prj > 0:
            print(f"  ✓ Removed {removed_non_prj} rows with accession not starting with 'PRJ'")
        
        print(f"\nAfter filtering: {len(combined_df)} records")
        
        # 3. 按 accession 去重
        before_dedup = len(combined_df)
        combined_df = combined_df.drop_duplicates(subset=['accession'], keep='first')
        removed_duplicates = before_dedup - len(combined_df)
        if removed_duplicates > 0:
            print(f"  ✓ Removed {removed_duplicates} duplicate accessions")
        
        print(f"\nFinal dataset: {len(combined_df)} unique records")
        print(f"  NCBI: {len(combined_df[combined_df['source']=='NCBI'])}")
        print(f"  CNCB: {len(combined_df[combined_df['source']=='CNCB'])}")
        
        # 重新排列列
        column_order = ['accession', 'title', 'description', 'organism', 'source']
        combined_df = combined_df[column_order]
        
        # 保存
        dirs = all_results.get('dirs', {})
        if dirs and 'results' in dirs:
            output_file = dirs['results'] / output_filename
            combined_df.to_csv(output_file, index=False)
            print(f"\n✓ Saved to: {output_file}")
        
        return combined_df


# ============ 测试脚本 ============

def test_full_keywords():
    """完整测试"""
    
    print("=" * 80)
    print("BioProject Downloader - Full Test with Filtering")
    print("=" * 80)
    
    downloader = BioProjectDownloader(email="your.email@example.com")
    output_dir = Path("/Users/a1-6/Desktop/GetMetaFromNCBI/test1/")
    
    # 搜索
    results = downloader.search_and_download_batch(
        field=["microbiome", "microbiota", "bacteria"],
        organism=["bee", "pollinator"],
        opt=["gut"],
        output_dir=output_dir / "full_test_filtered",
        databases=["ncbi", "cncb"],
        delay=2.0
    )
    
    # 合并（包含过滤逻辑）
    combined_df = downloader.combine_batch_results(
        all_results=results,
        output_filename="combined_results.csv"
    )
    
    # 显示结果
    if not combined_df.empty:
        print("\n" + "=" * 80)
        print("Final Results Preview:")
        print("=" * 80)
        print(combined_df.head(10))
        
        # 验证所有 accession 都以 PRJ 开头
        non_prj = combined_df[~combined_df['accession'].str.startswith('PRJ', na=False)]
        if len(non_prj) > 0:
            print(f"\n⚠ Warning: Found {len(non_prj)} accessions not starting with PRJ")
        else:
            print(f"\n✓ All {len(combined_df)} accessions start with 'PRJ'")
        
        # 检查是否有空值
        empty_accessions = combined_df['accession'].isna().sum()
        if empty_accessions > 0:
            print(f"⚠ Warning: Found {empty_accessions} empty accessions")
        else:
            print(f"✓ No empty accessions")


if __name__ == "__main__":
    test_full_keywords()
