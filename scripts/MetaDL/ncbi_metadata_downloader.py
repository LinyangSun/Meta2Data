def get_sequencing_platform(srr_id):
    """
    通过SRA accession获取测序平台信息
    """
    from Bio import Entrez
    import xml.etree.ElementTree as ET
    
    # 设置你的email（NCBI要求）
    Entrez.email = "your_email@example.com"
    
    try:
        # 搜索获取UID
        search_handle = Entrez.esearch(db="sra", term=srr_id)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        if not search_results['IdList']:
            return None
        
        uid = search_results['IdList'][0]
        
        # 获取详细信息
        fetch_handle = Entrez.efetch(db="sra", id=uid, retmode="xml")
        xml_data = fetch_handle.read()
        fetch_handle.close()
        
        # 解析XML
        root = ET.fromstring(xml_data)
        platform = root.find('.//PLATFORM')
        
        if platform is not None and len(platform) > 0:
            return platform[0].tag
        
        return None
        
    except Exception as e:
        print(f"Error querying {srr_id}: {e}")
        return None


# 使用示例
if __name__ == "__main__":
    test_sras = ['SRR35007785', 'SRR36832824', 'SRR32509170', 'SRR20818414']
    
    for sra in test_sras:
        platform = get_sequencing_platform(sra)
        print(f"{sra}: {platform}")