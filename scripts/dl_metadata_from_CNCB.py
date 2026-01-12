#!/usr/bin/env python3
"""
从CNCB/GSA下载元数据
支持单个和批量下载
使用Selenium自动化浏览器
"""
import time
import os
import sys
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options


def download_single_metadata(cra_id, download_dir, driver):
    """
    下载单个CRA项目的元数据

    Args:
        cra_id: CRA项目ID (如 CRA005036)
        download_dir: 下载目录路径
        driver: Selenium WebDriver实例

    Returns:
        bool: 下载是否成功
    """
    print(f"\n{'='*60}")
    print(f"Downloading: {cra_id}")
    print('='*60)

    expected_file = os.path.join(download_dir, f"{cra_id}.xlsx")

    # 如果文件已存在且大小>0，跳过
    if os.path.exists(expected_file) and os.path.getsize(expected_file) > 0:
        print(f"✓ Already exists: {expected_file} ({os.path.getsize(expected_file)} bytes)")
        return True

    try:
        # 访问GSA页面
        url = f"https://ngdc.cncb.ac.cn/gsa/browse/{cra_id}"
        print(f"Navigating to: {url}")
        driver.get(url)

        # 等待页面加载
        time.sleep(2)

        # 查找下载链接
        print("Looking for download button...")
        selectors = [
            f"//a[contains(text(), '{cra_id}.xlsx')]",
            "//a[contains(@href, 'exportExcelFile')]",
        ]

        element = None
        for selector in selectors:
            try:
                element = driver.find_element(By.XPATH, selector)
                if element:
                    break
            except:
                continue

        if not element:
            print(f"✗ Could not find download button for {cra_id}")
            return False

        # 点击下载
        print("Clicking download...")
        element.click()

        # 等待下载完成
        print("Waiting for download...")
        for i in range(60):
            if os.path.exists(expected_file) and os.path.getsize(expected_file) > 0:
                print(f"✓ Downloaded: {expected_file} ({os.path.getsize(expected_file)} bytes)")
                return True
            time.sleep(1)

        print(f"✗ Download timeout for {cra_id}")
        return False

    except Exception as e:
        print(f"✗ Error downloading {cra_id}: {e}")
        return False


def init_driver(download_dir):
    """
    初始化Selenium WebDriver

    Args:
        download_dir: 下载目录路径

    Returns:
        WebDriver实例
    """
    download_dir = os.path.abspath(download_dir)
    os.makedirs(download_dir, exist_ok=True)

    chrome_options = Options()
    chrome_options.add_argument('--headless')  # 无头模式
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--disable-dev-shm-usage')

    # 设置下载目录
    prefs = {
        "download.default_directory": download_dir,
        "download.prompt_for_download": False,
        "download.directory_upgrade": True,
        "safebrowsing.enabled": True
    }
    chrome_options.add_experimental_option("prefs", prefs)

    print("Initializing Chrome browser...")
    driver = webdriver.Chrome(options=chrome_options)
    return driver


def download_from_bioproject_list(bioproject_list, output_dir="."):
    """
    从BioProject/CRA ID列表批量下载元数据

    Args:
        bioproject_list: BioProject ID列表或文件路径
        output_dir: 输出目录

    Returns:
        tuple: (成功数量, 失败数量)
    """
    # 如果是文件路径，读取文件
    if isinstance(bioproject_list, str) and os.path.isfile(bioproject_list):
        with open(bioproject_list, 'r') as f:
            ids = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    elif isinstance(bioproject_list, (list, tuple)):
        ids = bioproject_list
    else:
        ids = [bioproject_list]

    # 转换PRJCA到CRA
    cra_ids = []
    for id_str in ids:
        if id_str.startswith('PRJCA'):
            cra_id = id_str.replace('PRJCA', 'CRA')
            print(f"Converting {id_str} → {cra_id}")
            cra_ids.append(cra_id)
        elif id_str.startswith('CRA'):
            cra_ids.append(id_str)
        else:
            print(f"Warning: Unknown ID format: {id_str}")
            cra_ids.append(id_str)

    print(f"\nFound {len(cra_ids)} projects to download")
    print(f"Output directory: {output_dir}")
    print()

    # 初始化浏览器
    driver = init_driver(output_dir)

    # 统计
    success_count = 0
    failed_list = []

    try:
        for idx, cra_id in enumerate(cra_ids, 1):
            print(f"\n[{idx}/{len(cra_ids)}] Processing {cra_id}...")

            success = download_single_metadata(cra_id, output_dir, driver)

            if success:
                success_count += 1
            else:
                failed_list.append(cra_id)

            # 每10个下载后稍作休息
            if idx % 10 == 0:
                print("\nTaking a short break...")
                time.sleep(5)

    finally:
        driver.quit()

    # 输出统计
    print("\n" + "="*60)
    print("Download Summary")
    print("="*60)
    print(f"Total projects: {len(cra_ids)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {len(failed_list)}")

    if failed_list:
        print("\nFailed projects:")
        for cra in failed_list:
            print(f"  - {cra}")

        # 保存失败列表
        failed_file = os.path.join(output_dir, "failed_downloads.txt")
        with open(failed_file, 'w') as f:
            for cra in failed_list:
                f.write(f"{cra}\n")
        print(f"\nFailed list saved to: {failed_file}")

    print()
    return success_count, len(failed_list)


def main():
    """命令行入口"""
    if len(sys.argv) < 2:
        print("Usage: python3 dl_metadata_from_CNCB.py <bioproject_id|list_file> [output_dir]")
        print()
        print("Examples:")
        print("  # Download single project")
        print("  python3 dl_metadata_from_CNCB.py CRA005036 ./metadata")
        print("  python3 dl_metadata_from_CNCB.py PRJCA005036 ./metadata")
        print()
        print("  # Download from list file")
        print("  python3 dl_metadata_from_CNCB.py projects.txt ./metadata")
        print()
        print("List file format (one ID per line):")
        print("  CRA005036")
        print("  PRJCA004523")
        print("  CRA000613")
        print("  # Comments start with #")
        sys.exit(1)

    input_arg = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "."

    success, failed = download_from_bioproject_list(input_arg, output_dir)

    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
