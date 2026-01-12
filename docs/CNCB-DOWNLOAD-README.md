# CNCB/GSA元数据下载工具使用说明

## 快速开始

### 安装依赖

```bash
cd /Users/a1-6/Project-helsinki/Meta2Data

# 安装Python依赖
pip install -r requirements_gsa.txt

# 安装ChromeDriver (macOS)
brew install --cask chromedriver
```

### 基本使用

```bash
# 下载单个项目
python3 scripts/dl_metadata_from_CNCB.py CRA005036 ./metadata

# 支持PRJCA自动转换
python3 scripts/dl_metadata_from_CNCB.py PRJCA005036 ./metadata

# 批量下载
python3 scripts/dl_metadata_from_CNCB.py projects.txt ./metadata
```

## 文件说明

### 核心脚本

1. **dl_metadata_from_CNCB.py** - 主下载脚本
   - 支持单个/批量下载
   - 自动PRJCA→CRA转换
   - 失败记录和重试

2. **download_gsa_fastq.sh** - FASTQ数据下载
   - 通过FTP下载原始测序数据

### 依赖文件

- **requirements_gsa.txt** - Python依赖列表

## 使用示例

### 示例1: 下载单个项目

```bash
python3 scripts/dl_metadata_from_CNCB.py CRA005036 ./metadata
```

输出:
```
Found 1 projects to download
Output directory: ./metadata
Initializing Chrome browser...
[1/1] Processing CRA005036...
✓ Downloaded: ./metadata/CRA005036.xlsx (24658 bytes)
Total projects: 1
Successful: 1
Failed: 0
```

### 示例2: 批量下载

创建项目列表文件 `projects.txt`:
```
CRA005036
PRJCA004523
CRA003789
```

运行批量下载:
```bash
python3 scripts/dl_metadata_from_CNCB.py projects.txt ./metadata
```

### 示例3: 重新下载失败的项目

```bash
# 第一次下载后会生成 failed_downloads.txt
python3 scripts/dl_metadata_from_CNCB.py ./metadata/failed_downloads.txt ./metadata
```

## 函数说明

### download_from_bioproject_list()

主下载函数，可在Python脚本中导入使用：

```python
from scripts.dl_metadata_from_CNCB import download_from_bioproject_list

# 下载单个项目
success, failed = download_from_bioproject_list("CRA005036", "./metadata")

# 下载多个项目
projects = ["CRA005036", "PRJCA004523", "CRA003789"]
success, failed = download_from_bioproject_list(projects, "./metadata")

# 从文件下载
success, failed = download_from_bioproject_list("projects.txt", "./metadata")
```

参数:
- `bioproject_list`: 项目ID、ID列表或文件路径
- `output_dir`: 输出目录 (默认当前目录)

返回:
- `(成功数量, 失败数量)`

## 性能

- 单个项目: 3-5秒
- 100个项目: 5-10分钟
- 1000个项目: 1-2小时

## 故障排除

### ChromeDriver版本不匹配

```bash
# 更新ChromeDriver
brew upgrade chromedriver
```

### 下载失败

检查网络连接，使用失败列表重试:
```bash
python3 scripts/dl_metadata_from_CNCB.py ./metadata/failed_downloads.txt ./metadata
```

### Chrome未安装

```bash
brew install --cask google-chrome
brew install --cask chromedriver
```

## 配置Conda环境

```bash
# 创建环境
conda create -n gsa_download python=3.10

# 激活环境
conda activate gsa_download

# 安装依赖
pip install -r requirements_gsa.txt

# 安装ChromeDriver
brew install --cask chromedriver
```

## 与Meta2Data集成

下载的XLSX元数据可以与Meta2Data流程集成使用。

---

**最后更新:** 2025-12-09
