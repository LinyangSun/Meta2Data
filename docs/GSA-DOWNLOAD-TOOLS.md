# GSA数据下载工具集

中国国家基因库(CNCB/GSA)数据下载完整解决方案

## 快速开始

### 批量下载元数据（推荐）

```bash
# 1. 创建CRA列表文件
cat > cra_list.txt << EOF
CRA005036
CRA004523
CRA003789
EOF

# 2. 批量下载
python3 scripts/batch_download_gsa_metadata.py cra_list.txt ./metadata

# 3. 检查结果
ls -lh ./metadata/*.xlsx
```

## 工具总览

| 工具 | 功能 | 稳定性 | 推荐度 | 用途 |
|------|------|--------|--------|------|
| `batch_download_gsa_metadata.py` | 批量下载元数据(XLSX) | ⭐⭐⭐⭐⭐ | ✅ **强烈推荐** | 大规模批量下载 |
| `download_gsa_selenium.py` | 单个元数据下载 | ⭐⭐⭐⭐⭐ | ✅ 推荐 | 单个项目下载 |
| `download_gsa_fastq.sh` | 下载原始FASTQ数据 | ⭐⭐⭐⭐ | ✅ 推荐 | 下载测序数据 |
| `generate_gsa_metadata_from_ftp.py` | FTP生成元数据CSV | ⭐⭐⭐ | 备选 | 备用方案 |
| `download_gsa_metadata.py` | API下载（Python） | ⭐ | ❌ 不推荐 | API不稳定 |
| `download_gsa_direct.sh` | API下载（Bash） | ⭐ | ❌ 不推荐 | API不稳定 |

## 详细说明

### 1. batch_download_gsa_metadata.py ⭐⭐⭐⭐⭐

**最佳选择 - 批量下载元数据**

```bash
python3 scripts/batch_download_gsa_metadata.py cra_list.txt ./output_dir
```

**特点:**
- ✅ 自动化浏览器（Selenium）
- ✅ 批量处理
- ✅ 断点续传（跳过已下载）
- ✅ 错误处理和统计
- ✅ 失败列表保存
- ✅ 下载完整XLSX元数据

**依赖:**
- `pip3 install selenium`
- `brew install --cask chromedriver`

**参数:**
```
python3 batch_download_gsa_metadata.py <cra_list_file> [output_dir]
```

**示例:**
```bash
# 基本用法
python3 scripts/batch_download_gsa_metadata.py my_cra_list.txt ./metadata

# 重新下载失败的
python3 scripts/batch_download_gsa_metadata.py ./metadata/failed_downloads.txt ./metadata
```

---

### 2. download_gsa_selenium.py ⭐⭐⭐⭐⭐

**单个项目下载**

```bash
python3 scripts/download_gsa_selenium.py CRA005036 ./output_dir
```

**特点:**
- ✅ 简单直接
- ✅ 下载单个CRA
- ✅ 自动化浏览器
- ✅ 快速测试

**依赖:** 同上

**参数:**
```
python3 download_gsa_selenium.py <CRA_ID> [output_dir]
```

---

### 3. download_gsa_fastq.sh ⭐⭐⭐⭐

**下载原始测序数据（FASTQ文件）**

```bash
./scripts/download_gsa_fastq.sh CRA005036 ./data 10
```

**特点:**
- ✅ 通过FTP下载（稳定）
- ✅ 支持限制样本数量
- ✅ 自动生成样本CSV
- ✅ 显示进度

**依赖:** 仅需 `curl`

**参数:**
```bash
./download_gsa_fastq.sh <CRA_ID> [output_dir] [sample_limit]
```

**示例:**
```bash
# 下载前5个样本
./scripts/download_gsa_fastq.sh CRA005036 ./test 5

# 下载所有样本
./scripts/download_gsa_fastq.sh CRA005036 ./data all
```

---

### 4. generate_gsa_metadata_from_ftp.py ⭐⭐⭐

**从FTP生成基础元数据CSV**

```bash
python3 scripts/generate_gsa_metadata_from_ftp.py CRA005036 ./output_dir
```

**特点:**
- ✅ 不依赖API
- ✅ 通过FTP稳定获取
- ⚠️  信息较少（基础信息）
- ⚠️  生成CSV非XLSX

**依赖:** 仅需 `curl`

**用途:** 当Selenium不可用时的备选方案

---

### 5. download_gsa_metadata.py ⭐ (不推荐)

**通过API下载 - 不稳定**

```bash
python3 scripts/download_gsa_metadata.py CRA005036
```

**问题:**
- ❌ API经常超时
- ❌ 连接频繁中断
- ❌ 下载成功率低

**仅用于:** 学习和测试

---

### 6. download_gsa_direct.sh ⭐ (不推荐)

**Bash版API下载 - 不稳定**

```bash
./scripts/download_gsa_direct.sh CRA005036 ./output
```

**问题:** 同上

**仅用于:** 学习和测试

---

## 使用场景

### 场景1: 我需要下载100个项目的元数据

```bash
# 创建列表
cat > projects.txt << EOF
CRA001
CRA002
CRA003
...
CRA100
EOF

# 批量下载
python3 scripts/batch_download_gsa_metadata.py projects.txt ./metadata
```

### 场景2: 我只需要一个项目的元数据

```bash
python3 scripts/download_gsa_selenium.py CRA005036 ./
```

### 场景3: 我需要下载原始测序数据

```bash
# 下载前10个样本的FASTQ文件
./scripts/download_gsa_fastq.sh CRA005036 ./data 10
```

### 场景4: 我只需要知道有哪些样本

```bash
# 生成样本列表CSV
python3 scripts/generate_gsa_metadata_from_ftp.py CRA005036 ./
```

### 场景5: Selenium无法使用

```bash
# 方案A: 使用FTP方案
python3 scripts/generate_gsa_metadata_from_ftp.py CRA005036 ./

# 方案B: 手动下载
# 1. 访问 https://ngdc.cncb.ac.cn/gsa/browse/CRA005036
# 2. 点击下载 CRA005036.xlsx
```

## 环境要求

### 最小要求（FTP方案）
```bash
# 只需要curl
curl --version
```

### 推荐配置（Selenium方案）
```bash
# Python 3.x
python3 --version

# Selenium
pip3 install selenium

# ChromeDriver
brew install --cask chromedriver

# 验证安装
chromedriver --version
```

## 性能对比

| 方法 | 100个项目 | 1000个项目 | 稳定性 |
|------|-----------|-----------|---------|
| Selenium批量 | 5-10分钟 | 1-2小时 | 高 |
| FTP生成CSV | 10-15分钟 | 2-3小时 | 中 |
| API直接下载 | ❌ 失败率高 | ❌ 不可行 | 极低 |

## 故障排除

### Q: ChromeDriver报错？
```bash
# 更新ChromeDriver
brew upgrade chromedriver

# 或手动下载匹配版本
# https://chromedriver.chromium.org/downloads
```

### Q: 下载失败怎么办？
```bash
# 使用失败列表重新下载
python3 scripts/batch_download_gsa_metadata.py failed_downloads.txt ./metadata
```

### Q: 网络不稳定？
```bash
# 方法1: 分批下载（每次20-50个）
# 方法2: 使用VPN连接中国节点
# 方法3: 选择夜间或周末下载
```

### Q: 没有Chrome浏览器？
```bash
# 安装Chrome
brew install --cask google-chrome

# 然后安装ChromeDriver
brew install --cask chromedriver
```

## 完整工作流

### 步骤1: 准备项目列表
```bash
# 从PRJCA转换到CRA
# PRJCA005036 → CRA005036

cat > my_projects.txt << EOF
CRA005036
CRA004523
CRA003789
EOF
```

### 步骤2: 批量下载元数据
```bash
mkdir -p ./metadata_downloads
python3 scripts/batch_download_gsa_metadata.py my_projects.txt ./metadata_downloads
```

### 步骤3: 检查结果
```bash
ls -lh ./metadata_downloads/
cat ./metadata_downloads/failed_downloads.txt  # 如果有失败的
```

### 步骤4: 重试失败项目
```bash
python3 scripts/batch_download_gsa_metadata.py \
    ./metadata_downloads/failed_downloads.txt \
    ./metadata_downloads
```

### 步骤5: 验证数据
```bash
# 检查文件大小
find ./metadata_downloads -name "*.xlsx" -size 0 -ls

# 统计
echo "Total downloaded: $(ls ./metadata_downloads/*.xlsx 2>/dev/null | wc -l)"
```

## 技术原理

### Selenium方案（推荐）
1. 启动无头Chrome浏览器
2. 访问GSA browse页面
3. 模拟点击下载按钮
4. 等待文件下载完成
5. 验证文件完整性

**优点:** 绕过API限制，模拟真实用户行为，稳定性高

### FTP方案（备选）
1. 通过FTP列出样本目录
2. 解析文件名识别配对关系
3. 生成基础元数据CSV
4. 可选下载FASTQ文件

**优点:** 不依赖浏览器，轻量级

### API方案（已弃用）
1. 直接调用GSA API
2. 下载XLSX/JSON

**问题:** API极不稳定，经常超时

## 相关文档

- [GSA-BATCH-DOWNLOAD.md](GSA-BATCH-DOWNLOAD.md) - 批量下载详细教程
- [GSA-METADATA-DOWNLOAD.md](GSA-METADATA-DOWNLOAD.md) - 元数据下载指南
- [CLAUDE.MD](../CLAUDE.MD) - Meta2Data项目总览

## 支持

如有问题，请检查：
1. ChromeDriver版本是否匹配Chrome版本
2. 网络连接是否稳定
3. GSA网站是否正常访问
4. 失败列表中的项目是否真实存在

---

**最后更新:** 2025-12-09
