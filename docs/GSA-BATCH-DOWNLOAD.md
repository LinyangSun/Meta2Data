# GSA元数据批量下载指南

## 问题解决方案

通过使用 **Selenium自动化浏览器**，我们成功解决了GSA API不稳定的问题，实现了大规模批量下载。

## 前置要求

### 1. 安装Selenium
```bash
pip3 install selenium
```

### 2. 安装ChromeDriver
```bash
# macOS
brew install --cask chromedriver

# 或者手动下载
# https://chromedriver.chromium.org/downloads
```

## 使用方法

### 方法1: 单个CRA下载

```bash
cd /Users/a1-6/Project-helsinki/Meta2Data

python3 scripts/download_gsa_selenium.py CRA005036 ./test
```

**输出:**
```
Downloading metadata for CRA005036...
Initializing Chrome browser...
Navigating to: https://ngdc.cncb.ac.cn/gsa/browse/CRA005036
Looking for download button...
Clicking download button...
Waiting for download to complete...
✓ Downloaded: ./test/CRA005036.xlsx (24658 bytes)
```

### 方法2: 批量下载（推荐）

#### 步骤1: 创建CRA列表文件

创建一个文本文件 `cra_list.txt`，每行一个CRA ID：

```txt
# GSA项目列表
CRA005036
CRA004523
CRA003789
CRA002156
# 更多CRA ID...
```

#### 步骤2: 运行批量下载

```bash
python3 scripts/batch_download_gsa_metadata.py cra_list.txt ./metadata_output
```

**输出示例:**
```
Found 4 CRA IDs to download
Output directory: ./metadata_output

Initializing Chrome browser...

[1/4] Processing CRA005036...
============================================================
Downloading: CRA005036
============================================================
✓ Downloaded: ./metadata_output/CRA005036.xlsx (24656 bytes)

[2/4] Processing CRA004523...
============================================================
Downloading: CRA004523
============================================================
✓ Downloaded: ./metadata_output/CRA004523.xlsx (18234 bytes)

...

============================================================
Download Summary
============================================================
Total CRAs: 4
Successful: 3
Failed: 1

Failed CRAs:
  - CRA002156

Failed list saved to: ./metadata_output/failed_downloads.txt
```

## 功能特点

1. **自动跳过已下载文件**: 如果文件已存在且大小>0，自动跳过
2. **失败重试**: 可以使用失败列表文件重新下载失败的项目
3. **进度显示**: 实时显示下载进度和状态
4. **无头模式**: 后台运行，不打开浏览器窗口
5. **批量休息**: 每10个下载后暂停5秒，避免请求过于频繁

## 高级选项

### 重新下载失败的项目

```bash
# 第一次批量下载后，会生成 failed_downloads.txt
# 使用这个文件重新下载失败的项目
python3 scripts/batch_download_gsa_metadata.py ./metadata_output/failed_downloads.txt ./metadata_output
```

### 从PRJCA转换到CRA

如果您有PRJCA编号，通常对应的CRA编号是：
- PRJCA000613 → CRA000613
- PRJCA004523 → CRA004523
- PRJCA005036 → CRA005036

只需将 `PRJCA` 替换为 `CRA` 即可。

### 查看下载的元数据

```bash
# 使用Excel或LibreOffice打开
open CRA005036.xlsx

# 或者在Python中处理
python3 -c "
import pandas as pd
df = pd.read_excel('CRA005036.xlsx')
print(df.head())
"
```

## 性能估计

- **下载速度**: 约每个CRA 3-5秒
- **100个CRA项目**: 约5-10分钟
- **1000个CRA项目**: 约1-2小时

## 故障排除

### 问题1: ChromeDriver版本不匹配

**错误:** `This version of ChromeDriver only supports Chrome version...`

**解决:**
```bash
# 更新ChromeDriver
brew upgrade chromedriver

# 或下载匹配的版本
# https://chromedriver.chromium.org/downloads
```

### 问题2: 下载超时

**原因:** 网络连接不稳定或GSA服务器繁忙

**解决:** 使用失败列表重新下载：
```bash
python3 scripts/batch_download_gsa_metadata.py failed_downloads.txt ./metadata_output
```

### 问题3: Chrome未安装

**错误:** `selenium.common.exceptions.WebDriverException: Message: 'chromedriver' executable needs to be in PATH`

**解决:**
```bash
# 安装Chrome浏览器
brew install --cask google-chrome

# 安装ChromeDriver
brew install --cask chromedriver
```

## 脚本说明

### 1. download_gsa_selenium.py
- **功能**: 下载单个CRA的元数据
- **优点**: 简单直接
- **用途**: 测试或下载少量项目

### 2. batch_download_gsa_metadata.py
- **功能**: 批量下载多个CRA的元数据
- **优点**: 自动化、可恢复、统计报告
- **用途**: 大规模数据采集

### 3. generate_gsa_metadata_from_ftp.py
- **功能**: 从FTP生成基础元数据CSV
- **优点**: 不需要浏览器，仅依赖FTP
- **限制**: 信息较少，不如XLSX完整
- **用途**: 备用方案

## 最佳实践

1. **分批下载**: 建议每次下载100-200个项目
2. **定期保存**: 使用版本控制跟踪下载的元数据
3. **验证文件**: 下载后检查文件大小和格式
4. **网络时段**: 选择网络稳定的时间段进行大批量下载

## 与Meta2Data流程集成

下载元数据后，可以使用Meta2Data的AmpliconPIP流程处理：

```bash
# 1. 下载元数据
python3 scripts/batch_download_gsa_metadata.py cra_list.txt ./metadata

# 2. 转换XLSX到CSV（如果需要）
# Meta2Data可能需要CSV格式

# 3. 运行AmpliconPIP
./bin/Meta2Data AmpliconPIP -m metadata/CRA005036_converted.csv -t 8
```

## 总结

**成功的批量下载方案:**
- ✅ 使用Selenium自动化浏览器
- ✅ 绕过不稳定的GSA API
- ✅ 支持大规模批量下载
- ✅ 自动错误处理和重试
- ✅ 完整的元数据信息（XLSX格式）

**不推荐的方法:**
- ❌ 直接调用GSA API（极不稳定）
- ❌ 手动逐个下载（效率太低）
