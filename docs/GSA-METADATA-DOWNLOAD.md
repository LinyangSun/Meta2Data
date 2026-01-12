# GSA (中国国家基因库) 元数据下载指南

## 问题说明

CNCB/GSA数据库的元数据导出API (`https://ngdc.cncb.ac.cn/gsa/file/exportExcelFile`) 经常不稳定，连接会被服务器中断。

## 手动下载方法（推荐）

### 方法1: 通过网页下载

1. 访问GSA项目页面：`https://ngdc.cncb.ac.cn/gsa/browse/<CRA_ID>`
   - 例如：https://ngdc.cncb.ac.cn/gsa/browse/CRA005036

2. 在页面上找到"元数据信息"部分

3. 点击 `<CRA_ID>.xlsx` 链接下载元数据文件

4. 将下载的文件保存到您的项目目录

### 方法2: 使用浏览器开发者工具捕获下载链接

1. 打开浏览器开发者工具 (F12)
2. 切换到 Network (网络) 标签
3. 访问 GSA browse 页面
4. 点击元数据下载按钮
5. 在 Network 标签中找到实际的下载请求
6. 右键该请求 → Copy → Copy as cURL
7. 在终端粘贴并执行该 cURL 命令

### 方法3: 从FTP目录下载样本列表

虽然无法直接获取XLSX格式的元数据，但可以从FTP列出所有样本：

```bash
# 列出所有样本ID
curl -s --list-only "ftp://download.big.ac.cn/gsa/CRA005036/" | grep "^CRR"

# 输出示例：
# CRR325340
# CRR325341
# CRR325342
# ...
```

## 自动化下载（部分功能）

我们提供了几个脚本尝试自动下载：

### Python脚本 (不稳定)

```bash
cd /Users/a1-6/Project-helsinki/Meta2Data
python3 scripts/download_gsa_metadata.py CRA005036
```

**问题**: API经常返回空响应或连接中断

### Bash脚本 (不稳定)

```bash
./scripts/download_gsa_direct.sh CRA005036 ./test
```

**问题**: 同样依赖不稳定的API

## 下载原始测序数据

如果需要下载FASTQ文件（原始测序数据）：

```bash
# 下载指定数量的样本
./scripts/download_gsa_fastq.sh CRA005036 ./output 5

# 下载所有样本
./scripts/download_gsa_fastq.sh CRA005036 ./output all
```

## 元数据文件格式

GSA元数据XLSX文件通常包含以下信息：

- **Run信息**: CRR编号、样本名称、测序平台
- **样本信息**: Biosample ID、样本描述
- **测序信息**: Library策略、Layout (single/paired)
- **数据文件**: FASTQ文件名和下载链接

## 替代方案

如果CNCB网站持续无法访问，可以考虑：

1. **使用NCBI SRA数据库**:
   - 很多GSA项目也会同步到NCBI SRA
   - 使用 `iseq` 工具下载 (更稳定)
   - 例如: `iseq -i PRJNA661613 -m`

2. **联系数据提供者**:
   - 直接向论文作者索要元数据

3. **使用Aspera高速下载** (如果可用):
   ```bash
   ascp -P33001 -i <key_file> -QT -l100m -k1 -d \
     aspera01@download.cncb.ac.cn:gsa/CRA005036 ./
   ```

## 网络建议

- GSA/CNCB服务器位于中国，从中国大陆访问会更稳定
- 如果在海外访问，建议：
  - 多次重试
  - 使用VPN连接到中国大陆节点
  - 选择访问流量较少的时间段

## 总结

**最可靠的方法**: 手动通过网页下载元数据XLSX文件

由于GSA的API不稳定，目前无法完全自动化下载元数据，建议您：
1. 手动下载元数据XLSX文件
2. 使用我们的脚本自动下载FASTQ数据（如需要）
