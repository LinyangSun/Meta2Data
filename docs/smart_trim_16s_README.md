# smart_trim_16s - 智能16S引物检测与切除

## 功能简介

`smart_trim_16s` 是一个自动化的 16S rRNA 基因引物检测和切除工具。它通过以下方式智能判断序列是否包含引物：

1. **频率分析**: 统计样本前20bp序列的出现频率
2. **坐标验证**: 使用BLAST将高频序列比对到 E. coli 16S 参考序列
3. **锚点校验**: 验证比对位置是否落在已知保守区（V1-V9区域）
4. **智能处理**: 仅当确认为引物时才执行切除，否则直接分发原始数据

## 为什么需要这个工具？

在16S测序数据中，有时前20bp的高频序列可能是：
- **真正的引物**: 来自PCR扩增的保守区引物
- **高丰度菌种**: 样本中某个优势菌株的保守序列

传统方法会误将高丰度菌种序列当作引物切除，导致数据损失。`smart_trim_16s` 通过坐标验证避免这个问题。

## 安装依赖

### 必需软件

```bash
# BLAST (用于序列比对)
conda install -c bioconda blast

# cutadapt (用于序列切除)
conda install -c bioconda cutadapt

# BioPython (用于序列解析)
pip install biopython
```

### 参考文件

需要 E. coli 16S rRNA 基因参考序列：
- 默认路径: `./Meta2Data/docs/J01859.1.fna`
- GenBank ID: J01859.1
- 序列长度: 1542 bp (完整16S基因)

## 使用方法

### 基本用法

```bash
python3 py_16s.py smart_trim_16s \
    --input_path /path/to/fastq_dir \
    --output_path /path/to/output_dir
```

### 完整参数

```bash
python3 py_16s.py smart_trim_16s \
    --input_path /path/to/fastq_dir \      # 输入FASTQ目录
    --output_path /path/to/output_dir \    # 输出目录
    --tmp_path /tmp \                       # 临时文件目录 (可选)
    --ref_path ./Meta2Data/docs/J01859.1.fna  # E.coli参考序列 (可选)
```

### 参数说明

| 参数 | 必需 | 默认值 | 说明 |
|------|------|--------|------|
| `--input_path` | ✓ | - | 包含 *_R1*.fastq* 和 *_R2*.fastq* 文件的目录 |
| `--output_path` | ✓ | - | 输出目录，会自动创建 |
| `--tmp_path` | ✗ | `/tmp` | 临时文件存储目录 |
| `--ref_path` | ✗ | `./Meta2Data/docs/J01859.1.fna` | E. coli 16S 参考序列路径 |

## 工作流程

### 1. 初始化
- 检查参考文件是否存在
- 构建 BLAST 数据库（如果不存在）
- 扫描输入目录中的FASTQ文件

### 2. 标杆样本探测
- 读取第一个样本文件的前1000条序列
- 统计前20bp序列的出现频率
- 提取出现频率最高的序列

### 3. 坐标验证（仅当高频序列占比 > 50%）

使用 BLAST 将高频序列比对到 E. coli 16S 参考序列：

```bash
blastn -query top_candidate.fa \
       -db J01859.1.fna \
       -outfmt '6 sstart' \
       -task blastn-short
```

获取比对位置（sstart），并检查是否落在已知锚点附近（±20bp容错）：

**合法锚点列表** (基于E. coli 16S基因坐标):
- **V1区**: 8, 27
- **V3区**: 338, 341
- **V4区**: 515, 518, 534
- **V5区**: 785, 799, 806
- **V6-V7区**: 907, 926, 1046
- **V9区**: 1099, 1100
- **常用反向引物**: 1391, 1492

### 4. 决策与处理

#### 情况A: 确认为引物
```
高频序列占比 > 50% ✓
BLAST 比对成功 ✓
坐标落在合法锚点附近 ✓
→ 执行切除
```

使用 cutadapt 切除双端前20bp：
```bash
cutadapt -u 20 -U 20 \
    -o output_R1.fastq \
    -p output_R2.fastq \
    input_R1.fastq input_R2.fastq
```

#### 情况B: 非引物序列
```
高频序列占比 < 50% ✗
或 BLAST 比对失败 ✗
或 坐标不在合法锚点 ✗
→ 软链接分发
```

直接创建软链接，不修改原始数据：
```bash
ln -s /path/to/input/*.fastq /path/to/output/
```

## 输出说明

### 日志信息

```
[INFO] --- 标杆样本探测: SampleA_R1.fastq.gz ---
[INFO] 高频前20bp序列: GTGCCAGCMGCCGCGGTAA (出现 850/1000 次)
[INFO] 【验证通过】坐标 515 属于合法保守区。确认为 16S 引物。
[INFO] 检测到有效引物信号。启动全量物理切除 (双端前 20bp)...
[INFO] ✓ 成功处理 24 对样本文件
[INFO] 处理完成。结果存储于: /path/to/output
```

### 输出文件

**检测到引物时**:
- 所有样本的 R1/R2 文件都会被切除前20bp
- 文件名保持不变
- FASTQ 格式保持不变（包括压缩格式）

**未检测到引物时**:
- 所有文件通过软链接指向原始文件
- 不占用额外磁盘空间
- 保持原始数据完整性

## 示例

### 示例1: 标准 V4 区引物

```bash
# 输入文件
input/
├── Sample1_R1.fastq.gz  # 前20bp: GTGCCAGCMGCCGCGGTAA (515F引物)
├── Sample1_R2.fastq.gz
├── Sample2_R1.fastq.gz
└── Sample2_R2.fastq.gz

# 运行命令
python3 py_16s.py smart_trim_16s \
    --input_path input/ \
    --output_path output/

# 结果
[INFO] 【验证通过】坐标 515 属于合法保守区。确认为 16S 引物。
[INFO] ✓ 成功处理 2 对样本文件

# 输出文件（已切除前20bp）
output/
├── Sample1_R1.fastq.gz
├── Sample1_R2.fastq.gz
├── Sample2_R1.fastq.gz
└── Sample2_R2.fastq.gz
```

### 示例2: 高丰度菌种误判

```bash
# 输入文件
input/
├── Sample1_R1.fastq.gz  # 前20bp: ACGTACGTACGTACGTACGT (某优势菌保守区)
├── Sample1_R2.fastq.gz

# 运行命令
python3 py_16s.py smart_trim_16s \
    --input_path input/ \
    --output_path output/

# 结果
[INFO] 【验证拦截】坐标 1250 属于高变区内部。判定为高丰度菌种序列，非引物。
[INFO] ✓ 创建 2 个软链接

# 输出文件（软链接，保持原始数据）
output/
├── Sample1_R1.fastq.gz -> /full/path/to/input/Sample1_R1.fastq.gz
└── Sample1_R2.fastq.gz -> /full/path/to/input/Sample1_R2.fastq.gz
```

## 常见问题

### Q1: 为什么我的引物没有被检测到？

**A**: 可能原因：
1. 引物出现频率 < 50% → 样本中可能包含无引物的序列
2. 引物序列与 E. coli 16S 差异较大 → BLAST 比对失败
3. 使用了非标准引物位点 → 不在预定义的锚点列表中

**解决方法**: 检查日志中的频率和坐标信息，必要时手动添加锚点

### Q2: 可以用于非16S的数据吗？

**A**: 不推荐。该工具专门针对16S rRNA基因设计，使用 E. coli 16S作为参考。对于ITS、18S等其他标记基因，需要修改参考序列和锚点列表。

### Q3: 为什么固定切除20bp？

**A**: 20bp 是大多数16S引物的标准长度（如515F, 806R等）。如果你的引物长度不同，需要修改源代码中的 `cutadapt -u 20 -U 20` 参数。

### Q4: 支持单端测序数据吗？

**A**: 当前版本仅支持双端测序（paired-end）。单端数据会被跳过并警告。

## 技术细节

### 锚点设计原理

锚点列表基于以下原则设计：
1. **保守区边界**: V1-V9 各区域的起始/结束位置
2. **常用引物**: 文献中高频使用的引物位点
3. **容错范围**: ±20bp 允许不同研究中的坐标差异

### BLAST 参数说明

```bash
-task blastn-short   # 短序列优化模式（适用于 < 50bp）
-outfmt '6 sstart'   # 仅输出目标序列起始坐标
```

### 性能优化

- 仅分析第一个样本的前1000条序列（约2秒）
- BLAST 数据库只构建一次，缓存复用
- cutadapt 自动多线程并行处理

## 集成到Pipeline

在 AmpliconPIP 中集成：

```bash
# 在 Illumina 流程的 PrimerDetection 步骤前调用
cd $dataset_path
python3 py_16s.py smart_trim_16s \
    --input_path ori_fastq/ \
    --output_path temp/step_01_trimmed/ \
    --ref_path ${SCRIPTS}/../docs/J01859.1.fna
```

## 更新日志

- **v1.0** (2025-01-14): 初始版本，支持双端测序数据的智能引物检测与切除

## 相关工具

- **cutadapt**: 序列adapter切除工具
- **BLAST+**: 序列比对工具套件
- **BioPython**: Python生物信息学库

## 引用

如果使用此工具，请引用：

```
Meta2Data Project (2025). smart_trim_16s: Intelligent primer detection 
and trimming for 16S rRNA gene sequencing data.
```

