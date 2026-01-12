# Changelog: Migration from Dataset-ID to BioProject ID

## 变更日期: 2025-12-09

## 变更概述

将流程从使用自定义的 `Datasets-ID` 列改为使用 `Data-Bioproject` 列作为主要标识符。

## 变更原因

1. **标准化**: BioProject ID是公共数据库的标准标识符
2. **唯一性**: BioProject ID全球唯一，避免命名冲突
3. **可追溯性**: 可以直接通过BioProject ID追溯到原始数据源
4. **兼容性**: 支持多种数据库格式（NCBI, ENA, DDBJ, GSA/CNCB）

## 支持的BioProject格式

- **PRJNA** + 数字: NCBI SRA (美国)
- **PRJEB** + 数字: ENA (欧洲)
- **PRJDB** + 数字: DDBJ (日本)
- **PRJCA** + 数字: GSA/CNCB (中国)

示例:
- PRJNA123456
- PRJEB567890
- PRJDB234567
- PRJCA012345

## 主要变更

### 1. MetaProcessingRun.sh

#### 变更前:
```bash
# 使用 Datasets-ID 列
py_16s.py GenerateDatasetsIDsFile \
    --FilePath "$metafile" \
    --Datasets_ID 'Datasets-ID' \
    --SequencingPlatform "Data-SequencingPlatform"

# 查找 wild-* 文件夹
all_folders=("${BASE_DIR}"/wild-*/)
```

#### 变更后:
```bash
# 使用 Data-Bioproject 列
py_16s.py GenerateDatasetsIDsFile \
    --FilePath "$metafile" \
    --Datasets_ID 'Data-Bioproject' \
    --SequencingPlatform "Data-SequencingPlatform"

# 查找 PRJ* 文件夹 (支持所有BioProject格式)
all_folders=()
for folder in "${BASE_DIR}"/PRJ*/; do
    if [ -d "$folder" ]; then
        all_folders+=("$folder")
    fi
done
```

### 2. 文件夹命名

#### 变更前:
```
BASE_DIR/
├── wild-dataset1/
├── wild-dataset2/
└── wild-dataset3/
```

#### 变更后:
```
BASE_DIR/
├── PRJNA123456/
├── PRJEB567890/
├── PRJCA012345/
└── PRJDB234567/
```

### 3. 文件命名

所有输出文件现在使用BioProject ID命名:

- `PRJNA123456_sra.txt`
- `PRJNA123456-final-rep-seqs.qza`
- `PRJNA123456-final-table.qza`

## CSV格式要求

### 必需列:
- `Data-Bioproject`: BioProject ID (现在用作主要标识符)
- `Data-SequencingPlatform`: 测序平台 (illumina/454/pacbio)
- `Data-SRA`: SRA运行编号
- `Data-Biosample`: BioSample ID
- `Data-Forward`: 正向引物
- `Data-Reverse`: 反向引物
- `Data-Region`: 目标区域

### 不再使用的列:
- ~~`Datasets-ID`~~: 已被 `Data-Bioproject` 替代

### CSV示例:
```csv
Data-Bioproject,Data-SequencingPlatform,Data-SRA,Data-Biosample,Data-Forward,Data-Reverse,Data-Region
PRJNA123456,illumina,SRR789012,SAMN123456,GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT,V4
PRJNA123456,illumina,SRR789013,SAMN123457,GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT,V4
PRJEB567890,illumina,ERR345678,SAMEA123456,CCTACGGGNGGCWGCAG,GACTACHVGGGTATCTAATCC,V3-V4
PRJCA012345,illumina,CRR123456,SAMC123456,ACTCCTACGGGAGGCAGCAG,GGACTACHVGGGTWTCTAAT,V3-V4
```

## 迁移步骤

如果您有旧的数据和CSV文件：

### 1. 更新CSV文件
```bash
# 如果原CSV有 Datasets-ID 列，需要确保也有 Data-Bioproject 列
# Data-Bioproject 列将成为新的标识符
```

### 2. 重命名文件夹
```bash
# 将旧的 wild-* 文件夹重命名为 BioProject ID
cd /your/base/directory
mv wild-dataset1 PRJNA123456
mv wild-dataset2 PRJEB567890
```

### 3. 重命名文件
```bash
# 在每个文件夹内
cd PRJNA123456
mv wild-dataset1_sra.txt PRJNA123456_sra.txt
mv wild-dataset1-final-rep-seqs.qza PRJNA123456-final-rep-seqs.qza
mv wild-dataset1-final-table.qza PRJNA123456-final-table.qza
```

## 向后兼容性

- ✅ `datasets_ID.txt` 文件名保持不变（为了兼容性）
- ✅ 文件内容现在包含BioProject ID而非自定义Dataset ID
- ✅ Python脚本 (`py_16s.py`) 的 `--Datasets_ID` 参数现在指向 `Data-Bioproject` 列

## 测试验证

运行流程前，验证：

```bash
# 1. 检查BioProject文件夹是否存在
ls -d PRJ*/

# 2. 验证CSV文件包含 Data-Bioproject 列
head -1 metadata.csv | grep "Data-Bioproject"

# 3. 检查输出文件命名
find . -name "*-final-*.qza"
```

## 错误排查

### 错误: "No BioProject folders (PRJ*) found"
**原因**: 文件夹命名不符合BioProject格式
**解决**: 确保文件夹名为 PRJNA*/PRJEB*/PRJCA*/PRJDB* 格式

### 错误: "No valid datasets were collected"
**原因**: 文件夹内缺少 `*-final-rep-seqs.qza` 或 `*-final-table.qza`
**解决**:
1. 检查流程是否成功完成
2. 验证输出文件命名是否正确
3. 确保文件使用BioProject ID命名

## 文档更新

相关文档已更新:
- ✅ CLAUDE.MD - 已说明使用BioProject作为标识符
- ✅ MetaProcessingRun.sh - 代码注释已更新
- ✅ 本CHANGELOG文档

---

**维护者**: Claude AI
**审核**: 待用户确认
