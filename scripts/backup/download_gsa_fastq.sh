#!/bin/bash
# GSA FASTQ数据批量下载脚本
# 使用方法: ./download_gsa_fastq.sh CRA005036 output_dir [sample_limit]

set -e

CRA_ID="$1"
OUTPUT_DIR="${2:-.}"
SAMPLE_LIMIT="${3:-all}"  # 默认下载所有样本，可以指定数量限制

if [ -z "$CRA_ID" ]; then
    echo "Usage: $0 <CRA_ID> [output_dir] [sample_limit]"
    echo "Example: $0 CRA005036 ./test 5"
    echo "         (downloads first 5 samples to ./test directory)"
    exit 1
fi

FTP_BASE="ftp://download.big.ac.cn/gsa"
CRA_URL="$FTP_BASE/$CRA_ID"

echo "========================================="
echo "GSA FASTQ Batch Downloader"
echo "========================================="
echo "CRA ID: $CRA_ID"
echo "Output: $OUTPUT_DIR"
echo "Sample limit: $SAMPLE_LIMIT"
echo ""

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 获取所有样本列表
echo "Fetching sample list..."
SAMPLES=$(curl -s --list-only "$CRA_URL/" | grep "^CRR")

if [ -z "$SAMPLES" ]; then
    echo "Error: No samples found for $CRA_ID"
    exit 1
fi

SAMPLE_COUNT=$(echo "$SAMPLES" | wc -l | tr -d ' ')
echo "Found $SAMPLE_COUNT samples"
echo ""

# 限制下载数量
if [ "$SAMPLE_LIMIT" != "all" ]; then
    SAMPLES=$(echo "$SAMPLES" | head -n "$SAMPLE_LIMIT")
    echo "Limiting to first $SAMPLE_LIMIT samples"
    echo ""
fi

# 创建元数据CSV文件
METADATA_FILE="$OUTPUT_DIR/${CRA_ID}_samples.csv"
echo "CRA_ID,CRR_ID,Forward_File,Reverse_File" > "$METADATA_FILE"

# 下载每个样本
COUNT=0
for SAMPLE in $SAMPLES; do
    COUNT=$((COUNT + 1))
    echo "[$COUNT/$SAMPLE_COUNT] Processing: $SAMPLE"

    SAMPLE_DIR="$OUTPUT_DIR/$SAMPLE"
    mkdir -p "$SAMPLE_DIR"

    # 获取该样本的文件列表
    FILES=$(curl -s --list-only "$CRA_URL/$SAMPLE/")

    FORWARD=""
    REVERSE=""

    # 下载FASTQ文件
    for FILE in $FILES; do
        if [[ $FILE == *.gz ]]; then
            echo "  Downloading: $FILE"
            curl -s "$CRA_URL/$SAMPLE/$FILE" -o "$SAMPLE_DIR/$FILE"

            # 记录forward和reverse文件
            if [[ $FILE == *_f1.gz ]] || [[ $FILE == *_1.fq.gz ]] || [[ $FILE == *_R1.fastq.gz ]]; then
                FORWARD="$SAMPLE_DIR/$FILE"
            elif [[ $FILE == *_r2.gz ]] || [[ $FILE == *_2.fq.gz ]] || [[ $FILE == *_R2.fastq.gz ]]; then
                REVERSE="$SAMPLE_DIR/$FILE"
            fi

            # 显示文件大小
            SIZE=$(ls -lh "$SAMPLE_DIR/$FILE" | awk '{print $5}')
            echo "    → $FILE ($SIZE)"
        fi
    done

    # 记录到元数据文件
    echo "$CRA_ID,$SAMPLE,$FORWARD,$REVERSE" >> "$METADATA_FILE"
    echo ""
done

echo "========================================="
echo "Download completed!"
echo "Total samples downloaded: $COUNT"
echo "Metadata saved to: $METADATA_FILE"
echo "========================================="
echo ""
echo "File structure:"
ls -lh "$OUTPUT_DIR" | head -20
