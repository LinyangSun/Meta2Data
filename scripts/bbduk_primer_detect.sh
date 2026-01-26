#!/bin/bash
#
# bbduk_primer_detect.sh - Frequency-based 16S Primer/Adapter Detection
#
# Based on the "Single-Unit Technical Extraction" principle:
# Technical sequences (Adapters + Linkers + Primers) form a unified block
# of high-frequency, low-entropy sequences at the 5' end.
#
# Usage:
#   Single-end: bbduk_primer_detect.sh -i input.fastq.gz -o output_dir
#   Paired-end: bbduk_primer_detect.sh -1 R1.fastq.gz -2 R2.fastq.gz -o output_dir
#
# Output:
#   - Cleaned FASTQ files in output_dir/
#   - Detection report: output_dir/primer_detection_report.txt
#   - Detected sequences: output_dir/detected_primers.fa

set -e

################################################################################
# Default Parameters
################################################################################

INPUT_R1=""
INPUT_R2=""
OUTPUT_DIR="bbduk_output"
THREADS=4
SAMPLE_SIZE=50000          # Number of reads to sample for detection
FREQ_THRESHOLD=1           # Frequency threshold (%)
PREFIX_LENGTH=40           # Length of 5' prefix to analyze
MIN_TECH_LENGTH=15         # Minimum technical sequence length to consider
BBDUK_K=15                 # BBDuk kmer length
BBDUK_MINK=11              # BBDuk minimum kmer length
BBDUK_HDIST=1              # Hamming distance tolerance

################################################################################
# Argument Parsing
################################################################################

show_help() {
    cat << EOF
Universal 16S Primer/Adapter Detection using Frequency Analysis

Usage:
  Single-end:  $0 -i <input.fastq.gz> -o <output_dir>
  Paired-end:  $0 -1 <R1.fastq.gz> -2 <R2.fastq.gz> -o <output_dir>

Required:
  -i FILE         Single-end input FASTQ (gzipped or not)
  -1 FILE         Paired-end R1 FASTQ (gzipped or not)
  -2 FILE         Paired-end R2 FASTQ (gzipped or not)
  -o DIR          Output directory

Optional:
  -t INT          Number of threads (default: 4)
  -s INT          Sample size for detection (default: 50000)
  -f INT          Frequency threshold in % (default: 1)
  -p INT          Prefix length to analyze (default: 40)
  -h              Show this help

EOF
}

while getopts "i:1:2:o:t:s:f:p:h" opt; do
    case $opt in
        i) INPUT_R1="$OPTARG" ;;
        1) INPUT_R1="$OPTARG" ;;
        2) INPUT_R2="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        s) SAMPLE_SIZE="$OPTARG" ;;
        f) FREQ_THRESHOLD="$OPTARG" ;;
        p) PREFIX_LENGTH="$OPTARG" ;;
        h) show_help; exit 0 ;;
        *) echo "Error: Invalid option -$OPTARG" >&2; show_help; exit 1 ;;
    esac
done

# Validate inputs
if [[ -z "$INPUT_R1" ]]; then
    echo "Error: Input file required (-i or -1)" >&2
    show_help
    exit 1
fi

if [[ ! -f "$INPUT_R1" ]]; then
    echo "Error: Input file not found: $INPUT_R1" >&2
    exit 1
fi

if [[ -n "$INPUT_R2" ]] && [[ ! -f "$INPUT_R2" ]]; then
    echo "Error: R2 file not found: $INPUT_R2" >&2
    exit 1
fi

# Check for bbduk.sh
if ! command -v bbduk.sh >/dev/null 2>&1; then
    echo "Error: bbduk.sh not found in PATH" >&2
    echo "Please install BBTools: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
REPORT="${OUTPUT_DIR}/primer_detection_report.txt"
DETECTED_FA="${OUTPUT_DIR}/detected_primers.fa"

# Use OUTPUT_DIR directly for temp files (don't create nested /temp)
TEMP_DIR="$OUTPUT_DIR"

################################################################################
# Step 1: Pre-clean Noise (Poly-G, Quality Trim)
################################################################################

echo "============================================================" | tee "$REPORT"
echo "16S PRIMER DETECTION - FREQUENCY-BASED METHOD" | tee -a "$REPORT"
echo "============================================================" | tee -a "$REPORT"
echo "" | tee -a "$REPORT"
echo "Strategy: K-mer frequency profiling to identify technical blocks" | tee -a "$REPORT"
echo "  1. Sample reads and extract 5' prefixes" | tee -a "$REPORT"
echo "  2. Find high-frequency sequences (>${FREQ_THRESHOLD}%)" | tee -a "$REPORT"
echo "  3. Use BBDuk to trim matching sequences" | tee -a "$REPORT"
echo "" | tee -a "$REPORT"

echo "[Step 1] Pre-cleaning: Removing poly-G tails and low-quality bases..." | tee -a "$REPORT"

if [[ -z "$INPUT_R2" ]]; then
    # Single-end
    bbduk.sh in="$INPUT_R1" out="${TEMP_DIR}/temp_clean.fq.gz" \
        trimpolyg=10 qtrim=rl trimq=10 minlen=50 threads="$THREADS" \
        2>&1 | grep -E "(Input:|Result:|Reads processed:)" | tee -a "$REPORT"
    CLEAN_R1="${TEMP_DIR}/temp_clean.fq.gz"
    CLEAN_R2=""
else
    # Paired-end
    bbduk.sh in="$INPUT_R1" in2="$INPUT_R2" \
        out="${TEMP_DIR}/temp_clean_R1.fq.gz" out2="${TEMP_DIR}/temp_clean_R2.fq.gz" \
        trimpolyg=10 qtrim=rl trimq=10 minlen=50 threads="$THREADS" \
        2>&1 | grep -E "(Input:|Result:|Reads processed:)" | tee -a "$REPORT"
    CLEAN_R1="${TEMP_DIR}/temp_clean_R1.fq.gz"
    CLEAN_R2="${TEMP_DIR}/temp_clean_R2.fq.gz"
fi

echo "" | tee -a "$REPORT"

################################################################################
# Step 2: Frequency Profiling - R1 (Forward)
################################################################################

echo "[Step 2] Profiling R1 5' ends for technical sequences..." | tee -a "$REPORT"
echo "  Sampling: First $SAMPLE_SIZE reads" | tee -a "$REPORT"
echo "  Prefix length: ${PREFIX_LENGTH}bp" | tee -a "$REPORT"
echo "  Frequency threshold: >${FREQ_THRESHOLD}%" | tee -a "$REPORT"
echo "" | tee -a "$REPORT"

# Convert FASTQ to FASTA and extract prefixes
SAMPLE_FASTA_R1="${TEMP_DIR}/sample_R1.fa"
zcat "$CLEAN_R1" 2>/dev/null || cat "$CLEAN_R1" | \
    head -n $((SAMPLE_SIZE * 4)) | \
    awk 'NR%4==1 {printf(">read%d\n", (NR-1)/4+1)} NR%4==2 {print}' > "$SAMPLE_FASTA_R1"

# Extract prefixes and find high-frequency sequences
LITERAL_REFS_R1=$(awk -v len="$PREFIX_LENGTH" '
    NR%2==0 {print substr($0,1,len)}
' "$SAMPLE_FASTA_R1" | \
    sort | uniq -c | sort -rn | \
    awk -v limit=$((SAMPLE_SIZE * FREQ_THRESHOLD / 100)) -v minlen="$MIN_TECH_LENGTH" \
        '$1 > limit && length($2) >= minlen {print $2}')

if [[ -z "$LITERAL_REFS_R1" ]] || [[ $(echo "$LITERAL_REFS_R1" | wc -l) -eq 0 ]]; then
    echo "  ✓ STATUS: No high-frequency technical sequences detected in R1" | tee -a "$REPORT"
    echo "  → Data appears CLEAN (primers already removed)" | tee -a "$REPORT"
    R1_NEEDS_TRIM=false
else
    echo "  ✓ Detected $(echo "$LITERAL_REFS_R1" | wc -l) high-frequency sequences in R1:" | tee -a "$REPORT"
    echo "$LITERAL_REFS_R1" | head -5 | nl -w2 -s'. ' | tee -a "$REPORT"
    if [[ $(echo "$LITERAL_REFS_R1" | wc -l) -gt 5 ]]; then
        echo "  ... and $(($(echo "$LITERAL_REFS_R1" | wc -l) - 5)) more" | tee -a "$REPORT"
    fi
    R1_NEEDS_TRIM=true

    # Save to FASTA
    echo "$LITERAL_REFS_R1" | awk '{print ">R1_tech_seq_"NR"\n"$0}' > "$DETECTED_FA"
fi

echo "" | tee -a "$REPORT"

################################################################################
# Step 3: Frequency Profiling - R2 (Reverse, if paired)
################################################################################

if [[ -n "$CLEAN_R2" ]]; then
    echo "[Step 3] Profiling R2 5' ends for technical sequences..." | tee -a "$REPORT"

    # Convert FASTQ to FASTA and extract prefixes
    SAMPLE_FASTA_R2="${TEMP_DIR}/sample_R2.fa"
    zcat "$CLEAN_R2" 2>/dev/null || cat "$CLEAN_R2" | \
        head -n $((SAMPLE_SIZE * 4)) | \
        awk 'NR%4==1 {printf(">read%d\n", (NR-1)/4+1)} NR%4==2 {print}' > "$SAMPLE_FASTA_R2"

    # Extract prefixes and find high-frequency sequences
    LITERAL_REFS_R2=$(awk -v len="$PREFIX_LENGTH" '
        NR%2==0 {print substr($0,1,len)}
    ' "$SAMPLE_FASTA_R2" | \
        sort | uniq -c | sort -rn | \
        awk -v limit=$((SAMPLE_SIZE * FREQ_THRESHOLD / 100)) -v minlen="$MIN_TECH_LENGTH" \
            '$1 > limit && length($2) >= minlen {print $2}')

    if [[ -z "$LITERAL_REFS_R2" ]] || [[ $(echo "$LITERAL_REFS_R2" | wc -l) -eq 0 ]]; then
        echo "  ✓ STATUS: No high-frequency technical sequences detected in R2" | tee -a "$REPORT"
        echo "  → R2 appears CLEAN" | tee -a "$REPORT"
        R2_NEEDS_TRIM=false
    else
        echo "  ✓ Detected $(echo "$LITERAL_REFS_R2" | wc -l) high-frequency sequences in R2:" | tee -a "$REPORT"
        echo "$LITERAL_REFS_R2" | head -5 | nl -w2 -s'. ' | tee -a "$REPORT"
        if [[ $(echo "$LITERAL_REFS_R2" | wc -l) -gt 5 ]]; then
            echo "  ... and $(($(echo "$LITERAL_REFS_R2" | wc -l) - 5)) more" | tee -a "$REPORT"
        fi
        R2_NEEDS_TRIM=true

        # Append to FASTA
        echo "$LITERAL_REFS_R2" | awk '{print ">R2_tech_seq_"NR"\n"$0}' >> "$DETECTED_FA"
    fi

    echo "" | tee -a "$REPORT"
fi

################################################################################
# Step 4: Execute BBDuk Trimming
################################################################################

echo "[Step 4] Executing primer/adapter removal..." | tee -a "$REPORT"

if [[ "$R1_NEEDS_TRIM" == false ]] && [[ "${R2_NEEDS_TRIM:-false}" == false ]]; then
    echo "  → No trimming needed. Copying files to output..." | tee -a "$REPORT"

    if [[ -z "$CLEAN_R2" ]]; then
        cp "$CLEAN_R1" "${OUTPUT_DIR}/$(basename "$INPUT_R1" .gz | sed 's/\.fastq/_clean.fastq/').gz"
    else
        cp "$CLEAN_R1" "${OUTPUT_DIR}/$(basename "$INPUT_R1" .gz | sed 's/\.fastq/_clean.fastq/').gz"
        cp "$CLEAN_R2" "${OUTPUT_DIR}/$(basename "$INPUT_R2" .gz | sed 's/\.fastq/_clean.fastq/').gz"
    fi
else
    echo "  → Running BBDuk with detected sequences..." | tee -a "$REPORT"

    if [[ -z "$CLEAN_R2" ]]; then
        # Single-end trimming
        bbduk.sh in="$CLEAN_R1" out="${OUTPUT_DIR}/$(basename "$INPUT_R1" .gz | sed 's/\.fastq/_clean.fastq/').gz" \
            ref="$DETECTED_FA" \
            ktrim=l k="$BBDUK_K" mink="$BBDUK_MINK" hdist="$BBDUK_HDIST" \
            restrictleft="$PREFIX_LENGTH" \
            minlen=50 threads="$THREADS" \
            stats="${OUTPUT_DIR}/bbduk_stats.txt" \
            2>&1 | grep -E "(Input:|Result:|Reads processed:|Trimmed:)" | tee -a "$REPORT"
    else
        # Paired-end trimming with overlap handling
        bbduk.sh in="$CLEAN_R1" in2="$CLEAN_R2" \
            out="${OUTPUT_DIR}/$(basename "$INPUT_R1" .gz | sed 's/\.fastq/_clean.fastq/').gz" \
            out2="${OUTPUT_DIR}/$(basename "$INPUT_R2" .gz | sed 's/\.fastq/_clean.fastq/').gz" \
            ref="$DETECTED_FA" \
            ktrim=l k="$BBDUK_K" mink="$BBDUK_MINK" hdist="$BBDUK_HDIST" \
            restrictleft="$PREFIX_LENGTH" \
            tpe tbo \
            minlen=50 threads="$THREADS" \
            stats="${OUTPUT_DIR}/bbduk_stats.txt" \
            2>&1 | grep -E "(Input:|Result:|Reads processed:|Trimmed:)" | tee -a "$REPORT"
    fi
fi

echo "" | tee -a "$REPORT"

################################################################################
# Cleanup and Summary
################################################################################

rm -rf "$TEMP_DIR"

echo "============================================================" | tee -a "$REPORT"
echo "PROCESSING COMPLETE" | tee -a "$REPORT"
echo "============================================================" | tee -a "$REPORT"
echo "" | tee -a "$REPORT"
echo "Output files:" | tee -a "$REPORT"
if [[ -z "$INPUT_R2" ]]; then
    echo "  - Cleaned: ${OUTPUT_DIR}/$(basename "$INPUT_R1" .gz | sed 's/\.fastq/_clean.fastq/').gz" | tee -a "$REPORT"
else
    echo "  - Cleaned R1: ${OUTPUT_DIR}/$(basename "$INPUT_R1" .gz | sed 's/\.fastq/_clean.fastq/').gz" | tee -a "$REPORT"
    echo "  - Cleaned R2: ${OUTPUT_DIR}/$(basename "$INPUT_R2" .gz | sed 's/\.fastq/_clean.fastq/').gz" | tee -a "$REPORT"
fi
echo "  - Report: $REPORT" | tee -a "$REPORT"
if [[ -f "$DETECTED_FA" ]]; then
    echo "  - Detected sequences: $DETECTED_FA" | tee -a "$REPORT"
fi
if [[ -f "${OUTPUT_DIR}/bbduk_stats.txt" ]]; then
    echo "  - BBDuk stats: ${OUTPUT_DIR}/bbduk_stats.txt" | tee -a "$REPORT"
fi
echo "" | tee -a "$REPORT"

exit 0
