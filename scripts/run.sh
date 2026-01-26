#!/bin/bash
set -e

################################################################################
#                          MULTI-VSEARCH PIPELINE                              #
################################################################################
#
# Purpose: Complete multi-platform amplicon data processing pipeline
# Phases: 1) Dataset preparation, 2) Smart trimming + 454 processing, 3) Merging
#
################################################################################

# Find scripts directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS="${SCRIPT_DIR}/scripts"

show_help() {
    cat << EOF
Usage: multi-vsearch.sh [options]

Complete multi-platform amplicon processing pipeline with smart primer trimming.

Required options:
    -m, --metadata FILE        Input metadata CSV file
    -o, --output DIR           Output base directory

Optional options:
    -t, --threads INT          Number of CPU threads (default: 4)
    --col-bioproject NAME      Column for BioProject/Dataset ID (default: 'Data-Bioproject')
    --col-sra NAME             Column for SRA accession (default: 'Data-SRA')
    -h, --help                 Show this help message

EOF
}

################################################################################
#                          DEFAULT PARAMETERS                                  #
################################################################################

METADATA=""
OUTPUT=""
THREADS=4
COL_BIOPROJECT="Data-Bioproject"
COL_SRA="Data-SRA"

################################################################################
#                          ARGUMENT PARSING                                    #
################################################################################

while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--metadata) METADATA="$2"; shift 2 ;;
        -o|--output) OUTPUT="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --col-bioproject) COL_BIOPROJECT="$2"; shift 2 ;;
        --col-sra) COL_SRA="$2"; shift 2 ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Error: Unknown option '$1'"; show_help; exit 1 ;;
    esac
done

# Validate input
if [[ -z "$METADATA" ]]; then
    echo "Error: --metadata is required"
    exit 1
fi

if [[ "$METADATA" != /* ]]; then
    echo "Error: --metadata must be a full absolute path (starting with '/')"
    exit 1
fi

if [[ ! -f "$METADATA" ]]; then
    echo "Error: Metadata file '$METADATA' does not exist"
    exit 1
fi

if [[ -z "$OUTPUT" ]]; then
    OUTPUT=$(dirname "$METADATA")
fi

export cpu=$THREADS

# Source function library
if [[ -f "${SCRIPTS}/Function_Import.sh" ]]; then
    source "${SCRIPTS}/Function_Import.sh"
else
    echo "Error: Function_Import.sh not found"
    exit 1
fi

cd "$OUTPUT" || exit 1
mkdir -p "${OUTPUT}/analysis/"

################################################################################
#                          PHASE 1: DATASET PREPARATION                        #
################################################################################

echo "========================================="
echo "PHASE 1: Dataset Preparation"
echo "Started: $(date)"
echo "========================================="

# Step 1: Generate dataset ID list
if ! py_16s.py GenerateDatasetsIDsFile --FilePath "$METADATA" --Bioproject "$COL_BIOPROJECT"; then
    echo "❌ ERROR: Failed to generate dataset IDs"
    exit 1
fi

mapfile -t Dataset_ID_sets < <(awk '{print $1}' "${OUTPUT}/datasets_ID.txt")

if [ ${#Dataset_ID_sets[@]} -eq 0 ]; then
    echo "❌ ERROR: No datasets found"
    exit 1
fi

# Step 2: Generate SRA file list
if ! py_16s.py GenerateSRAsFile --FilePath "$METADATA" --Bioproject "$COL_BIOPROJECT" --SRA_Number "$COL_SRA"; then
    echo "❌ ERROR: Failed to generate SRA file lists"
    exit 1
fi

################################################################################
#                   PHASE 2: INDIVIDUAL DATASET PROCESSING                     #
################################################################################

echo "========================================="
echo "PHASE 2: Individual Dataset Processing"
echo "Started: $(date)"
echo "========================================="

failed_log="${OUTPUT}/failed_datasets.log"
success_log="${OUTPUT}/success_datasets.log"
skipped_log="${OUTPUT}/skipped_datasets.log"

: > "$failed_log"
: > "$success_log"
: > "$skipped_log"

for i in "${!Dataset_ID_sets[@]}"; do
    dataset_ID="${Dataset_ID_sets[$i]}"
    dataset_path="${OUTPUT}/${dataset_ID}/"
    sra_file_name="${dataset_ID}_sra.txt"
    platform="Unknown" 

    echo "----------------------------------------"
    echo "Dataset $((i+1))/${#Dataset_ID_sets[@]}: $dataset_ID"
    
    # Check if already processed
    if [ -f "${dataset_path}${dataset_ID}-final-rep-seqs.qza" ]; then
        echo "✓ Already processed. Skipping."
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - ALREADY_DONE" >> "$skipped_log"
        continue
    fi

    {
        cd "$dataset_path" || exit 1
        
        # 1. Download data
        echo ">>> Downloading SRA data..."
        Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}"

        # 2. Dynamic Platform Detection
        echo ">>> Detecting sequencing platform..."
        first_srr=$(awk 'NR==1 {print $1}' "${sra_file_name}")
        platform=$(py_16s.py get_sequencing_platform --srr_id "$first_srr")
        echo "Detected platform: $platform"

        # 3. Analyze sequence characteristics
        ori_fastq_path="${dataset_path}ori_fastq/"
        line_count=$(wc -l < "${dataset_path}${sra_file_name}")
        file_count=$(find "$ori_fastq_path" -type f 2>/dev/null | wc -l)
        
        if [ $((line_count * 2)) -eq $file_count ]; then
            sequence_type="paired"
        else
            sequence_type="single"
        fi
        echo "Sequence type: ${sequence_type^^}"

        # 4. Primer detection & Trimming
        echo ">>> Handling primers..."
        working_fastq_path="${dataset_path}working_fastq/"
        mkdir -p "$working_fastq_path"

        # Call primer detection and capture output
        primer_result=$(python "${SCRIPTS}/py_16s.py" detect_primers_16s \
            --input_path "$ori_fastq_path" \
            --tmp_path "${dataset_path}temp/" \
            --ref_path "${SCRIPT_DIR}/docs/J01859.1.fna" 2>&1)

        # Display the full output for debugging
        echo "$primer_result"
        echo ""

        # Parse result: output format is "TRIM:number"
        if [[ "$primer_result" =~ TRIM:([0-9]+) ]]; then
            trim_length="${BASH_REMATCH[1]}"

            if [[ "$trim_length" -gt 0 ]]; then
                echo "✓ Primers detected, trimming ${trim_length}bp from both ends..."
                if [[ "$sequence_type" == "single" ]]; then
                    for f in "${ori_fastq_path}"*.fastq*; do
                        echo "  Processing: $(basename "$f")"
                        fastp -i "$f" -o "${working_fastq_path}$(basename "$f")" \
                            -f "$trim_length" -w "$THREADS" -j /dev/null -h /dev/null || {
                            echo "  ✗ fastp failed for $(basename "$f")"
                            exit 1
                        }
                    done
                else
                    for r1 in "${ori_fastq_path}"*_R1*.fastq*; do
                        r2="${r1/_R1/_R2}"
                        echo "  Processing: $(basename "$r1") and $(basename "$r2")"
                        fastp -i "$r1" -I "$r2" \
                            -o "${working_fastq_path}$(basename "$r1")" \
                            -O "${working_fastq_path}$(basename "$r2")" \
                            -f "$trim_length" -F "$trim_length" -w "$THREADS" \
                            -j /dev/null -h /dev/null || {
                            echo "  ✗ fastp failed for $(basename "$r1")/$(basename "$r2")"
                            exit 1
                        }
                    done
                fi
                echo "✓ Fastp trimming completed"
            else
                echo "✓ Primers already removed, copying files..."
                cp "${ori_fastq_path}"*.fastq* "$working_fastq_path" 2>/dev/null || true
            fi
        else
            echo "✓ No primers detected or ambiguous, copying files..."
            cp "${ori_fastq_path}"*.fastq* "$working_fastq_path" 2>/dev/null || true
        fi

        # 5. Pipeline execution
        fastq_path="$working_fastq_path"
        export dataset_path sequence_type fastq_path
        export dataset_name="$dataset_ID"

        if [[ "$platform" == "OXFORD_NANOPORE" ]]; then
            echo "⚠️ Nanopore not supported yet."
            echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - SKIPPED - Platform: $platform" >> "$skipped_log"
            continue
        elif [[ "$platform" == "ILLUMINA" || "$platform" == "ION_TORRENT" ]]; then
            Amplicon_Common_MakeManifestFileForQiime2
            Amplicon_Common_ImportFastqToQiime2
            Amplicon_Illumina_QualityControl
            Amplicon_Illumina_DenosingDada2
            Amplicon_Common_FinalFilesCleaning
        elif [[ "$platform" == "PACBIO_SMRT" ]]; then
            Amplicon_Common_MakeManifestFileForQiime2
            Amplicon_Common_ImportFastqToQiime2
            Amplicon_Pacbio_QualityControlForQZA
            Amplicon_Pacbio_DenosingDada2
            Amplicon_Common_FinalFilesCleaning
        else
            echo "❌ Unknown platform: $platform"
            exit 1
        fi

        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - SUCCESS - Platform: $platform" >> "$success_log"

    } || {
        echo "❌ Pipeline failed for $dataset_ID"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - FAILED - Platform: $platform" >> "$failed_log"
    }
done

echo "========================================="
echo "ALL DONE. Check logs for summary."
echo "========================================="