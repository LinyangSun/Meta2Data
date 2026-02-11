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

# Source function library (AmpliconFunction.sh contains all processing functions)
if [[ -f "${SCRIPTS}/AmpliconFunction.sh" ]]; then
    source "${SCRIPTS}/AmpliconFunction.sh"
else
    echo "Error: AmpliconFunction.sh not found at ${SCRIPTS}/AmpliconFunction.sh"
    exit 1
fi

cd "$OUTPUT" || exit 1
mkdir -p "${OUTPUT}/analysis/"

################################################################################
#                          FIND REFERENCE FILE                                 #
################################################################################

# Find the 16S reference file in multiple possible locations
find_reference_file() {
    local ref_file="J01859.1.fna"
    local search_paths=(
        "${SCRIPT_DIR}/docs/${ref_file}"                    # Development/git repo
        "/usr/local/share/Meta2Data/docs/${ref_file}"       # System install
        "${HOME}/.local/share/Meta2Data/docs/${ref_file}"   # User install
        "${CONDA_PREFIX}/share/Meta2Data/docs/${ref_file}"  # Conda install
        "${PREFIX}/share/Meta2Data/docs/${ref_file}"        # Alternative prefix
    )

    for path in "${search_paths[@]}"; do
        if [[ -f "$path" ]]; then
            echo "$path"
            return 0
        fi
    done

    echo "ERROR: 16S reference file '${ref_file}' not found in any location:" >&2
    for path in "${search_paths[@]}"; do
        echo "  - $path" >&2
    done
    return 1
}

# Find and export reference file path
REF_16S_PATH=$(find_reference_file)
if [[ $? -ne 0 ]]; then
    echo "❌ ERROR: Cannot find 16S reference sequence file"
    exit 1
fi
echo "Using 16S reference: $REF_16S_PATH"

################################################################################
#                          PHASE 1: DATASET PREPARATION                        #
################################################################################

echo "========================================="
echo "PHASE 1: Dataset Preparation"
echo "Started: $(date)"
echo "========================================="

# Step 1: Generate dataset ID list
if ! python "${SCRIPTS}/py_16s.py" GenerateDatasetsIDsFile --FilePath "$METADATA" --Bioproject "$COL_BIOPROJECT" --OutputDir "$OUTPUT"; then
    echo "❌ ERROR: Failed to generate dataset IDs"
    exit 1
fi

mapfile -t Dataset_ID_sets < <(awk '{print $1}' "${OUTPUT}/datasets_ID.txt")

if [ ${#Dataset_ID_sets[@]} -eq 0 ]; then
    echo "❌ ERROR: No datasets found"
    exit 1
fi

# Step 2: Generate SRA file list
if ! python "${SCRIPTS}/py_16s.py" GenerateSRAsFile --FilePath "$METADATA" --Bioproject "$COL_BIOPROJECT" --SRA_Number "$COL_SRA" --OutputDir "$OUTPUT"; then
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

        # 1. Dynamic Platform Detection (BEFORE downloading)
        echo ">>> Detecting sequencing platform..."
        first_srr=$(awk 'NR==1 {print $1}' "${sra_file_name}")

        # Query API for platform (pass BioProject ID for CNCB/CRR accessions)
        if [[ "$first_srr" =~ ^CRR ]]; then
            # CRR accession - need to pass BioProject ID to CNCB API
            echo "  CNCB accession detected, using BioProject: $dataset_ID"
            platform=$(python "${SCRIPTS}/py_16s.py" get_sequencing_platform --srr_id "$first_srr" --bioproject_id "$dataset_ID")
        else
            # NCBI accession (SRR/ERR/DRR)
            platform=$(python "${SCRIPTS}/py_16s.py" get_sequencing_platform --srr_id "$first_srr")
        fi

        echo "Detected platform: $platform"

        # Skip unsupported platforms BEFORE downloading
        if [[ "$platform" == "OXFORD_NANOPORE" ]]; then
            echo "⚠️ Nanopore not supported yet. Skipping download and processing."
            echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - SKIPPED - Platform: $platform" >> "$skipped_log"
            continue
        fi

        # 2. Download data (only for supported platforms)
        echo ">>> Downloading SRA data..."
        Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}"

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

        # 4. Platform-specific processing
        export dataset_path sequence_type
        export dataset_name="$dataset_ID"

        if [[ "$platform" == "ILLUMINA" || "$platform" == "ION_TORRENT" ]]; then
            # Illumina/IonTorrent: Entropy-based primer detection & trimming
            echo ">>> Handling adapters and primers (entropy method)..."
            fastp_path="${dataset_path}temp/step_02_fastp/"
            mkdir -p "$fastp_path"

            python3 "${SCRIPTS}/entropy_primer_detect.py" \
                -i "$ori_fastq_path" \
                -o "$fastp_path" || {
                echo "  ✗ Entropy primer detection failed"
                exit 1
            }

            echo "✓ Entropy primer detection/removal completed"

            # Delete original fastq files to save space
            echo ">>> Cleaning up original files..."
            rm -rf "$ori_fastq_path"
            echo "✓ Removed ori_fastq/ to save space"

            # Run Illumina pipeline
            fastq_path="$fastp_path"
            export fastq_path
            Amplicon_Common_MakeManifestFileForQiime2
            Amplicon_Common_ImportFastqToQiime2
            Amplicon_Illumina_QualityControlForQZA
            Amplicon_Illumina_DenosingDada2
            Amplicon_Common_FinalFilesCleaning

        elif [[ "$platform" == "PACBIO_SMRT" ]]; then
            # PacBio: Entropy-based primer detection & trimming
            echo ">>> Handling adapters and primers (entropy method)..."
            fastp_path="${dataset_path}temp/step_02_fastp/"
            mkdir -p "$fastp_path"

            python3 "${SCRIPTS}/entropy_primer_detect.py" \
                -i "$ori_fastq_path" \
                -o "$fastp_path" || {
                echo "  ✗ Entropy primer detection failed"
                exit 1
            }

            echo "✓ Entropy primer detection/removal completed"

            # Delete original fastq files to save space
            echo ">>> Cleaning up original files..."
            rm -rf "$ori_fastq_path"
            echo "✓ Removed ori_fastq/ to save space"

            # Run PacBio pipeline
            fastq_path="$fastp_path"
            export fastq_path
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