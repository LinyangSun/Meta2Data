#!/bin/bash
set -e

################################################################################
#                          MULTI-VSEARCH PIPELINE                              #
################################################################################
#
# Purpose: Complete multi-platform amplicon data processing pipeline
# Phases: 1) Dataset preparation, 2) Smart trimming + 454 processing, 3) Merging
# Called by: Meta2Data-AmpliconPIP with --multi-vsearch flag
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
    -m, --metadata FILE       Input metadata CSV file
    -o, --output DIR          Output base directory

Optional options:
    -t, --threads INT         Number of CPU threads (default: 4)
    --skip-download BOOL      Skip SRA download (default: false)
    --keep-tmp BOOL           Keep temporary files (default: false)
    --col-bioproject NAME     Column for BioProject/Dataset ID (default: 'Data-Bioproject')
    --col-platform NAME       Column for platform (default: 'Data-SequencingPlatform')
    --col-sra NAME            Column for SRA accession (default: 'Data-SRA')
    --col-biosample NAME      Column for BioSample (default: 'Data-Biosample')
    --col-forward NAME        Column for forward primer (default: 'Data-Forward')
    --col-reverse NAME        Column for reverse primer (default: 'Data-Reverse')
    --col-region NAME         Column for target region (default: 'Data-Region')
    --db-dir DIR              GreenGenes2 database directory (auto-detected if not specified)
    -h, --help                Show this help message

Pipeline workflow:
    Phase 1: Dataset preparation
      - Generate dataset ID list
      - Generate SRA file lists

    Phase 2: Individual dataset processing
      - Download SRA data
      - Smart primer detection and trimming (auto-detects V1-V9 regions)
      - Platform-based pipeline selection:
        * ILLUMINA / ION_TORRENT -> 454 short-read pipeline
        * PACBIO_SMRT -> PacBio long-read pipeline
        * OXFORD_NANOPORE -> TBD (skipped)

    Phase 3: Multi-platform merging
      - Merge feature tables and sequences
      - GreenGenes2 annotation
      - Taxonomy assignment
      - Filter phylogenetic tree

Example:
    bash multi-vsearch.sh -m metadata.csv -o output/ -t 16

EOF
}

################################################################################
#                          DEFAULT PARAMETERS                                  #
################################################################################

METADATA=""
OUTPUT=""
THREADS=4
SKIP_DOWNLOAD=false
KEEP_TMP=false
COL_BIOPROJECT="Data-Bioproject"
COL_PLATFORM="Data-SequencingPlatform"
COL_SRA="Data-SRA"
COL_BIOSAMPLE="Data-Biosample"
COL_FORWARD="Data-Forward"
COL_REVERSE="Data-Reverse"
COL_REGION="Data-Region"
DB_DIR=""

################################################################################
#                          ARGUMENT PARSING                                    #
################################################################################

while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--metadata)
            METADATA="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --skip-download)
            SKIP_DOWNLOAD="$2"
            shift 2
            ;;
        --keep-tmp)
            KEEP_TMP="$2"
            shift 2
            ;;
        --col-bioproject)
            COL_BIOPROJECT="$2"
            shift 2
            ;;
        --col-platform)
            COL_PLATFORM="$2"
            shift 2
            ;;
        --col-sra)
            COL_SRA="$2"
            shift 2
            ;;
        --col-biosample)
            COL_BIOSAMPLE="$2"
            shift 2
            ;;
        --col-forward)
            COL_FORWARD="$2"
            shift 2
            ;;
        --col-reverse)
            COL_REVERSE="$2"
            shift 2
            ;;
        --col-region)
            COL_REGION="$2"
            shift 2
            ;;
        --db-dir)
            DB_DIR="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Error: Unknown option '$1'"
            show_help
            exit 1
            ;;
    esac
done

# Validate metadata/output input
if [[ -z "$METADATA" ]]; then
    echo "Error: --metadata is required"
    show_help
    exit 1
fi

# Check if metadata is a full/absolute path
if [[ "$METADATA" != /* ]]; then
    echo "Error: --metadata must be a full absolute path (starting with '/')"
    echo "       Provided: $METADATA"
    echo "       Example: /home/user/data/metadata.csv"
    exit 1
fi

# Check if metadata file exists
if [[ ! -f "$METADATA" ]]; then
    echo "Error: Metadata file '$METADATA' does not exist"
    exit 1
fi

# If OUTPUT is not specified, use METADATA's directory
if [[ -z "$OUTPUT" ]]; then
    OUTPUT=$(dirname "$METADATA")
    echo "Note: --output not specified, using metadata directory: $OUTPUT"
fi

# Export environment variables
export cpu=$THREADS

# Source function library
if [[ -f "${SCRIPTS}/Function_Import.sh" ]]; then
    source "${SCRIPTS}/Function_Import.sh"
else
    echo "Error: Function_Import.sh not found at ${SCRIPTS}/Function_Import.sh"
    exit 1
fi

cd "$OUTPUT" || { echo "Error: Cannot access $OUTPUT, it must be an existing directory"; exit 1; }
Analysis_path="${OUTPUT}/analysis/"
mkdir -p "${OUTPUT}/analysis/"

################################################################################
#                       PHASE 1: DATASET PREPARATION                           #
################################################################################

echo "========================================="
echo "PHASE 1: Dataset Preparation"
echo "Started: $(date)"
echo "========================================="
echo ""

# Step 1: Generate dataset ID list
echo ">>> Step 1: Generating dataset IDs..."
echo "Using column names:"
echo "  BioProject: $COL_BIOPROJECT"
echo "  Platform:   $COL_PLATFORM"

if ! py_16s.py GenerateDatasetsIDsFile \
    --FilePath "$METADATA" \
    --Bioproject "$COL_BIOPROJECT"; then
    echo "❌ ERROR: Failed to generate dataset IDs based on BioProject"
    exit 1
fi

# Validate datasets_ID.txt exists
if [ ! -f "${OUTPUT}/datasets_ID.txt" ]; then
    echo "❌ ERROR: datasets_ID.txt was not created, it was expected at ${OUTPUT}/datasets_ID.txt"
    exit 1
fi

# Load dataset arrays
mapfile -t Dataset_ID_sets < <(awk '{print $1}' "${OUTPUT}/datasets_ID.txt")

if [ ${#Dataset_ID_sets[@]} -eq 0 ]; then
    echo "❌ ERROR: No datasets found"
    exit 1
fi

echo "✓ Found ${#Dataset_ID_sets[@]} datasets"
echo ""

# Step 2: Generate SRA file list
echo ">>> Step 2: Generating SRA file lists..."
echo "Using column names:"
echo "  BioProject: $COL_BIOPROJECT"
echo "  SRA:        $COL_SRA"
echo "  BioSample:  $COL_BIOSAMPLE"
echo "  Forward:    $COL_FORWARD"
echo "  Reverse:    $COL_REVERSE"
echo "  Region:     $COL_REGION"

if ! py_16s.py GenerateSRAsFile \
    --FilePath "$METADATA" \
    --Bioproject "$COL_BIOPROJECT" \
    --SRA_Number "$COL_SRA"; then
    echo "❌ ERROR: Failed to generate SRA file lists"
    exit 1
fi

echo "✓ SRA file lists generated"
echo ""

################################################################################
#                  PHASE 2: INDIVIDUAL DATASET PROCESSING                      #
################################################################################

echo "========================================="
echo "PHASE 2: Individual Dataset Processing"
echo "  (Smart trimming + Platform-based pipeline)"
echo "Started: $(date)"
echo "========================================="
echo ""

# Initialize log files
failed_log="${OUTPUT}/failed_datasets.log"
success_log="${OUTPUT}/success_datasets.log"
skipped_log="${OUTPUT}/skipped_datasets.log"

: > "$failed_log"
: > "$success_log"
: > "$skipped_log"

echo ""

# Process each dataset
for i in "${!Dataset_ID_sets[@]}"; do

    dataset_ID="${Dataset_ID_sets[$i]}"
    original_platform="${SequencingPlatform_sets[$i]}"
    dataset_path="${OUTPUT}/${dataset_ID}/"
    sra_file_name="${dataset_ID}_sra.txt"

    echo "========================================"
    echo "Dataset $((i+1))/${#Dataset_ID_sets[@]}: $dataset_ID"
    echo "Original platform: $original_platform"
    echo "Path: $dataset_path"
    echo "========================================"

    # Check if already processed
    if [ -f "${dataset_path}${dataset_ID}-final-rep-seqs.qza" ]; then
        echo "✓ Already processed. Skipping."
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - ALREADY_DONE" >> "$skipped_log"
        echo ""
        continue
    fi

    # Clean previous incomplete outputs
    echo "Cleaning previous outputs..."
    find "$dataset_path" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} + 2>/dev/null || true

    # Process dataset
    {
        cd "$dataset_path" || { echo "❌ ERROR: Cannot access $dataset_path"; exit 1; }
        
        ##########################################################################
        # Step 1: Download data
        if [[ "$SKIP_DOWNLOAD" == false ]]; then
            echo ""
            echo ">>> Downloading SRA data..."
            if ! Common_SRADownloadToFastq_Common_SraToolkits -d "$dataset_path" -a "${sra_file_name}"; then
                echo "❌ ERROR: SRA download failed"
                exit 1
            fi
            echo "✓ Download completed"
        fi

        ##########################################################################
        # Detect whether paired-end or single-end
        echo ""
        echo ">>> Analyzing sequence characteristics..."
        sra_file_path="${dataset_path}${sra_file_name}"
        ori_fastq_path="${dataset_path}ori_fastq/"

        line_count=$(wc -l < "$sra_file_path")
        file_count=$(find "$ori_fastq_path" -type f 2>/dev/null | wc -l)
        
        # Determine if paired-end or single-end
        if [ $((line_count * 2)) -eq $file_count ]; then
            echo "Sequence type: PAIRED-END"
            sequence_type="paired"
        elif [ $line_count -eq $file_count ]; then
            echo "Sequence type: SINGLE-END"
            sequence_type="single"
        else
            echo "❌ ERROR: Mismatch between SRA lines ($line_count) and FASTQ files ($file_count)"
            exit 1
        fi
        ##########################################################################

        # Detect sequencing platform
        first_element=$(awk 'NR==1 {print $1}' "${sra_file_name}")
        platform=$(py_16s.py get_sequencing_platform \
            --srr_id "$first_element")
        
        echo "Detected platform: $platform"

        ##########################################################################
        # Step 3: Detect primers
        echo ""
        echo ">>> Detecting primers..."
        if py_16s.py detect_primers_16s \
                --input_path "$ori_fastq_path" \
                --tmp_path "${dataset_path}temp/" \
                --ref_path "${SCRIPT_DIR}/docs/J01859.1.fna" 2>&1; then
            PRIMERS_DETECTED="yes"
            echo "✓ Primers detected"
        else
            PRIMERS_DETECTED="no"
            echo "✓ No primers detected"
        fi

        ##########################################################################
        # Step 4: Create working_fastq directory
        # This ensures a consistent input path for downstream pipeline
        ##########################################################################
        echo ""
        echo ">>> Preparing working fastq directory..."
        
        working_fastq_path="${dataset_path}working_fastq/"
        mkdir -p "$working_fastq_path"

        if [[ "$PRIMERS_DETECTED" == "yes" ]]; then
            # Primers detected: trim and output to working_fastq
            echo ">>> Trimming primers with fastp (20bp from 5' end)..."

            if [[ "$sequence_type" == "single" ]]; then
                # Single-end trimming with fastp
                echo "  Processing single-end data..."
                processed_count=0
                
                for fastq_file in "${ori_fastq_path}"*.fastq*; do
                    if [ ! -f "$fastq_file" ]; then
                        continue
                    fi

                    base_name=$(basename "$fastq_file")
                    out_file="${working_fastq_path}${base_name}"

                    # fastp: trim 20bp from front (-f 20)
                    if fastp -i "$fastq_file" -o "$out_file" \
                        -f 20 \
                        -w "$THREADS" \
                        -j "${working_fastq_path}${base_name}.json" \
                        -h "${working_fastq_path}${base_name}.html" 2>/dev/null; then
                        ((processed_count++))
                    else
                        echo "  ⚠️  Trimming failed for $base_name"
                    fi
                done

                echo "✓ Trimmed $processed_count single-end files to working_fastq/"

            elif [[ "$sequence_type" == "paired" ]]; then
                # Paired-end trimming with fastp
                echo "  Processing paired-end data..."
                processed_count=0

                for r1 in "${ori_fastq_path}"*_R1*.fastq*; do
                    if [ ! -f "$r1" ]; then
                        continue
                    fi

                    r2="${r1/_R1/_R2}"
                    if [ ! -f "$r2" ]; then
                        echo "  ⚠️  No R2 for $(basename "$r1"), skipping"
                        continue
                    fi

                    base_name=$(basename "$r1")
                    out1="${working_fastq_path}$(basename "$r1")"
                    out2="${working_fastq_path}$(basename "$r2")"

                    # fastp: trim 20bp from front of both reads (-f 20 -F 20)
                    if fastp -i "$r1" -I "$r2" -o "$out1" -O "$out2" \
                        -f 20 -F 20 \
                        -w "$THREADS" \
                        -j "${working_fastq_path}${base_name}.json" \
                        -h "${working_fastq_path}${base_name}.html" 2>/dev/null; then
                        ((processed_count++))
                    else
                        echo "  ⚠️  Trimming failed for $base_name"
                    fi
                done

                echo "✓ Trimmed $processed_count paired-end files to working_fastq/"
            fi

        else
            # No primers detected: copy original files to working_fastq
            echo ">>> No primer trimming needed - copying original files to working_fastq/..."
            
            cp "${ori_fastq_path}"*.fastq* "$working_fastq_path" 2>/dev/null || true
            
            copied_count=$(find "$working_fastq_path" -type f -name "*.fastq*" 2>/dev/null | wc -l)
            echo "✓ Copied $copied_count files to working_fastq/"
        fi

        # Verify working_fastq has files
        working_file_count=$(find "$working_fastq_path" -type f -name "*.fastq*" 2>/dev/null | wc -l)
        if [ "$working_file_count" -eq 0 ]; then
            echo "❌ ERROR: No fastq files in working_fastq directory"
            exit 1
        fi
        echo "✓ Working directory ready with $working_file_count fastq files"

        ##########################################################################
        # Update fastq_path to use working_fastq for all downstream steps
        ##########################################################################
        fastq_path="$working_fastq_path"

        # Export variables for pipeline
        export dataset_path F_primer R_primer REVRC sequence_type
        export dataset_name="$dataset_ID"

        ##########################################################################
        # Step 5: Run appropriate pipeline based on platform
        ##########################################################################
        echo ""
        echo ">>> Selecting pipeline based on platform..."
        echo "  Platform: ${platform}"

        # Skip OXFORD_NANOPORE (TBD)
        if [[ "$platform" == "OXFORD_NANOPORE" ]]; then
            echo "  Status: TBD"
            echo ""
            echo "⚠️  Oxford Nanopore pipeline not yet implemented"
            echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - SKIPPED - Platform: OXFORD_NANOPORE (TBD)" >> "$skipped_log"
            continue
        fi

        # ILLUMINA or ION_TORRENT -> 454 short-read pipeline
        if [[ "$platform" == "ILLUMINA" || "$platform" == "ION_TORRENT" ]]; then
            echo "  Selected pipeline: 454 (short-read)"
            echo ""
            echo ">>> Running 454 pipeline..."

            Amplicon_Common_MakeManifestFileForQiime2 || { echo "❌ Failed at MakeManifestFileForQiime2"; exit 1; }
            Amplicon_Common_ImportFastqToQiime2 || { echo "❌ Failed at ImportFastqToQiime2"; exit 1; }
            Amplicon_Illumina_QualityControl || { echo "❌ Failed at QualityControl"; exit 1; }
            Amplicon_Illumina_DenosingDada2 || { echo "❌ Failed at DenosingDada"; exit 1; }
             Amplicon_Common_FinalFilesCleaning || { echo "Failed at FinalFilesCleaning"; exit 1; }
        fi

        # PACBIO_SMRT -> PacBio long-read pipeline
        if [[ "$platform" == "PACBIO_SMRT" ]]; then
            echo "  Selected pipeline: PacBio (long-read)"
            echo ""
            echo ">>> Running PacBio pipeline..."

            Amplicon_Common_MakeManifestFileForQiime2 || { echo "❌ Failed at MakeManifestFileForQiime2"; exit 1; }
            Amplicon_Common_ImportFastqToQiime2 || { echo "❌ Failed at ImportFastqToQiime2"; exit 1; }
            Amplicon_Pacbio_QualityControlForQZA || { echo "❌ Failed at QualityControlForQZA"; exit 1; }
            Amplicon_Pacbio_DenosingDada2 || { echo "❌ Failed at DenosingDada2"; exit 1; }
            Amplicon_Common_FinalFilesCleaning || { echo "Failed at FinalFilesCleaning"; exit 1; }
        fi

        # Cleanup temporary files
        if [[ "$KEEP_TMP" == false ]]; then
            echo ""
            echo ">>> Cleaning temporary files..."
            rm -rf "${dataset_path}ori_fastq/" 2>/dev/null || true
            rm -rf "${dataset_path}working_fastq/" 2>/dev/null || true
            rm -rf "${dataset_path}temp"* 2>/dev/null || true
        fi

        echo ""
        echo "✓ Pipeline completed for $dataset_ID"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - SUCCESS - Platform: $platform" >> "$success_log"

    } || {
        echo ""
        echo "❌ ERROR: Pipeline failed for $dataset_ID"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - FAILED - Platform: $original_platform" >> "$failed_log"
    }

    echo ""

done



# Phase 2 summary
echo "================================"
echo "PHASE 2 COMPLETE"
echo "================================"
echo "Total datasets: ${#Dataset_ID_sets[@]}"
echo "Successful: $(wc -l < "$success_log" 2>/dev/null || echo 0)"
echo "Failed: $(wc -l < "$failed_log" 2>/dev/null || echo 0)"
echo "Skipped: $(wc -l < "$skipped_log" 2>/dev/null || echo 0)"
echo ""

################################################################################
#                          DATABASE RESOLUTION                                 #
################################################################################

echo "========================================="
echo "Resolving GreenGenes2 database..."
echo "========================================="

# Database auto-detection function
find_database_dir() {
    if [[ -n "$DB_DIR" && -d "$DB_DIR" ]]; then
        echo "$DB_DIR"
        return 0
    fi

    local COMMON_PATHS=(
        "$CONDA_PREFIX/share/qiime2/data/greengenes2"
        "$HOME/.qiime2/db/greengenes2"
        "/scratch/project_2009135/db/gg2"
        "/usr/local/share/qiime2/data/greengenes2"
    )

    for path in "${COMMON_PATHS[@]}"; do
        if [[ -d "$path" ]] && [[ -n "$(ls -A "$path" 2>/dev/null)" ]]; then
            echo "$path"
            return 0
        fi
    done

    echo ""
    return 1
}

# Find specific database file
find_db_file() {
    local pattern="$1"
    local db_dir="$2"
    local files=("${db_dir}/"*${pattern}*.qza)

    if [[ -f "${files[0]}" ]]; then
        echo "${files[0]}"
        return 0
    fi

    echo ""
    return 1
}

# Resolve database
DB_DIR=$(find_database_dir)

if [[ -z "$DB_DIR" ]]; then
    echo "❌ ERROR: GreenGenes2 database not found"
    echo ""
    echo "Searched locations:"
    echo "  - \$CONDA_PREFIX/share/qiime2/data/greengenes2"
    echo "  - \$HOME/.qiime2/db/greengenes2"
    echo "  - /scratch/project_2009135/db/gg2"
    echo ""
    echo "Install with: conda install -c conda-forge -c bioconda q2-greengenes2"
    echo "Or specify: --db-dir /path/to/greengenes2"
    exit 1
fi

echo "✓ Database found: $DB_DIR"

# Resolve individual files
BACKBONE=$(find_db_file "backbone" "$DB_DIR")
TAXONOMY=$(find_db_file "taxonomy" "$DB_DIR")
PHYLO=$(find_db_file "phylogeny" "$DB_DIR")

if [[ -z "$BACKBONE" ]] || [[ -z "$TAXONOMY" ]] || [[ -z "$PHYLO" ]]; then
    echo "❌ ERROR: Required database files not found"
    exit 1
fi

echo "✓ Backbone: $(basename "$BACKBONE")"
echo "✓ Taxonomy: $(basename "$TAXONOMY")"
echo "✓ Phylogeny: $(basename "$PHYLO")"
echo ""

################################################################################
#                    PHASE 3: MERGE AND TAXONOMY ASSIGNMENT                    #
################################################################################

echo "========================================="
echo "PHASE 3: Merge & Taxonomy Assignment"
echo "Started: $(date)"
echo "========================================="
echo ""

# Define output directories
FINAL_DIR="${OUTPUT}/final"
TMP_DIR="${FINAL_DIR}/tmp"
MERGED_DIR="${FINAL_DIR}/merged"

# Create directories
echo ">>> Step 1: Preparing output directories..."
mkdir -p "${FINAL_DIR}" "${TMP_DIR}" "${MERGED_DIR}"
echo "✓ Directories ready"
echo ""

# Collect dataset folders
echo ">>> Step 2: Collecting dataset outputs..."

all_folders=()
for folder in "${OUTPUT}"/PRJ*/; do
    if [ -d "$folder" ]; then
        all_folders+=("$folder")
    fi
done

if [ ${#all_folders[@]} -eq 0 ]; then
    echo "❌ ERROR: No dataset folders found"
    exit 1
fi

echo "Found ${#all_folders[@]} dataset folder(s)"
echo ""

# Collect tables and sequences
echo ">>> Step 3: Collecting QZA files..."

ALL_TABLES=()
ALL_REP_SEQS=()
SUCCESSFUL_DATASETS=()

for folder in "${all_folders[@]}"; do
    folder_name=$(basename "$folder")
    echo "Processing: $folder_name"

    rep_files=("${folder}"/*-final-rep-seqs.qza)
    table_files=("${folder}"/*-final-table.qza)

    if [ ${#rep_files[@]} -eq 0 ] || [ ! -f "${rep_files[0]}" ]; then
        echo "  ⚠️  No rep-seqs file, skipping..."
        continue
    fi

    if [ ${#table_files[@]} -eq 0 ] || [ ! -f "${table_files[0]}" ]; then
        echo "  ⚠️  No table file, skipping..."
        continue
    fi

    ALL_TABLES+=("${table_files[0]}")
    ALL_REP_SEQS+=("${rep_files[0]}")
    SUCCESSFUL_DATASETS+=("$folder_name")

    echo "  ✓ Collected"
done

if [ ${#ALL_TABLES[@]} -eq 0 ]; then
    echo "❌ ERROR: No valid datasets collected"
    exit 1
fi

echo ""
echo "Successfully collected ${#ALL_TABLES[@]} dataset(s)"
echo ""

# Merge feature tables
echo ">>> Step 4: Merging feature tables..."
MERGED_TABLE="${MERGED_DIR}/merged-table.qza"

if ! qiime feature-table merge \
    --i-tables "${ALL_TABLES[@]}" \
    --o-merged-table "$MERGED_TABLE" --verbose; then
    echo "❌ ERROR: Table merging failed"
    exit 1
fi

echo "✓ Tables merged"
echo ""

# Merge sequences
echo ">>> Step 5: Merging representative sequences..."
MERGED_REP_SEQS="${MERGED_DIR}/merged-rep-seqs.qza"

if ! qiime feature-table merge-seqs \
    --i-data "${ALL_REP_SEQS[@]}" \
    --o-merged-data "$MERGED_REP_SEQS" --verbose; then
    echo "❌ ERROR: Sequence merging failed"
    exit 1
fi

echo "✓ Sequences merged"
echo ""

# Generate summary
echo ">>> Step 6: Generating merged table summary..."
if qiime feature-table summarize \
    --i-table "$MERGED_TABLE" \
    --o-visualization "${MERGED_DIR}/merged-table-summary.qzv" --verbose; then
    echo "✓ Summary generated"
else
    echo "⚠️  WARNING: Summary generation failed"
fi
echo ""

# GreenGenes2 annotation
echo ">>> Step 7: Running GreenGenes2 annotation..."
echo "Using $cpu CPU threads"

ANNOTATED_TABLE="${MERGED_DIR}/merged-table-gg2.qza"
ANNOTATED_REP_SEQS="${MERGED_DIR}/merged-rep-seqs-gg2.qza"

if ! qiime greengenes2 non-v4-16s \
    --i-table "$MERGED_TABLE" \
    --i-sequences "$MERGED_REP_SEQS" \
    --p-threads "$cpu" \
    --i-backbone "${BACKBONE}" \
    --o-mapped-table "$ANNOTATED_TABLE" \
    --o-representatives "$ANNOTATED_REP_SEQS" --verbose; then
    echo "❌ ERROR: GreenGenes2 annotation failed"
    exit 1
fi

echo "✓ GreenGenes2 mapping completed"
echo ""

# Taxonomy assignment
echo ">>> Step 8: Assigning taxonomy..."
MERGED_TAXONOMY="${MERGED_DIR}/merged-taxonomy.qza"

if ! qiime greengenes2 taxonomy-from-table \
    --i-reference-taxonomy "${TAXONOMY}" \
    --i-table "$ANNOTATED_TABLE" \
    --o-classification "$MERGED_TAXONOMY" --verbose; then
    echo "❌ ERROR: Taxonomy assignment failed"
    exit 1
fi

echo "✓ Taxonomy assigned"
echo ""

# Generate annotated summary
echo ">>> Step 9: Generating annotated table summary..."
if qiime feature-table summarize \
    --i-table "$ANNOTATED_TABLE" \
    --o-visualization "${MERGED_DIR}/merged-table-gg2-summary.qzv" --verbose; then
    echo "✓ Annotated summary generated"
else
    echo "⚠️  WARNING: Annotated summary generation failed"
fi
echo ""

# Filter tree
echo ">>> Step 10: Filtering phylogenetic tree..."
FILTERED_TREE="${MERGED_DIR}/final-tree.qza"

if ! qiime phylogeny filter-tree \
    --i-tree "$PHYLO" \
    --i-table "$ANNOTATED_TABLE" \
    --o-filtered-tree "$FILTERED_TREE" --verbose; then
    echo "❌ ERROR: Tree filtering failed"
    exit 1
fi

echo "✓ Phylogenetic tree filtered"
echo ""

################################################################################
#                           PIPELINE COMPLETION                                #
################################################################################

echo "========================================="
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "Finished: $(date)"
echo "========================================="
echo ""
echo "Output files: ${MERGED_DIR}"
echo ""
echo "Final outputs:"
echo "  - Merged table:          merged-table.qza"
echo "  - Merged sequences:      merged-rep-seqs.qza"
echo "  - Annotated table:       merged-table-gg2.qza"
echo "  - Annotated sequences:   merged-rep-seqs-gg2.qza"
echo "  - Taxonomy:              merged-taxonomy.qza"
echo "  - Filtered tree:         final-tree.qza"
echo ""
echo "Processed datasets: ${#SUCCESSFUL_DATASETS[@]}"
for dataset in "${SUCCESSFUL_DATASETS[@]}"; do
    echo "  - $dataset"
done
echo ""
echo "Log files:"
echo "  - Success: $success_log"
echo "  - Failed:  $failed_log"
echo "  - Skipped: $skipped_log"
echo ""
echo "========================================="