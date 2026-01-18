#!/usr/bin/env bash
#SBATCH --job-name=MetaAnalysis
#SBATCH --output=MetaAnalysis-all-%j.out
#SBATCH --error=MetaAnalysis-all-%j.err
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --partition=small
#SBATCH --account=project_2009135
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=linyang.sun@helsinki.fi
#SBATCH --gres=nvme:100


################################################################################
#                          PARAMETER CONFIGURATION                             #
################################################################################

# Base directories
BASE_DIR="/scratch/project_2009135/linyang/meta-analysis"

# Input metadata file
metafile="${BASE_DIR}/combined_OriMetaCSV.csv"

# CPU configuration
cpu=$SLURM_CPUS_PER_TASK
FINAL_DIR="${BASE_DIR}/final1"
TMP_DIR="${FINAL_DIR}/tmp"
MERGED_DIR="${FINAL_DIR}/merged"

# Database directories
DB_DIR="/scratch/project_2009135/db/gg2"
BACKBONE="${DB_DIR}/2024.09.backbone.full-length.fna.qza"
TAXONOMY="${DB_DIR}/2024.09.taxonomy.asv.nwk.qza"
PHYLO="${DB_DIR}/2024.09.phylogeny.id.nwk.qza"


################################################################################
#                         ENVIRONMENT INITIALIZATION                           #
################################################################################

echo "========================================="
echo "Initializing environment..."
echo "========================================="

# Load required modules and environment
source /scratch/project_2009135/linyang/common_code/Function_Import.sh
export PATH="/projappl/project_2009135/linyang/qiime2_gg2/bin:$PATH"
module load sratoolkit/3.0.0

echo "Environment activated successfully"
echo "CPU cores allocated: $cpu"
echo ""


################################################################################
#                         VALIDATION AND SETUP                                 #
################################################################################

echo "========================================="
echo "Validating required files..."
echo "========================================="

# Validate required database files
for required_file in "$BACKBONE" "$TAXONOMY" "$PHYLO"; do
    if [ ! -f "$required_file" ]; then
        echo "❌ ERROR: Required database file not found: $required_file"
        exit 1
    fi
    echo "✓ Found: $(basename "$required_file")"
done

# Validate metadata file
if [ ! -f "$metafile" ]; then
    echo "❌ ERROR: Metadata file not found: $metafile"
    exit 1
fi
echo "✓ Found: $(basename "$metafile")"
echo ""

# Setup cleanup trap
cleanup() {
    if [ $? -ne 0 ]; then
        echo ""
        echo "========================================="
        echo "Pipeline interrupted or failed"
        echo "Check log files for details"
        echo "========================================="
    fi
}
trap cleanup EXIT


################################################################################
#                       DATASET PREPARATION PHASE                              #
################################################################################

echo "========================================="
echo "PHASE 1: Dataset Preparation"
echo "Started: $(date)"
echo "========================================="
echo ""

# Navigate to base directory
cd "$BASE_DIR" || { echo "❌ ERROR: Cannot access $BASE_DIR"; exit 1; }

# Step 1: Generate BioProject list (replacing dataset ID)
echo ">>> Step 1: Generating BioProject list..."
if ! py_16s.py GenerateDatasetsIDsFile \
    --FilePath "$metafile" \
    --Datasets_ID 'Data-Bioproject' \
    --SequencingPlatform "Data-SequencingPlatform"; then
    echo "❌ ERROR: Failed to generate BioProject list"
    exit 1
fi

# Validate datasets_ID.txt exists (still using this filename for compatibility)
if [ ! -f "${BASE_DIR}/datasets_ID.txt" ]; then
    echo "❌ ERROR: datasets_ID.txt was not created"
    exit 1
fi

# Read BioProject IDs and platforms
mapfile -t Dataset_ID_sets < <(awk '{print $1}' "${BASE_DIR}/datasets_ID.txt")
mapfile -t SequencingPlatform_sets < <(awk '{print $2}' "${BASE_DIR}/datasets_ID.txt")

# Validate array sizes match
if [ ${#Dataset_ID_sets[@]} -ne ${#SequencingPlatform_sets[@]} ]; then
    echo "❌ ERROR: Mismatch between BioProject IDs (${#Dataset_ID_sets[@]}) and platforms (${#SequencingPlatform_sets[@]})"
    exit 1
fi

if [ ${#Dataset_ID_sets[@]} -eq 0 ]; then
    echo "❌ ERROR: No BioProjects found in datasets_ID.txt"
    exit 1
fi

echo "Found ${#Dataset_ID_sets[@]} BioProjects to process"
echo ""

# Step 2: Generate SRA file list
echo ">>> Step 2: Generating SRA file list..."
if ! py_16s.py GenerateSRAsFile \
    --FilePath "$metafile" \
    --Datasets_ID 'Data-Bioproject' \
    --Bioproject 'Data-Bioproject' \
    --SRA_Number 'Data-SRA' \
    --Biosample 'Data-Biosample' \
    --Forward 'Data-Forward' \
    --Reverse 'Data-Reverse' \
    --region 'Data-Region'; then
    echo "❌ ERROR: Failed to generate SRA file list"
    exit 1
fi

echo "SRA file list generated successfully"
echo ""


################################################################################
#                       INDIVIDUAL DATASET PROCESSING                          #
################################################################################

echo "========================================="
echo "PHASE 2: Processing Individual Datasets"
echo "Started: $(date)"
echo "========================================="
echo ""

# Initialize log files
failed_log="${BASE_DIR}/failed_datasets.log"
success_log="${BASE_DIR}/success_datasets.log"
skipped_log="${BASE_DIR}/skipped_datasets.log"

: > "$failed_log"
: > "$success_log"
: > "$skipped_log"

# Process each dataset
for i in "${!Dataset_ID_sets[@]}"; do
    
    # -----------------------------------------------------------------------
    # Dataset initialization (using BioProject ID as identifier)
    # -----------------------------------------------------------------------
    dataset_ID="${Dataset_ID_sets[$i]}"  # BioProject ID (e.g., PRJNA123456)
    platform="${SequencingPlatform_sets[$i]}"
    platform="${platform,,}"  # Convert to lowercase
    dataset_path="${BASE_DIR}/${dataset_ID}/"  # Folder named by BioProject ID
    sra_file_name="${dataset_ID}_sra.txt"

    echo "================================"
    echo "BioProject $((i+1))/${#Dataset_ID_sets[@]}: $dataset_ID"
    echo "Platform: $platform"
    echo "Path: $dataset_path"
    echo "================================"
    
    # -----------------------------------------------------------------------
    # Check if already processed
    # -----------------------------------------------------------------------
    if [ -f "${dataset_path}${dataset_ID}-final-rep-seqs.qza" ]; then
        echo "✓ Dataset already processed. Skipping."
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - ALREADY_DONE - Platform: $platform" >> "$skipped_log"
        echo ""
        continue
    fi
    
    # -----------------------------------------------------------------------
    # Clean previous outputs
    # -----------------------------------------------------------------------
    echo "Cleaning previous outputs..."
    find "$dataset_path" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} + 2>/dev/null || true
    
    # -----------------------------------------------------------------------
    # Process dataset (with error handling)
    # -----------------------------------------------------------------------
    {
        cd "$dataset_path" || { echo "❌ ERROR: Cannot access $dataset_path"; exit 1; }
        
        # Step 1: Download data
        echo ""
        echo ">>> Downloading SRA data..."
        if ! Common_SRADownloadToFastq_Common_SraToolkits -d "$dataset_path" -a "${sra_file_name}"; then
            echo "❌ ERROR: SRA download failed"
            exit 1
        fi
        echo "✓ Download completed"
        
        # Step 2: Analyze sequence type and primers
        echo ""
        echo ">>> Analyzing sequence characteristics..."
        sra_file_path="${dataset_path}${sra_file_name}"
        fastq_path="${dataset_path}ori_fastq/"
        
        # Validate SRA file exists
        if [ ! -f "$sra_file_path" ]; then
            echo "❌ ERROR: SRA file not found: $sra_file_path"
            exit 1
        fi
        
        F_primer=$(cut -d$'\t' -f4 "$sra_file_path" | uniq)
        R_primer=$(cut -d$'\t' -f5 "$sra_file_path" | uniq)
        REVRC=$(echo "$R_primer" | tr ACGTacgt TGCAtgca | rev)
        
        line_count=$(wc -l < "$sra_file_path")
        file_count=$(find "$fastq_path" -type f 2>/dev/null | wc -l)
        
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
        
        echo "Forward primer: $F_primer"
        echo "Reverse primer: $R_primer"
        echo "Files found: $file_count"
        
        # Export variables for downstream processing
        export dataset_path F_primer R_primer REVRC sequence_type
        export dataset_name="$dataset_ID"
        
        # Step 3: Platform-specific processing pipeline
        echo ""
        echo ">>> Running $platform pipeline..."
        
        if [[ "$platform" == "454" ]]; then
            # ---- 454 SEQUENCING PIPELINE ----
            echo "Pipeline: 454 Pyrosequencing"
            Amplicon_454_ImportToQiime2 || { echo "Failed at ImportToQiime2"; exit 1; }
            Amplicon_454_QualityControl || { echo "Failed at QualityControl"; exit 1; }
            Amplicon_454_Deduplication || { echo "Failed at Deduplication"; exit 1; }
            Amplicon_454_ChimerasRemoval || { echo "Failed at ChimerasRemoval"; exit 1; }
            Amplicon_454_ClusterDenovo || { echo "Failed at ClusterDenovo"; exit 1; }
            Amplicon_Common_FinalFilesCleaning || { echo "Failed at FinalFilesCleaning"; exit 1; }
            
        elif [[ "$platform" == "pacbio" ]]; then
            # ---- PACBIO SEQUENCING PIPELINE ----
            echo "Pipeline: PacBio Long-Read"
            Amplicon_Pacbio_PrimerDetectionAndQualityControl || { echo "Failed at PrimerDetection"; exit 1; }
            Amplicon_Common_MakeManifestFileForQiime2 || { echo "Failed at MakeManifest"; exit 1; }
            Amplicon_Common_ImportFastqToQiime2 || { echo "Failed at ImportFastq"; exit 1; }
            Amplicon_Pacbio_QualityControlForQZA || { echo "Failed at QualityControl"; exit 1; }
            Amplicon_Pacbio_DenosingDada2 || { echo "Failed at Denoising"; exit 1; }
            Amplicon_Common_FinalFilesCleaning || { echo "Failed at FinalFilesCleaning"; exit 1; }
            
        else
            # ---- ILLUMINA SEQUENCING PIPELINE (DEFAULT) ----
            echo "Pipeline: Illumina Short-Read"
            Amplicon_Illumina_PrimerDetectionAndQualityControl || { echo "Failed at PrimerDetection"; exit 1; }
            Amplicon_Common_MakeManifestFileForQiime2 || { echo "Failed at MakeManifest"; exit 1; }
            Amplicon_Common_ImportFastqToQiime2 || { echo "Failed at ImportFastq"; exit 1; }
            Amplicon_Illumina_QualityControlForQZA || { echo "Failed at QualityControl"; exit 1; }
            Amplicon_Illumina_DenosingDeblur || { echo "Failed at Denoising"; exit 1; }
            Amplicon_Common_FinalFilesCleaning || { echo "Failed at FinalFilesCleaning"; exit 1; }
        fi
        
        echo ""
        echo "✅ Dataset $dataset_ID completed successfully"
        
    } && {
        # Log success with timestamp
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - SUCCESS - Platform: $platform" >> "$success_log"
        
    } || {
        # Log failure and continue
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - FAILED - Platform: $platform" >> "$failed_log"
        echo "❌ Dataset $dataset_ID failed. Continuing to next dataset..."
    }
    
    echo ""
    
done

# -----------------------------------------------------------------------
# Processing summary
# -----------------------------------------------------------------------
echo "================================"
echo "PHASE 2 COMPLETE"
echo "================================"
echo "Total datasets: ${#Dataset_ID_sets[@]}"
echo "Successful: $(wc -l < "$success_log" 2>/dev/null || echo 0)"
echo "Failed: $(wc -l < "$failed_log" 2>/dev/null || echo 0)"
echo "Skipped: $(wc -l < "$skipped_log" 2>/dev/null || echo 0)"
echo ""


################################################################################
#                    MERGE AND TAXONOMY ASSIGNMENT PHASE                       #
################################################################################

echo "========================================="
echo "PHASE 3: Merge & Taxonomy Assignment"
echo "Started: $(date)"
echo "========================================="
echo ""

# -----------------------------------------------------------------------
# Directory setup
# -----------------------------------------------------------------------
echo ">>> Step 1: Preparing output directories..."
mkdir -p "${FINAL_DIR}"
mkdir -p "${TMP_DIR}"
mkdir -p "${MERGED_DIR}"
echo "Directories ready"
echo ""

# -----------------------------------------------------------------------
# Collect all dataset folders
# -----------------------------------------------------------------------
echo ">>> Step 2: Identifying dataset folders..."

# Find all BioProject folders (PRJ* pattern)
# BioProject IDs: PRJNA*, PRJEB*, PRJCA*, PRJDB*
all_folders=()
for folder in "${BASE_DIR}"/PRJ*/; do
    if [ -d "$folder" ]; then
        all_folders+=("$folder")
    fi
done

if [ ${#all_folders[@]} -eq 0 ]; then
    echo "❌ ERROR: No BioProject folders (PRJ*) found in ${BASE_DIR}!"
    echo "    Expected folders: PRJNA*, PRJEB*, PRJCA*, PRJDB*"
    exit 1
fi

echo "Found ${#all_folders[@]} dataset folder(s)"
echo ""

# -----------------------------------------------------------------------
# Collect tables and representative sequences
# -----------------------------------------------------------------------
echo ">>> Step 3: Collecting tables and representative sequences..."

ALL_TABLES=()
ALL_REP_SEQS=()
SUCCESSFUL_DATASETS=()

for folder in "${all_folders[@]}"; do
    folder_name=$(basename "$folder")
    
    echo "Processing: $folder_name"
    
    # Find output files
    rep_files=("${folder}"/*-final-rep-seqs.qza)
    table_files=("${folder}"/*-final-table.qza)
    
    # Validate rep-seqs file
    if [ ${#rep_files[@]} -eq 0 ] || [ ! -f "${rep_files[0]}" ]; then
        echo "  ⚠️  WARNING: No rep-seqs file found, skipping..."
        continue
    fi
    
    # Validate table file
    if [ ${#table_files[@]} -eq 0 ] || [ ! -f "${table_files[0]}" ]; then
        echo "  ⚠️  WARNING: No table file found, skipping..."
        continue
    fi
    
    REP_SEQS="${rep_files[0]}"
    TABLE="${table_files[0]}"
    
    echo "  Rep-seqs: $(basename "$REP_SEQS")"
    echo "  Table: $(basename "$TABLE")"
    
    # Add to collection
    ALL_TABLES+=("$TABLE")
    ALL_REP_SEQS+=("$REP_SEQS")
    SUCCESSFUL_DATASETS+=("$folder_name")
    
    echo "  ✓ Collected"
    echo ""
done

# Validate collection
if [ ${#ALL_TABLES[@]} -eq 0 ]; then
    echo "❌ ERROR: No valid datasets were collected!"
    exit 1
fi

echo "Successfully collected ${#ALL_TABLES[@]} dataset(s)"
echo ""

# -----------------------------------------------------------------------
# Merge feature tables
# -----------------------------------------------------------------------
echo ">>> Step 4: Merging feature tables..."

MERGED_TABLE="${MERGED_DIR}/merged-table.qza"

if ! qiime feature-table merge \
    --i-tables "${ALL_TABLES[@]}" \
    --o-merged-table "$MERGED_TABLE" --verbose; then
    echo "❌ ERROR: Table merging failed!"
    exit 1
fi

echo "✓ Feature tables merged successfully"
echo ""

# -----------------------------------------------------------------------
# Merge representative sequences
# -----------------------------------------------------------------------
echo ">>> Step 5: Merging representative sequences..."

MERGED_REP_SEQS="${MERGED_DIR}/merged-rep-seqs.qza"

if ! qiime feature-table merge-seqs \
    --i-data "${ALL_REP_SEQS[@]}" \
    --o-merged-data "$MERGED_REP_SEQS" --verbose; then
    echo "❌ ERROR: Rep-seqs merging failed!"
    exit 1
fi

echo "✓ Representative sequences merged successfully"
echo ""

# -----------------------------------------------------------------------
# Generate summary statistics
# -----------------------------------------------------------------------
echo ">>> Step 6: Generating merged table summary..."

if qiime feature-table summarize \
    --i-table "$MERGED_TABLE" \
    --o-visualization "${MERGED_DIR}/merged-table-summary.qzv" --verbose; then
    echo "✓ Summary generated"
else
    echo "⚠️  WARNING: Failed to generate summary visualization"
fi
echo ""

# -----------------------------------------------------------------------
# GreenGenes2 annotation
# -----------------------------------------------------------------------
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
    echo "❌ ERROR: GreenGenes2 annotation failed!"
    exit 1
fi

echo "✓ GreenGenes2 mapping completed"
echo ""

# -----------------------------------------------------------------------
# Taxonomy assignment
# -----------------------------------------------------------------------
echo ">>> Step 8: Assigning taxonomy..."

MERGED_TAXONOMY="${MERGED_DIR}/merged-taxonomy.qza"

if ! qiime greengenes2 taxonomy-from-table \
    --i-reference-taxonomy "${TAXONOMY}" \
    --i-table "$ANNOTATED_TABLE" \
    --o-classification "$MERGED_TAXONOMY" --verbose; then
    echo "❌ ERROR: Taxonomy assignment failed!"
    exit 1
fi

echo "✓ Taxonomy assigned successfully"
echo ""

# -----------------------------------------------------------------------
# Generate annotated table summary
# -----------------------------------------------------------------------
echo ">>> Step 9: Generating annotated table summary..."

if qiime feature-table summarize \
    --i-table "$ANNOTATED_TABLE" \
    --o-visualization "${MERGED_DIR}/merged-table-gg2-summary.qzv" --verbose; then
    echo "✓ Annotated table summary generated"
else
    echo "⚠️  WARNING: Failed to generate annotated table summary"
fi
echo ""

# -----------------------------------------------------------------------
# Filter phylogenetic tree
# -----------------------------------------------------------------------
echo ">>> Step 10: Filtering phylogenetic tree..."

FILTERED_TREE="${MERGED_DIR}/final-tree.qza"

if ! qiime phylogeny filter-tree \
    --i-tree "$PHYLO" \
    --i-table "$ANNOTATED_TABLE" \
    --o-filtered-tree "$FILTERED_TREE" --verbose; then
    echo "❌ ERROR: Tree filtering failed!"
    exit 1
fi

echo "✓ Phylogenetic tree filtered successfully"
echo ""


################################################################################
#                           PIPELINE COMPLETION                                #
################################################################################

echo "========================================="
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "Finished: $(date)"
echo "========================================="
echo ""
echo "Output files located in: ${MERGED_DIR}"
echo ""
echo "Final outputs:"
echo "  - Merged table:          $(basename "$MERGED_TABLE")"
echo "  - Merged sequences:      $(basename "$MERGED_REP_SEQS")"
echo "  - Annotated table:       $(basename "$ANNOTATED_TABLE")"
echo "  - Annotated sequences:   $(basename "$ANNOTATED_REP_SEQS")"
echo "  - Taxonomy:              $(basename "$MERGED_TAXONOMY")"
echo "  - Filtered tree:         $(basename "$FILTERED_TREE")"
echo ""
echo "Processed datasets:        ${#SUCCESSFUL_DATASETS[@]}"
echo ""
echo "Log files:"
echo "  - Success log:  $success_log"
echo "  - Failed log:   $failed_log"
echo "  - Skipped log:  $skipped_log"
echo ""
echo "========================================="
echo "Review log files for detailed information"
echo "========================================="