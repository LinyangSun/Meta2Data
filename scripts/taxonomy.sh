################################################################################
#                          DATABASE RESOLUTION                                 #
################################################################################

echo "========================================="
echo "Resolving GreenGenes2 database..."
echo "========================================="

# If BACKBONE, TAXONOMY, PHYLO are already set (e.g. from ggCOMBO), use them directly
if [[ -n "$BACKBONE" && -f "$BACKBONE" ]] && \
   [[ -n "$TAXONOMY" && -f "$TAXONOMY" ]] && \
   [[ -n "$PHYLO" && -f "$PHYLO" ]]; then
    echo "Using pre-configured database paths:"
    echo "  Backbone:  $(basename "$BACKBONE")"
    echo "  Taxonomy:  $(basename "$TAXONOMY")"
    echo "  Phylogeny: $(basename "$PHYLO")"
    echo ""
else
    # Auto-detection fallback
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

    DB_DIR=$(find_database_dir)

    if [[ -z "$DB_DIR" ]]; then
        echo "ERROR: GreenGenes2 database not found"
        echo ""
        echo "Searched locations:"
        echo "  - \$CONDA_PREFIX/share/qiime2/data/greengenes2"
        echo "  - \$HOME/.qiime2/db/greengenes2"
        echo "  - /scratch/project_2009135/db/gg2"
        echo ""
        echo "Use 'Meta2Data ggCOMBO --db /path/to/db' to specify manually"
        exit 1
    fi

    echo "Database found: $DB_DIR"

    BACKBONE=$(find_db_file "backbone" "$DB_DIR")
    TAXONOMY=$(find_db_file "taxonomy" "$DB_DIR")
    PHYLO=$(find_db_file "phylogeny" "$DB_DIR")

    if [[ -z "$BACKBONE" ]] || [[ -z "$TAXONOMY" ]] || [[ -z "$PHYLO" ]]; then
        echo "ERROR: Required database files not found in $DB_DIR"
        exit 1
    fi

    echo "Backbone:  $(basename "$BACKBONE")"
    echo "Taxonomy:  $(basename "$TAXONOMY")"
    echo "Phylogeny: $(basename "$PHYLO")"
    echo ""
fi

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

# Clean previous run output to avoid QIIME2 "file already exists" errors
echo ">>> Step 1: Preparing output directories..."
if [ -d "${MERGED_DIR}" ]; then
    echo "Removing previous merged output: ${MERGED_DIR}"
    rm -rf "${MERGED_DIR}"
fi
mkdir -p "${FINAL_DIR}" "${TMP_DIR}" "${MERGED_DIR}"
echo "✓ Directories ready"
echo ""

# Collect dataset folders
echo ">>> Step 2: Collecting dataset outputs..."

all_folders=()
for folder in "${OUTPUT}"/*/; do
    [ -d "$folder" ] || continue
    folder_base=$(basename "$folder")
    # Match BioProject accessions (PRJ*) and CNCB/GSA accessions (CRA*)
    if [[ "$folder_base" == PRJ* || "$folder_base" == CRA* ]]; then
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