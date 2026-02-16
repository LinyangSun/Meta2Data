set -euo pipefail

# Guard required environment variables
if [[ -z "${INPUT_DIR:-}" ]]; then
    echo "ERROR: \$INPUT_DIR is not set. taxonomy.sh must be called from ggCOMBO or AmpliconPIP." >&2
    exit 1
fi
if [[ -z "${OUTPUT:-}" ]]; then
    echo "ERROR: \$OUTPUT is not set. taxonomy.sh must be called from ggCOMBO or AmpliconPIP." >&2
    exit 1
fi
if [[ -z "${cpu:-}" ]]; then
    echo "ERROR: \$cpu is not set. taxonomy.sh must be called from ggCOMBO or AmpliconPIP." >&2
    exit 1
fi

################################################################################
#                          DATABASE RESOLUTION                                 #
################################################################################

echo "========================================="
echo "Resolving database..."
echo "========================================="

# If NB_CLASSIFIER, SEPP_REF are already set (e.g. from ggCOMBO), use them directly
if [[ -n "${NB_CLASSIFIER:-}" && -f "${NB_CLASSIFIER:-}" ]] && \
   [[ -n "${SEPP_REF:-}" && -f "${SEPP_REF:-}" ]]; then
    echo "Using pre-configured database paths:"
    echo "  Classifier:    $(basename "$NB_CLASSIFIER")"
    echo "  SEPP Ref:      $(basename "$SEPP_REF")"
    echo ""
else
    # Auto-detection fallback
    find_database_dir() {
        if [[ -n "${DB_DIR:-}" && -d "${DB_DIR:-}" ]]; then
            echo "$DB_DIR"
            return 0
        fi

        local COMMON_PATHS=(
            "${CONDA_PREFIX:-}/share/qiime2/data/greengenes2"
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
        echo "ERROR: Database not found"
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

    NB_CLASSIFIER=$(find_db_file "backbone.full-length.nb" "$DB_DIR")
    SEPP_REF=$(find_db_file "sepp-refs" "$DB_DIR")

    if [[ -z "$NB_CLASSIFIER" ]] || [[ -z "$SEPP_REF" ]]; then
        echo "ERROR: Required database files not found in $DB_DIR"
        exit 1
    fi

    echo "Classifier:    $(basename "$NB_CLASSIFIER")"
    echo "SEPP Ref:      $(basename "$SEPP_REF")"
    echo ""
fi

################################################################################
#                    PHASE 3: MERGE AND TAXONOMY ASSIGNMENT                    #
################################################################################

echo "========================================="
echo "PHASE 3: Merge & Taxonomy Assignment"
echo "Started: $(date)"
echo "TMPDIR:  ${TMPDIR:-/tmp (default)}"
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
for folder in "${INPUT_DIR}"/PRJ*/; do
    [ -d "$folder" ] || continue
    all_folders+=("$folder")
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

# Taxonomy assignment via pre-trained Naive Bayes classifier (sklearn)
# Uses the GG2 full-length backbone classifier which has both reference
# sequences and taxonomy built into the trained model.
echo ">>> Step 7: Assigning taxonomy via pre-trained Naive Bayes classifier..."
echo "Using $cpu CPU threads"
MERGED_TAXONOMY="${MERGED_DIR}/merged-taxonomy.qza"

if ! qiime feature-classifier classify-sklearn \
    --i-classifier "${NB_CLASSIFIER}" \
    --i-reads "$MERGED_REP_SEQS" \
    --p-n-jobs "$cpu" \
    --o-classification "$MERGED_TAXONOMY" --verbose; then
    echo "❌ ERROR: Taxonomy assignment failed"
    exit 1
fi

echo "✓ Taxonomy assigned via sklearn classifier"
echo ""

# SEPP fragment insertion — builds the phylogenetic tree by inserting query
# sequences into a reference tree (e.g. Greengenes 13.8).  Unlike GG2
# non-v4-16s this does NOT rename feature IDs, so the tree, table, and
# taxonomy all share the same original ASV hashes.
echo ">>> Step 8: Building phylogenetic tree via SEPP fragment insertion..."

INSERTION_TREE="${MERGED_DIR}/insertion-tree.qza"
INSERTION_PLACEMENTS="${MERGED_DIR}/insertion-placements.qza"
SEPP_TREE_OK=false

if qiime fragment-insertion sepp \
    --i-representative-sequences "$MERGED_REP_SEQS" \
    --i-reference-database "$SEPP_REF" \
    --p-threads "$cpu" \
    --o-tree "$INSERTION_TREE" \
    --o-placements "$INSERTION_PLACEMENTS" --verbose; then

    SEPP_TREE_OK=true
    echo "✓ SEPP fragment insertion completed"

    # Filter table to features that were placed in the tree.
    # Features that could not be inserted are separated into a removed table.
    echo ">>> Step 9: Filtering table to tree-placed features..."
    if qiime fragment-insertion filter-features \
        --i-table "$MERGED_TABLE" \
        --i-tree "$INSERTION_TREE" \
        --o-filtered-table "${MERGED_DIR}/merged-table-tree.qza" \
        --o-removed-table "${MERGED_DIR}/merged-table-no-tree.qza" --verbose; then
        echo "✓ Table filtered by tree placement"
    else
        echo "⚠️  WARNING: Feature filtering by tree failed."
    fi
else
    echo "⚠️  WARNING: SEPP fragment insertion failed."
    echo "   Taxonomy was already assigned via sklearn classifier — no reads lost."
    echo "   Only UniFrac / phylogenetic diversity analyses will be unavailable."
fi
echo ""

# Generate summary visualizations
echo ">>> Step 10: Generating summary visualizations..."
if [[ "$SEPP_TREE_OK" == true ]]; then
    if qiime feature-table summarize \
        --i-table "${MERGED_DIR}/merged-table-tree.qza" \
        --o-visualization "${MERGED_DIR}/merged-table-tree-summary.qzv" --verbose; then
        echo "✓ Tree-placed feature table summary generated"
    else
        echo "⚠️  WARNING: Tree-placed table summary generation failed"
    fi
fi
echo ""