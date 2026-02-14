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
echo "Resolving GreenGenes2 database..."
echo "========================================="

# If BACKBONE, BACKBONE_TAX, TAXONOMY, PHYLO are already set (e.g. from ggCOMBO), use them directly
if [[ -n "$BACKBONE" && -f "$BACKBONE" ]] && \
   [[ -n "${BACKBONE_TAX:-}" && -f "${BACKBONE_TAX:-}" ]] && \
   [[ -n "$TAXONOMY" && -f "$TAXONOMY" ]] && \
   [[ -n "$PHYLO" && -f "$PHYLO" ]]; then
    echo "Using pre-configured database paths:"
    echo "  Backbone:      $(basename "$BACKBONE")"
    echo "  Backbone Tax:  $(basename "$BACKBONE_TAX")"
    echo "  Taxonomy:      $(basename "$TAXONOMY")"
    echo "  Phylogeny:     $(basename "$PHYLO")"
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

    BACKBONE=$(find_db_file "backbone.full-length" "$DB_DIR")
    BACKBONE_TAX=$(find_db_file "backbone.tax" "$DB_DIR")
    TAXONOMY=$(find_db_file "taxonomy" "$DB_DIR")
    PHYLO=$(find_db_file "phylogeny" "$DB_DIR")

    if [[ -z "$BACKBONE" ]] || [[ -z "$BACKBONE_TAX" ]] || [[ -z "$TAXONOMY" ]] || [[ -z "$PHYLO" ]]; then
        echo "ERROR: Required database files not found in $DB_DIR"
        exit 1
    fi

    echo "Backbone:      $(basename "$BACKBONE")"
    echo "Backbone Tax:  $(basename "$BACKBONE_TAX")"
    echo "Taxonomy:      $(basename "$TAXONOMY")"
    echo "Phylogeny:     $(basename "$PHYLO")"
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

# Taxonomy assignment via vsearch consensus classifier
# This replaces the previous GG2 non-v4-16s + taxonomy-from-table approach.
# vsearch classifies ALL features without dropping any reads, regardless of
# which variable region (V1-V9, full-length) the sequences come from.
echo ">>> Step 7: Assigning taxonomy via vsearch consensus classifier..."
echo "Using $cpu CPU threads"
MERGED_TAXONOMY="${MERGED_DIR}/merged-taxonomy.qza"

if ! qiime feature-classifier classify-consensus-vsearch \
    --i-query "$MERGED_REP_SEQS" \
    --i-reference-reads "${BACKBONE}" \
    --i-reference-taxonomy "${BACKBONE_TAX}" \
    --p-threads "$cpu" \
    --p-perc-identity 0.8 \
    --p-maxaccepts 3 \
    --o-classification "$MERGED_TAXONOMY" \
    --o-search-results "${MERGED_DIR}/vsearch-hits.qza" --verbose; then
    echo "❌ ERROR: Taxonomy assignment failed"
    exit 1
fi

echo "✓ Taxonomy assigned (all features retained)"
echo ""

# GG2 backbone mapping — used ONLY for building the phylogenetic tree.
# This step may drop features that cannot be placed on the GG2 backbone.
# That is acceptable: the tree is only needed for UniFrac / phylogenetic
# diversity analyses.  Taxonomy has already been assigned above.
echo ">>> Step 8: Mapping to GG2 backbone for phylogenetic tree..."

ANNOTATED_TABLE="${MERGED_DIR}/merged-table-gg2.qza"
ANNOTATED_REP_SEQS="${MERGED_DIR}/merged-rep-seqs-gg2.qza"
FILTERED_TREE="${MERGED_DIR}/final-tree.qza"
GG2_TREE_OK=false

if qiime greengenes2 non-v4-16s \
    --i-table "$MERGED_TABLE" \
    --i-sequences "$MERGED_REP_SEQS" \
    --p-threads "$cpu" \
    --i-backbone "${BACKBONE}" \
    --o-mapped-table "$ANNOTATED_TABLE" \
    --o-representatives "$ANNOTATED_REP_SEQS" --verbose; then

    echo "✓ GG2 backbone mapping completed"

    # Filter phylogenetic tree to mapped features
    echo ">>> Step 9: Filtering phylogenetic tree..."
    if qiime phylogeny filter-tree \
        --i-tree "$PHYLO" \
        --i-table "$ANNOTATED_TABLE" \
        --o-filtered-tree "$FILTERED_TREE" --verbose; then
        GG2_TREE_OK=true
        echo "✓ Phylogenetic tree filtered"
    else
        echo "⚠️  WARNING: Tree filtering failed. UniFrac analyses will not be available."
    fi
else
    echo "⚠️  WARNING: GG2 backbone mapping failed."
    echo "   Taxonomy was already assigned via vsearch — no reads lost."
    echo "   Only UniFrac / phylogenetic diversity analyses will be unavailable."
fi
echo ""

# Generate summary visualizations
echo ">>> Step 10: Generating summary visualizations..."
if [[ "$GG2_TREE_OK" == true ]]; then
    if qiime feature-table summarize \
        --i-table "$ANNOTATED_TABLE" \
        --o-visualization "${MERGED_DIR}/merged-table-gg2-summary.qzv" --verbose; then
        echo "✓ GG2-mapped table summary generated"
    else
        echo "⚠️  WARNING: GG2 summary generation failed"
    fi
fi
echo ""