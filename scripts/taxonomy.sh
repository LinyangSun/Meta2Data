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

# Required environment variables from ggCOMBO
for var in NB_CLASSIFIER REF_SEQS SEPP_REF DB_LABEL CONFIDENCE; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: \$$var is not set. taxonomy.sh must be called from ggCOMBO." >&2
        exit 1
    fi
done
for var in NB_CLASSIFIER REF_SEQS SEPP_REF; do
    if [[ ! -f "${!var}" ]]; then
        echo "ERROR: File not found for \$$var: ${!var}" >&2
        exit 1
    fi
done

echo "Using database paths:"
echo "  Classifier:    $(basename "$NB_CLASSIFIER")"
echo "  Ref Seqs:      $(basename "$REF_SEQS")"
echo "  SEPP Ref:      $(basename "$SEPP_REF")"
echo "  DB label:      ${DB_LABEL}"
echo "  Confidence:    ${CONFIDENCE}"
echo ""

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

# Clean previous run output for this DB to avoid QIIME2 "file already exists" errors
echo ">>> Step 1: Preparing output directories..."
mkdir -p "${FINAL_DIR}" "${TMP_DIR}" "${MERGED_DIR}"
# Remove only database-specific files from a previous run (preserve other DB's results)
rm -f "${MERGED_DIR}"/*-"${DB_LABEL}".qza "${MERGED_DIR}"/*-"${DB_LABEL}".qzv \
      "${MERGED_DIR}"/*-"${DB_LABEL}"-summary.qzv
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
rm -f "$MERGED_TABLE"

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
rm -f "$MERGED_REP_SEQS"

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
rm -f "${MERGED_DIR}/merged-table-summary.qzv"
if qiime feature-table summarize \
    --i-table "$MERGED_TABLE" \
    --o-visualization "${MERGED_DIR}/merged-table-summary.qzv" --verbose; then
    echo "✓ Summary generated"
else
    echo "⚠️  WARNING: Summary generation failed"
fi
echo ""

# Orient sequences against database reference sequences.
# This corrects reverse-complement sequences so all reads face the same
# direction before taxonomy assignment and tree building.
echo ">>> Step 7: Orienting sequences against reference..."
ORIENTED_REP_SEQS="${MERGED_DIR}/oriented-rep-seqs-${DB_LABEL}.qza"
UNMATCHED_REP_SEQS="${MERGED_DIR}/unmatched-rep-seqs-${DB_LABEL}.qza"

if ! qiime rescript orient-seqs \
    --i-sequences "$MERGED_REP_SEQS" \
    --i-reference-sequences "$REF_SEQS" \
    --p-threads "$cpu" \
    --o-oriented-seqs "$ORIENTED_REP_SEQS" \
    --o-unmatched-seqs "$UNMATCHED_REP_SEQS" --verbose; then
    echo "❌ ERROR: Sequence orientation failed"
    exit 1
fi

echo "✓ Sequences oriented"
echo ""

# Filter merged table to keep only features present in oriented sequences
echo ">>> Step 8: Filtering table to oriented features..."
ORIENTED_TABLE="${MERGED_DIR}/merged-table-oriented-${DB_LABEL}.qza"

if ! qiime feature-table filter-features \
    --i-table "$MERGED_TABLE" \
    --m-metadata-file "$ORIENTED_REP_SEQS" \
    --o-filtered-table "$ORIENTED_TABLE" --verbose; then
    echo "❌ ERROR: Table filtering by oriented sequences failed"
    exit 1
fi

echo "✓ Table filtered to oriented features"
echo ""

# Taxonomy assignment via pre-trained Naive Bayes classifier (sklearn)
# Uses the oriented sequences as input.
echo ">>> Step 9: Assigning taxonomy via pre-trained Naive Bayes classifier..."
echo "Using $cpu CPU threads, confidence=${CONFIDENCE}"
MERGED_TAXONOMY="${MERGED_DIR}/merged-taxonomy-${DB_LABEL}.qza"

if ! qiime feature-classifier classify-sklearn \
    --i-classifier "${NB_CLASSIFIER}" \
    --i-reads "$ORIENTED_REP_SEQS" \
    --p-n-jobs "$cpu" \
    --p-confidence "$CONFIDENCE" \
    --o-classification "$MERGED_TAXONOMY" --verbose; then
    echo "❌ ERROR: Taxonomy assignment failed"
    exit 1
fi

echo "✓ Taxonomy assigned via sklearn classifier (${DB_LABEL})"
echo ""

# SEPP fragment insertion — builds the phylogenetic tree by inserting query
# sequences into the GG2 reference tree (always GG2, regardless of --db-type).
# Uses oriented sequences as input.
echo ">>> Step 10: Building phylogenetic tree via SEPP fragment insertion..."

INSERTION_TREE="${MERGED_DIR}/insertion-tree-${DB_LABEL}.qza"
INSERTION_PLACEMENTS="${MERGED_DIR}/insertion-placements-${DB_LABEL}.qza"
SEPP_TREE_OK=false

if qiime fragment-insertion sepp \
    --i-representative-sequences "$ORIENTED_REP_SEQS" \
    --i-reference-database "$SEPP_REF" \
    --p-threads "$cpu" \
    --o-tree "$INSERTION_TREE" \
    --o-placements "$INSERTION_PLACEMENTS" --verbose; then

    SEPP_TREE_OK=true
    echo "✓ SEPP fragment insertion completed"

    # Filter oriented table to features that were placed in the tree.
    echo ">>> Step 11: Filtering table to tree-placed features..."
    if qiime fragment-insertion filter-features \
        --i-table "$ORIENTED_TABLE" \
        --i-tree "$INSERTION_TREE" \
        --o-filtered-table "${MERGED_DIR}/merged-table-tree-${DB_LABEL}.qza" \
        --o-removed-table "${MERGED_DIR}/merged-table-no-tree-${DB_LABEL}.qza" --verbose; then
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
echo ">>> Step 12: Generating summary visualizations..."
if [[ "$SEPP_TREE_OK" == true ]]; then
    if qiime feature-table summarize \
        --i-table "${MERGED_DIR}/merged-table-tree-${DB_LABEL}.qza" \
        --o-visualization "${MERGED_DIR}/merged-table-tree-${DB_LABEL}-summary.qzv" --verbose; then
        echo "✓ Tree-placed feature table summary generated"
    else
        echo "⚠️  WARNING: Tree-placed table summary generation failed"
    fi
fi
echo ""