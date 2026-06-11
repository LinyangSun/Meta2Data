set -euo pipefail

# Guard required environment variables (set by Meta2Data-AmpliconTAXA)
for var in MODE INPUT_DIR OUTPUT cpu NB_CLASSIFIER ORIENT_REF DB_LABEL CONFIDENCE; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: \$$var is not set. taxonomy.sh must be called from AmpliconTAXA." >&2
        exit 1
    fi
done

case "$MODE" in
    asv|otu) ;;
    *) echo "ERROR: \$MODE must be 'asv' or 'otu' (got '$MODE')" >&2; exit 1 ;;
esac

# Files that must exist on disk
NEEDED_DB=("$NB_CLASSIFIER" "$ORIENT_REF")
if [[ "$MODE" == "otu" ]]; then
    if [[ -z "${SEPP_REF:-}" ]]; then
        echo "ERROR: \$SEPP_REF is not set (required for --otu)." >&2
        exit 1
    fi
    NEEDED_DB+=("$SEPP_REF")
fi
for f in "${NEEDED_DB[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Database file not found: $f" >&2
        exit 1
    fi
done

# A .qza/.qzv is reusable only if it is a non-empty, valid zip archive.
qza_ok() { [[ -s "$1" ]] && unzip -tq "$1" >/dev/null 2>&1; }

# Reuse-or-build a feature-table summary visualization (non-fatal on failure).
# Args: <in-table> <out-viz> <warn-msg> [guard-qza]
# If guard-qza is given, the summary is only attempted when that artifact is OK.
summarize_table() {
    local in_table="$1" out_viz="$2" warn_msg="$3" guard="${4:-}"
    if qza_ok "$out_viz"; then
        echo "✓ reuse: $(basename "$out_viz")"
    elif { [[ -z "$guard" ]] || qza_ok "$guard"; } && qiime feature-table summarize \
            --i-table "$in_table" \
            --o-visualization "$out_viz" --verbose; then
        echo "✓ Summary generated"
    else
        echo "$warn_msg"
    fi
}

echo "========================================="
echo "PHASE 3: Merge, Orient, Classify & Tree (${MODE})"
echo "Started: $(date)"
echo "TMPDIR:  ${TMPDIR:-/tmp (default)}"
echo "Database (taxonomy): ${DB_LABEL}"
echo "Orient reference:    $(basename "$ORIENT_REF")  (GG2, always)"
echo "Confidence:          ${CONFIDENCE}"
echo "========================================="
echo ""

################################################################################
#                        OUTPUT LAYOUT (mode-scoped)                           #
################################################################################
# Final products live directly under final-<MODE>/ ; every intermediate and
# diagnostic artifact goes under final-<MODE>/tmp/ (safe to delete; deleting it
# only forces a recompute on the next run). Only the taxonomy artifact depends
# on --db-type, so re-running a different DB reuses everything else.

FINAL_DIR="${OUTPUT}/final-${MODE}"
TMP="${FINAL_DIR}/tmp"
mkdir -p "$FINAL_DIR" "$TMP"

# Shared (database-independent) reuse anchors — kept in tmp/
MERGED_TABLE="${TMP}/mergedTable.qza"
MERGED_REP_SEQS="${TMP}/mergedRepSeqs.qza"
MERGED_SUMMARY="${TMP}/mergedTableSummary.qzv"
UNMATCHED="${TMP}/unmatchedRepSeqs.qza"

# Oriented outputs: final products for --asv, intermediates for --otu
if [[ "$MODE" == "asv" ]]; then
    ORIENTED_REP_SEQS="${FINAL_DIR}/orientedRepSeqs.qza"
    ORIENTED_TABLE="${FINAL_DIR}/orientedTable.qza"
    ORIENTED_TABLE_SUMMARY="${TMP}/orientedTableSummary.qzv"
else
    ORIENTED_REP_SEQS="${TMP}/orientedRepSeqs.qza"
    ORIENTED_TABLE="${TMP}/orientedTable.qza"
fi

# Taxonomy artifact carries the DB label (only db-specific output)
TAXONOMY="${FINAL_DIR}/${DB_LABEL}Taxonomy.qza"

# OTU tree products
SEPP_TREE="${FINAL_DIR}/seppTree.qza"
SEPP_PLACEMENTS="${TMP}/seppPlacements.qza"
TREE_TABLE="${FINAL_DIR}/treeFilteredTable.qza"
TREE_REP_SEQS="${FINAL_DIR}/treeFilteredRepSeqs.qza"
NOTREE_TABLE="${TMP}/treeUnplacedTable.qza"
TREE_TABLE_SUMMARY="${TMP}/treeFilteredTableSummary.qzv"

# ASV de novo tree products
DENOVO_ALN="${TMP}/denovoAlignment.qza"
DENOVO_MASKED="${TMP}/denovoMaskedAlignment.qza"
DENOVO_UNROOTED="${TMP}/denovoUnrootedTree.qza"
DENOVO_ROOTED="${FINAL_DIR}/denovoRootedTree.qza"

################################################################################
#                    STEP 1: COLLECT DATASET OUTPUTS                            #
################################################################################

echo ">>> Step 1: Collecting ${MODE} dataset outputs..."

all_folders=()
for folder in "${INPUT_DIR}"/PRJ*/; do
    [ -d "$folder" ] || continue
    all_folders+=("$folder")
done

if [ ${#all_folders[@]} -eq 0 ]; then
    echo "❌ ERROR: No PRJ* dataset folders found in ${INPUT_DIR}"
    exit 1
fi

table_glob="*-${MODE}-final-table.qza"
rep_glob="*-${MODE}-final-rep-seqs.qza"

ALL_TABLES=()
ALL_REP_SEQS=()

for folder in "${all_folders[@]}"; do
    folder_name=$(basename "$folder")

    rep_files=("${folder}"/$rep_glob)
    table_files=("${folder}"/$table_glob)

    if [ ! -f "${rep_files[0]}" ] || [ ! -f "${table_files[0]}" ]; then
        echo "  ⚠️  ${folder_name}: no ${MODE} outputs, skipping"
        continue
    fi

    ALL_TABLES+=("${table_files[0]}")
    ALL_REP_SEQS+=("${rep_files[0]}")
    echo "  ✓ ${folder_name}"
done

if [ ${#ALL_TABLES[@]} -eq 0 ]; then
    echo "❌ ERROR: No ${MODE} datasets collected (looked for ${table_glob})"
    exit 1
fi

echo "Collected ${#ALL_TABLES[@]} ${MODE} dataset(s)"
echo ""

################################################################################
#               STEP 2: MERGE (database-independent, reusable)                 #
################################################################################

echo ">>> Step 2: Merging feature tables..."
if qza_ok "$MERGED_TABLE"; then
    echo "✓ reuse: $(basename "$MERGED_TABLE")"
else
    rm -f "$MERGED_TABLE"
    if ! qiime feature-table merge \
        --i-tables "${ALL_TABLES[@]}" \
        --o-merged-table "$MERGED_TABLE" --verbose; then
        echo "❌ ERROR: Table merging failed"
        exit 1
    fi
fi

echo ">>> Step 3: Merging representative sequences..."
if qza_ok "$MERGED_REP_SEQS"; then
    echo "✓ reuse: $(basename "$MERGED_REP_SEQS")"
else
    rm -f "$MERGED_REP_SEQS"
    if ! qiime feature-table merge-seqs \
        --i-data "${ALL_REP_SEQS[@]}" \
        --o-merged-data "$MERGED_REP_SEQS" --verbose; then
        echo "❌ ERROR: Sequence merging failed"
        exit 1
    fi
fi

echo ">>> Step 4: Merged table summary..."
summarize_table "$MERGED_TABLE" "$MERGED_SUMMARY" \
    "⚠️  WARNING: Merged summary generation failed (non-fatal)"
echo ""

################################################################################
#        STEP 5: ORIENT vs GG2 (database-independent, reusable)                #
################################################################################
# Orientation always uses the GG2 backbone (ORIENT_REF), regardless of
# --db-type. Sequences that cannot be oriented (no match) go to UNMATCHED and
# are dropped from the analysis.

echo ">>> Step 5: Orienting sequences against GG2 backbone..."
if qza_ok "$ORIENTED_REP_SEQS"; then
    echo "✓ reuse: $(basename "$ORIENTED_REP_SEQS")"
else
    rm -f "$ORIENTED_REP_SEQS" "$UNMATCHED"
    if ! qiime rescript orient-seqs \
        --i-sequences "$MERGED_REP_SEQS" \
        --i-reference-sequences "$ORIENT_REF" \
        --p-threads "$cpu" \
        --o-oriented-seqs "$ORIENTED_REP_SEQS" \
        --o-unmatched-seqs "$UNMATCHED" --verbose; then
        echo "❌ ERROR: Sequence orientation failed"
        exit 1
    fi
fi

echo ">>> Step 6: Filtering table to oriented features..."
if qza_ok "$ORIENTED_TABLE"; then
    echo "✓ reuse: $(basename "$ORIENTED_TABLE")"
else
    rm -f "$ORIENTED_TABLE"
    if ! qiime feature-table filter-features \
        --i-table "$MERGED_TABLE" \
        --m-metadata-file "$ORIENTED_REP_SEQS" \
        --o-filtered-table "$ORIENTED_TABLE" --verbose; then
        echo "❌ ERROR: Table filtering by oriented sequences failed"
        exit 1
    fi
fi
echo ""

################################################################################
#        STEP 6b (ASV only): single-region length-variance guard               #
################################################################################
# ASV mode assumes a single amplicon region. If feature lengths vary widely the
# merged set probably mixes regions (e.g. full-length + V-region), which makes
# the de novo alignment/tree unreliable. Warn but do not stop.

if [[ "$MODE" == "asv" ]]; then
    echo ">>> Step 6b: ASV single-region length check..."
    LEN_PROBE="${TMP}/.lenprobe"
    rm -rf "$LEN_PROBE"; mkdir -p "$LEN_PROBE"
    if qiime tools export --input-path "$ORIENTED_REP_SEQS" \
            --output-path "$LEN_PROBE" >/dev/null 2>&1 \
       && [[ -s "${LEN_PROBE}/dna-sequences.fasta" ]]; then
        awk '/^>/{if(s!=""){print length(s)}; s=""; next}{s=s $0}END{if(s!="")print length(s)}' \
            "${LEN_PROBE}/dna-sequences.fasta" | sort -n > "${LEN_PROBE}/lengths.txt"
        n=$(wc -l < "${LEN_PROBE}/lengths.txt")
        if [[ "$n" -ge 20 ]]; then
            p5=$(awk -v n="$n" 'NR==int(n*0.05)+1{print; exit}' "${LEN_PROBE}/lengths.txt")
            p95=$(awk -v n="$n" 'NR>=int(n*0.95){print; exit}' "${LEN_PROBE}/lengths.txt")
            ratio=$(awk -v a="$p95" -v b="$p5" 'BEGIN{if(b>0) printf "%.2f", a/b; else print "0"}')
            echo "    feature length p5=${p5}bp  p95=${p95}bp  (p95/p5=${ratio}, n=${n})"
            if awk -v r="$ratio" 'BEGIN{exit !(r>1.5)}'; then
                echo "⚠️  WARNING: ASV feature lengths vary widely (p95/p5=${ratio} > 1.5)."
                echo "   ASV mode assumes a SINGLE amplicon region; mixing regions"
                echo "   (e.g. full-length + V-region) makes the de novo tree unreliable."
                echo "   Proceeding anyway (non-fatal)."
            fi
        else
            echo "    too few features (${n}) for a meaningful length check, skipping"
        fi
    else
        echo "⚠️  WARNING: could not export sequences for length check (non-fatal)"
    fi
    rm -rf "$LEN_PROBE"
    echo ""
fi

################################################################################
#        STEP 7: TAXONOMY ASSIGNMENT (only database-dependent step)            #
################################################################################

echo ">>> Step 7: Assigning taxonomy (${DB_LABEL}, confidence=${CONFIDENCE})..."
if qza_ok "$TAXONOMY"; then
    echo "✓ reuse: $(basename "$TAXONOMY")"
else
    rm -f "$TAXONOMY"
    if ! qiime feature-classifier classify-sklearn \
        --i-classifier "$NB_CLASSIFIER" \
        --i-reads "$ORIENTED_REP_SEQS" \
        --p-n-jobs "$cpu" \
        --p-confidence "$CONFIDENCE" \
        --o-classification "$TAXONOMY" --verbose; then
        echo "❌ ERROR: Taxonomy assignment failed"
        exit 1
    fi
fi
echo "✓ Taxonomy: $(basename "$TAXONOMY")"
echo ""

################################################################################
#                    STEP 8: PHYLOGENETIC TREE                                  #
################################################################################

if [[ "$MODE" == "otu" ]]; then
    # ---- OTU: SEPP fragment insertion into the GG2 reference tree ----
    echo ">>> Step 8: Building phylogenetic tree via SEPP fragment insertion..."
    SEPP_TREE_OK=false
    if qza_ok "$SEPP_TREE"; then
        echo "✓ reuse: $(basename "$SEPP_TREE")"
        SEPP_TREE_OK=true
    else
        rm -f "$SEPP_TREE" "$SEPP_PLACEMENTS"
        if qiime fragment-insertion sepp \
            --i-representative-sequences "$ORIENTED_REP_SEQS" \
            --i-reference-database "$SEPP_REF" \
            --p-threads "$cpu" \
            --o-tree "$SEPP_TREE" \
            --o-placements "$SEPP_PLACEMENTS" --verbose; then
            SEPP_TREE_OK=true
            echo "✓ SEPP fragment insertion completed"
        else
            echo "⚠️  WARNING: SEPP fragment insertion failed."
            echo "   Taxonomy was already assigned — no reads lost."
            echo "   Only UniFrac / phylogenetic diversity will be unavailable."
        fi
    fi

    if [[ "$SEPP_TREE_OK" == true ]]; then
        echo ">>> Step 9: Filtering table to tree-placed features..."
        if qza_ok "$TREE_TABLE"; then
            echo "✓ reuse: $(basename "$TREE_TABLE")"
        else
            rm -f "$TREE_TABLE" "$NOTREE_TABLE"
            if qiime fragment-insertion filter-features \
                --i-table "$ORIENTED_TABLE" \
                --i-tree "$SEPP_TREE" \
                --o-filtered-table "$TREE_TABLE" \
                --o-removed-table "$NOTREE_TABLE" --verbose; then
                echo "✓ Table filtered by tree placement"
            else
                echo "⚠️  WARNING: Feature filtering by tree failed."
            fi
        fi

        echo ">>> Step 10: Filtering representative sequences to tree-placed features..."
        if qza_ok "$TREE_REP_SEQS"; then
            echo "✓ reuse: $(basename "$TREE_REP_SEQS")"
        elif qza_ok "$TREE_TABLE"; then
            rm -f "$TREE_REP_SEQS"
            if qiime feature-table filter-seqs \
                --i-data "$ORIENTED_REP_SEQS" \
                --i-table "$TREE_TABLE" \
                --o-filtered-data "$TREE_REP_SEQS" --verbose; then
                echo "✓ Representative sequences filtered by tree placement"
            else
                echo "⚠️  WARNING: Rep-seqs filtering by tree failed."
            fi
        fi

        echo ">>> Step 11: Tree-placed table summary..."
        summarize_table "$TREE_TABLE" "$TREE_TABLE_SUMMARY" \
            "⚠️  WARNING: Tree-placed table summary failed (non-fatal)" "$TREE_TABLE"
    fi
else
    # ---- ASV: de novo tree (mafft -> mask -> fasttree) ----
    echo ">>> Step 8: ASV final table summary..."
    summarize_table "$ORIENTED_TABLE" "$ORIENTED_TABLE_SUMMARY" \
        "⚠️  WARNING: ASV table summary failed (non-fatal)"

    echo ">>> Step 9: Building de novo tree (mafft + mask + fasttree)..."
    if qza_ok "$DENOVO_ROOTED"; then
        echo "✓ reuse: $(basename "$DENOVO_ROOTED")"
    else
        rm -f "$DENOVO_ALN" "$DENOVO_MASKED" "$DENOVO_UNROOTED" "$DENOVO_ROOTED"
        if ! qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences "$ORIENTED_REP_SEQS" \
            --p-n-threads "$cpu" \
            --o-alignment "$DENOVO_ALN" \
            --o-masked-alignment "$DENOVO_MASKED" \
            --o-tree "$DENOVO_UNROOTED" \
            --o-rooted-tree "$DENOVO_ROOTED" --verbose; then
            echo "❌ ERROR: De novo tree building failed"
            exit 1
        fi
        echo "✓ De novo tree built"
    fi
fi
echo ""

################################################################################
#                             FINAL SUMMARY                                     #
################################################################################

echo "========================================="
echo "Pipeline complete (${MODE}, ${DB_LABEL})"
echo "Finished: $(date)"
echo "========================================="
echo "Final products in: ${FINAL_DIR}/"
if [[ "$MODE" == "otu" ]]; then
    echo "  $(basename "$TREE_TABLE")"
    echo "  $(basename "$TREE_REP_SEQS")"
    echo "  $(basename "$SEPP_TREE")"
else
    echo "  $(basename "$ORIENTED_TABLE")"
    echo "  $(basename "$ORIENTED_REP_SEQS")"
    echo "  $(basename "$DENOVO_ROOTED")"
fi
echo "  $(basename "$TAXONOMY")"
echo "Intermediates/diagnostics in: ${TMP}/ (safe to delete)"
echo "========================================="
