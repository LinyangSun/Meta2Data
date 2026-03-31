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
    -t, --threads INT          Total CPU threads available (default: 4)
                               Automatically split across parallel datasets:
                               per-dataset threads = threads ÷ max-parallel
    --max-parallel INT         Number of datasets to process in parallel (default: 2)
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
MAX_PARALLEL=2
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
        --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
        --col-bioproject) COL_BIOPROJECT="$2"; shift 2 ;;
        --col-sra) COL_SRA="$2"; shift 2 ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Error: Unknown option '$1'"; show_help; exit 1 ;;
    esac
done

# Validate thread count is a positive integer
if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: --threads must be a positive integer, got '$THREADS'"
    exit 1
fi

# Validate max-parallel count
if ! [[ "$MAX_PARALLEL" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: --max-parallel must be a positive integer, got '$MAX_PARALLEL'"
    exit 1
fi

# Strip trailing slashes from paths
METADATA="${METADATA%/}"
OUTPUT="${OUTPUT%/}"

# Validate input
if [[ -z "$METADATA" ]]; then
    echo "Error: --metadata is required"
    exit 1
fi

if [[ ! -f "$METADATA" ]]; then
    echo "Error: Metadata file not found: '$METADATA'. Please check the path."
    exit 1
fi

if [[ -z "$OUTPUT" ]]; then
    OUTPUT=$(dirname "$METADATA")
fi

# Compute per-dataset thread count: total threads ÷ max parallel datasets
THREADS_PER_DATASET=$(( THREADS / MAX_PARALLEL ))
if [[ "$THREADS_PER_DATASET" -lt 1 ]]; then
    THREADS_PER_DATASET=1
fi
export THREADS_PER_DATASET
export cpu=$THREADS_PER_DATASET

# Source function library (AmpliconFunction.sh contains all processing functions)
if [[ -f "${SCRIPTS}/AmpliconFunction.sh" ]]; then
    source "${SCRIPTS}/AmpliconFunction.sh"
else
    echo "Error: AmpliconFunction.sh not found at ${SCRIPTS}/AmpliconFunction.sh"
    exit 1
fi

cd "$OUTPUT" || exit 1
mkdir -p "${OUTPUT}/logs"

################################################################################
#                          PHASE 1: DATASET PREPARATION                        #
################################################################################

echo "========================================="
echo "PHASE 1: Dataset Preparation"
echo "Started: $(date)"
echo "========================================="

# Step 1: Generate dataset ID list
if ! python "${SCRIPTS}/py_16s.py" GenerateDatasetsIDsFile --FilePath "$METADATA" --Bioproject "$COL_BIOPROJECT" --OutputDir "$OUTPUT"; then
    echo "❌ ERROR: Failed to generate dataset IDs, please check your metadata file and column names for BioProject."
    exit 1
fi

mapfile -t Dataset_ID_sets < <(awk '{print $1}' "${OUTPUT}/datasets_ID.txt")

if [ ${#Dataset_ID_sets[@]} -eq 0 ]; then
    echo "❌ ERROR: No datasets found, please check your metadata file and column names for BioProject."
    exit 1
fi

# Step 2: Generate SRA file list
if ! python "${SCRIPTS}/py_16s.py" GenerateSRAsFile --FilePath "$METADATA" --Bioproject "$COL_BIOPROJECT" --SRA_Number "$COL_SRA" --OutputDir "$OUTPUT"; then
    echo "❌ ERROR: Failed to generate SRA file lists, please check your metadata file and column names for SRA."
    exit 1
fi

################################################################################
#                   PHASE 1.5: PLATFORM PRE-DETECTION                          #
################################################################################
# Batch-detect sequencing platforms BEFORE parallel processing to avoid
# NCBI Entrez API rate limits (3 req/s without API key).

echo "========================================="
echo "PHASE 1.5: Platform Detection"
echo "Started: $(date)"
echo "========================================="

export PLATFORM_CACHE_FILE="${OUTPUT}/.platform_cache.txt"
_pairs_file="${OUTPUT}/.platform_query_pairs.txt"
: > "$_pairs_file"

# Collect dataset_id<TAB>first_srr[<TAB>bioproject_id] for each unprocessed dataset
for _ds_id in "${Dataset_ID_sets[@]}"; do
    _ds_path="${OUTPUT}/${_ds_id}"

    # Skip already processed
    [[ -f "${_ds_path}/${_ds_id}-final-rep-seqs.qza" ]] && continue

    _sra_file="${_ds_path}/${_ds_id}_sra.txt"
    if [[ ! -f "$_sra_file" ]]; then
        echo "  Warning: SRA file not found for ${_ds_id}, skipping"
        continue
    fi

    _first_srr=$(awk 'NR==1 {print $1}' "$_sra_file")
    # CRR accessions need bioproject_id for CNCB API
    if [[ "$_first_srr" =~ ^CRR ]]; then
        printf '%s\t%s\t%s\n' "$_ds_id" "$_first_srr" "$_ds_id" >> "$_pairs_file"
    else
        printf '%s\t%s\n' "$_ds_id" "$_first_srr" >> "$_pairs_file"
    fi
done

# Single Python call: batch Entrez for NCBI, serial for CNCB
: > "$PLATFORM_CACHE_FILE"
if [[ -s "$_pairs_file" ]]; then
    _n_queries=$(wc -l < "$_pairs_file" | tr -d ' ')
    echo "  Querying platforms for ${_n_queries} datasets..."
    python "${SCRIPTS}/py_16s.py" batch_get_sequencing_platforms --pairs_file "$_pairs_file" > "$PLATFORM_CACHE_FILE"
    echo "  Resolved $(wc -l < "$PLATFORM_CACHE_FILE" | tr -d ' ')/${_n_queries} datasets"
fi

rm -f "$_pairs_file"
echo ""

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
low_quality_log="${OUTPUT}/low_quality_datasets.log"
summary_csv="${OUTPUT}/summary.csv"

: > "$failed_log"
: > "$success_log"
: > "$skipped_log"
: > "$low_quality_log"

echo "Threads: $THREADS total, $MAX_PARALLEL parallel datasets, $THREADS_PER_DATASET threads per dataset"

running_jobs=0
exec 3>&1   # save console stdout for milestone messages during parallel mode
_pipeline_start=$(date +%s)

for i in "${!Dataset_ID_sets[@]}"; do
    dataset_ID="${Dataset_ID_sets[$i]}"
    dataset_path="${OUTPUT}/${dataset_ID}"
    sra_file_name="${dataset_ID}_sra.txt"
    log_file="${OUTPUT}/logs/${dataset_ID}.log"
    platform="Unknown"

    echo "----------------------------------------"
    echo "Dataset $((i+1))/${#Dataset_ID_sets[@]}: $dataset_ID"
    echo "  Log: logs/${dataset_ID}.log"

    # Check if already processed
    if [ -f "${dataset_path}/${dataset_ID}-final-rep-seqs.qza" ]; then
        echo "✓ Already processed. Skipping."
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - ALREADY_DONE" >> "$skipped_log"
        continue
    fi

    # Wait for a slot if running at max parallel capacity
    if [[ "$MAX_PARALLEL" -gt 1 ]] && [[ "$running_jobs" -ge "$MAX_PARALLEL" ]]; then
        wait -n 2>/dev/null || true
        running_jobs=$((running_jobs - 1))
    fi

    _process_one_dataset() {
    local _ds_start=$(date +%s)
    (
        set -e
        cd "$dataset_path"

        # 1. Platform Detection — use pre-detected cache, fallback to API query
        echo ">>> Detecting sequencing platform..."
        first_srr=$(awk 'NR==1 {print $1}' "${sra_file_name}")

        # Read from cache file (written by Phase 1.5)
        platform=""
        if [[ -f "$PLATFORM_CACHE_FILE" ]]; then
            platform=$(awk -v id="$dataset_ID" '$1 == id {print $2}' "$PLATFORM_CACHE_FILE")
        fi

        # Fallback: query API if cache miss
        if [[ -z "$platform" ]]; then
            echo "  Cache miss, querying API..."
            if [[ "$first_srr" =~ ^CRR ]]; then
                echo "  CNCB accession detected, using BioProject: $dataset_ID"
                platform=$(python "${SCRIPTS}/py_16s.py" get_sequencing_platform --srr_id "$first_srr" --bioproject_id "$dataset_ID")
            else
                platform=$(python "${SCRIPTS}/py_16s.py" get_sequencing_platform --srr_id "$first_srr")
            fi
        fi

        echo "Detected platform: $platform"
        echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [1/3] Platform: $platform" >&3

        # 2. Platform-specific pipeline (download + processing)
        export dataset_path
        export dataset_name="$dataset_ID"
        ori_fastq_path="${dataset_path}/ori_fastq"

        if [[ "$platform" == "ILLUMINA" ]]; then
            fastp_path="${dataset_path}/tmp/step_02_fastp"
            n_srr=$(wc -l < "${sra_file_name}" | tr -d ' ')
            quality_cache="${dataset_path}/${dataset_ID}_quality_status.txt"

            # ── Resume checkpoint: check if fastp data (primers already removed) is intact ──
            n_fastp_fq=0
            if [[ -d "$fastp_path" ]]; then
                n_fastp_fq=$(find "$fastp_path" -type f -name '*.fastq*' | wc -l)
            fi

            if [[ "$n_fastp_fq" -gt 0 ]] && \
               [[ "$n_fastp_fq" -eq "$n_srr" || "$n_fastp_fq" -eq $(( n_srr * 2 )) ]]; then
                # fastp data is intact — resume from here
                echo ">>> Resuming: found $n_fastp_fq fastp files for $n_srr SRR accessions"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Resuming from checkpoint" >&3

                # Clean up downstream intermediate directories (keep step_02c as checkpoint)
                rm -rf "${dataset_path}/tmp/step_03_qza_import"
                rm -rf "${dataset_path}/tmp/step_04_qza_import_QualityFilter"
                rm -rf "${dataset_path}/tmp/step_05_dedupicate"
                rm -rf "${dataset_path}/tmp/step_05_denoise"
                rm -rf "${dataset_path}/tmp/step_06_vsearch_cli"
                rm -rf "${dataset_path}/tmp/step_06_ChimerasRemoval"
                rm -rf "${dataset_path}/tmp/step_07_cluster"
                rm -rf "${dataset_path}/tmp/temp_file"

                # Detect PE/SE from fastp files (tolerate minor mismatches)
                r1_fastp=$(find "$fastp_path" -type f -name '*_1.fastq*' 2>/dev/null | wc -l)
                r2_fastp=$(find "$fastp_path" -type f -name '*_2.fastq*' 2>/dev/null | wc -l)
                if [ "$r1_fastp" -gt 0 ] && [ "$r2_fastp" -gt 0 ] && \
                   [ "$r2_fastp" -ge $(( r1_fastp * 9 / 10 )) ]; then
                    sequence_type="paired"
                else
                    sequence_type="single"
                fi
                echo "Sequence type: ${sequence_type^^}"
                export sequence_type
                original_sequence_type="$sequence_type"

                # Read cached quality status or re-check
                if [[ -f "$quality_cache" ]]; then
                    quality_status=$(cat "$quality_cache")
                    echo "Quality status (cached): $quality_status"
                else
                    echo ">>> Checking quality score diversity..."
                    quality_result=$(python3 "${SCRIPTS}/py_16s.py" check_quality_diversity \
                        --input_dir "$fastp_path" --n_samples 3 --n_reads 1000)
                    quality_status=$(echo "$quality_result" | grep "^QUALITY_STATUS=" | cut -d= -f2)
                    echo "$quality_status" > "$quality_cache"
                    echo "Quality status: $quality_status"
                fi
            else
                # No valid fastp checkpoint — full run from scratch
                if [[ -d "${dataset_path}/tmp" ]]; then
                    echo ">>> No valid fastp checkpoint ($n_fastp_fq files, expected $n_srr or $((n_srr*2))). Cleaning and re-running..."
                    rm -rf "${dataset_path}/tmp"
                    rm -rf "${dataset_path}/ori_fastq"
                fi

                # ── Step A: Download ──
                echo ">>> Downloading SRA data..."
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Downloading..." >&3
                if ! Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}" -b "$dataset_ID"; then
                    echo "Error: Download failed for dataset $dataset_ID" >&2
                    exit 1
                fi

                # Count raw reads before any processing
                Common_CountRawReads "$dataset_path" "$sra_file_name"

                # Per-sample layout detection: classify each sample as PE or SE
                local n_pe=0 n_se=0
                local -a pe_samples=() se_samples=()

                # Collect unique sample prefixes (strip _1/_2/_R1/_R2 and extension)
                for fq in "${ori_fastq_path}/"*.fastq*; do
                    [[ -f "$fq" ]] || continue
                    local bname
                    bname=$(basename "$fq")
                    # Skip R2/_2 files (will be found via R1/_1)
                    [[ "$bname" =~ _2\.fastq || "$bname" =~ _R2 ]] && continue

                    if [[ "$bname" =~ _1\.fastq ]]; then
                        local r2="${fq/_1.fastq/_2.fastq}"
                        if [[ -f "$r2" ]]; then
                            n_pe=$((n_pe + 1))
                            pe_samples+=("$fq")
                        else
                            echo "  Removing unpaired: $bname" >&2
                            rm -f "$fq"
                        fi
                    elif [[ "$bname" =~ _R1 ]]; then
                        local r2="${fq/_R1/_R2}"
                        if [[ -f "$r2" ]]; then
                            n_pe=$((n_pe + 1))
                            pe_samples+=("$fq")
                        else
                            echo "  Removing unpaired: $bname" >&2
                            rm -f "$fq"
                        fi
                    else
                        # Single-end file (no _1/_2 or _R1/_R2 suffix)
                        n_se=$((n_se + 1))
                        se_samples+=("$fq")
                    fi
                done

                # Determine majority layout
                if [[ "$n_pe" -ge "$n_se" ]]; then
                    sequence_type="paired"
                else
                    sequence_type="single"
                fi

                echo "  Layout: ${n_pe} PE + ${n_se} SE samples → majority ${sequence_type^^}"

                # Handle minority samples
                if [[ "$sequence_type" == "single" && "$n_pe" -gt 0 ]]; then
                    echo "  Merging ${n_pe} PE minority samples to SE..."
                    for r1 in "${pe_samples[@]}"; do
                        local r2
                        if [[ "$r1" =~ _1\.fastq ]]; then
                            r2="${r1/_1.fastq/_2.fastq}"
                        else
                            r2="${r1/_R1/_R2}"
                        fi
                        local merged_base="${r1%_[1R]*.*}"
                        local merged_tmp="${merged_base}_merged.fastq"
                        local merged_out="${merged_base}.fastq.gz"
                        # Use vsearch to merge PE reads
                        if vsearch --fastq_mergepairs "$r1" --reverse "$r2" \
                                   --fastqout "$merged_tmp" \
                                   --threads "$cpu" --quiet 2>/dev/null \
                           && [[ -s "$merged_tmp" ]]; then
                            gzip -c "$merged_tmp" > "$merged_out"
                            rm -f "$merged_tmp" "$r1" "$r2"
                            echo "    Merged: $(basename "$r1") → $(basename "$merged_out")"
                        else
                            echo "    Warning: merge failed for $(basename "$r1"), skipping sample" >&2
                            rm -f "$r1" "$r2" "$merged_tmp" "$merged_out"
                        fi
                    done
                elif [[ "$sequence_type" == "paired" && "$n_se" -gt 0 ]]; then
                    echo "  Warning: skipping ${n_se} SE-only samples (incompatible with PE pipeline):"
                    for se_fq in "${se_samples[@]}"; do
                        echo "    Skipped: $(basename "$se_fq")" >&2
                        rm -f "$se_fq"
                    done
                fi

                echo "Sequence type: ${sequence_type^^}"
                export sequence_type
                original_sequence_type="$sequence_type"

                # ── Quality Score Diversity Check (first 3 samples) ──
                echo ">>> Checking quality score diversity..."
                quality_result=$(python3 "${SCRIPTS}/py_16s.py" check_quality_diversity \
                    --input_dir "$ori_fastq_path" --n_samples 3 --n_reads 1000)
                quality_status=$(echo "$quality_result" | grep "^QUALITY_STATUS=" | cut -d= -f2)
                echo "$quality_status" > "$quality_cache"
                echo "Quality status: $quality_status"

                # ── Step B: Remove sequencing adapters with fastp ──
                adapter_removed_path="${dataset_path}/tmp/step_01_adapter_removed"
                mkdir -p "$adapter_removed_path"

                if [[ "$sequence_type" == "paired" ]]; then
                    pe_done=false
                    for r1 in "${ori_fastq_path}/"*_R1*.fastq*; do
                        [[ -f "$r1" ]] || continue
                        r2="${r1/_R1/_R2}"
                        [[ -f "$r2" ]] || continue
                        r1_out=$(basename "$r1"); r1_out="${r1_out/_R1/_1}"
                        r2_out=$(basename "$r2"); r2_out="${r2_out/_R2/_2}"
                        fastp -i "$r1" -I "$r2" \
                              -o "${adapter_removed_path}/${r1_out}" \
                              -O "${adapter_removed_path}/${r2_out}" \
                              --detect_adapter_for_pe \
                              --disable_quality_filtering \
                              --disable_length_filtering \
                              -w "$cpu" \
                              -j "${adapter_removed_path}/fastp.json" \
                              -h "${adapter_removed_path}/fastp.html"
                        pe_done=true
                    done
                    if [[ "$pe_done" == false ]]; then
                        for r1 in "${ori_fastq_path}/"*_1.fastq*; do
                            [[ -f "$r1" ]] || continue
                            r2="${r1/_1.fastq/_2.fastq}"
                            [[ -f "$r2" ]] || continue
                            fastp -i "$r1" -I "$r2" \
                                  -o "${adapter_removed_path}/$(basename "$r1")" \
                                  -O "${adapter_removed_path}/$(basename "$r2")" \
                                  --detect_adapter_for_pe \
                                  --disable_quality_filtering \
                                  --disable_length_filtering \
                                      -w "$cpu" \
                                  -j "${adapter_removed_path}/fastp.json" \
                                  -h "${adapter_removed_path}/fastp.html"
                        done
                    fi
                else
                    for fq in "${ori_fastq_path}/"*.fastq*; do
                        [[ -f "$fq" ]] || continue
                        fastp -i "$fq" \
                              -o "${adapter_removed_path}/$(basename "$fq")" \
                              --disable_quality_filtering \
                              --disable_length_filtering \
                              -w "$cpu" \
                              -j "${adapter_removed_path}/fastp.json" \
                              -h "${adapter_removed_path}/fastp.html"
                    done
                fi

                # ── Step C: Entropy-based primer detection & trimming ──
                mkdir -p "$fastp_path"

                python3 "${SCRIPTS}/entropy_primer_detect.py" \
                    -i "$adapter_removed_path" \
                    -o "$fastp_path" || {
                    echo "  ✗ Entropy primer detection failed"
                    exit 1
                }

                # Delete original and intermediate fastq files to save space
                rm -rf "$ori_fastq_path"
                rm -rf "$adapter_removed_path"
            fi

            # ── From here: same flow regardless of resume or fresh run ──
            _now=$(date +%s); _prep_elapsed=$(( _now - _ds_start ))
            echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [3/3] Processing... (prep: $(( _prep_elapsed / 60 ))m$(( _prep_elapsed % 60 ))s)" >&3
            _proc_start=$_now

            if [[ "$quality_status" == "degraded" ]]; then
                # ── Degraded Quality Branch: VSEARCH pipeline ──
                # Quality scores are binned/unreliable → DADA2 cannot learn error model
                echo ">>> DEGRADED quality scores. Using VSEARCH pipeline..."

                degraded_path="${dataset_path}/tmp/step_02c_degraded_preprocess"

                # ── Resume checkpoint: check if degraded preprocess is already done ──
                n_degraded_fq=0
                if [[ -d "$degraded_path" ]]; then
                    n_degraded_fq=$(find "$degraded_path" -type f -name '*.fastq*' | wc -l)
                fi

                if [[ "$n_degraded_fq" -gt 0 ]] && [[ "$n_degraded_fq" -eq "$n_srr" ]]; then
                    # Degraded preprocess already done — resume from manifest
                    echo ">>> Resuming: found $n_degraded_fq degraded-preprocessed files for $n_srr SRR accessions"
                    # Clean only downstream directories
                    rm -rf "${dataset_path}/tmp/step_06_vsearch_cli"
                    rm -rf "${dataset_path}/tmp/step_07_cluster"
                    rm -rf "${dataset_path}/tmp/temp_file"
                else
                    # Run degraded preprocess from fastp data
                    # Trim 15bp from 5' end + truncate to fixed length (auto-detect)
                    # + N filter (>1 → discard). PE → forward reads only.
                    rm -rf "$degraded_path"
                    mkdir -p "$degraded_path"
                    python3 "${SCRIPTS}/py_16s.py" degraded_quality_preprocess \
                        --input_dir "$fastp_path" \
                        --output_dir "$degraded_path" \
                        --trim_front 15 --truncate_length 0 \
                        --max_n 1 --min_length 50 \
                        --sequence_type "$sequence_type" \
                        --threads "$THREADS_PER_DATASET"
                fi

                # Force SE mode (PE already extracted forward reads only)
                sequence_type="single"
                fastq_path="$degraded_path"
                export fastq_path sequence_type

                # ── Manifest (needed for relabel_reads_for_mapping + import) ──
                Amplicon_Common_MakeManifestFileForQiime2

                # ── Direct vsearch pipeline (bypasses QIIME2 intermediate steps) ──
                # Replaces: ImportFastq → QualityControl → Deduplication → ExportForVsearch
                # Quality/length/N filtering already done in degraded_quality_preprocess
                Amplicon_DegradedQ_DirectDerep
                Amplicon_DegradedQ_VsearchDenoise
                Amplicon_DegradedQ_MapReadsToZotus
                Amplicon_DegradedQ_ImportResults

                # ── QIIME2: filter + cleanup ──
                Amplicon_LS454_FilterLowFreqOTUs
                Amplicon_Common_FinalFilesCleaning
            else
                # ── Normal Branch: DADA2 pipeline ──
                # Quality-filter is skipped for Illumina — DADA2's error model
                # handles quality directly, avoiding redundant double-filtering that
                # causes excessive read loss.
                local _skip_cleanup=false
                fastq_path="$fastp_path"
                export fastq_path
                Common_SanitizeFastq
                Amplicon_Common_MakeManifestFileForQiime2
                Amplicon_Common_ImportFastqToQiime2
                Amplicon_Illumina_DenosingDada2

                # PE fallback: if retention < 50%, retry with forward reads only (SE mode)
                if [[ "$sequence_type" == "paired" ]]; then
                    denoised_table="${dataset_path}/tmp/step_05_denoise/${dataset_ID}-table-denoising.qza"
                    raw_counts_file="${dataset_path}/${dataset_ID}_raw_read_counts.tsv"

                    if [[ -f "$denoised_table" ]]; then
                        final_reads=$(Count_Feature_Table_Reads "$denoised_table")
                        raw_reads=$(Count_Raw_Reads_Total "$raw_counts_file")
                        # raw_reads counts all individual reads (R1+R2), but feature table counts
                        # represent sequences (1 merged pair = 1 count). Use half for fair comparison.
                        raw_pairs=$((raw_reads / 2))

                        if [[ "$raw_pairs" -gt 0 ]]; then
                            retention_pct=$(python3 -c "print(f'{$final_reads / $raw_pairs * 100:.1f}')")
                            echo "  PE retention: ${retention_pct}% (${final_reads}/${raw_pairs})"

                            if python3 -c "import sys; sys.exit(0 if $final_reads / $raw_pairs < 0.5 else 1)"; then
                                echo "  ⚠️ PE retention < 50%. Falling back to SE (forward reads only)..."

                                # Clean up PE denoise outputs
                                rm -rf "${dataset_path}/tmp/step_05_denoise/"
                                rm -rf "${dataset_path}/tmp/temp_file/QualityFilter_vis/"

                                # Create SE manifest from forward reads, preserving PE sample names
                                temp_file_path="${dataset_path}/tmp/temp_file"
                                mkdir -p "$temp_file_path"
                                python3 -c "
import os, glob, sys
manifest_path = sys.argv[1]
fastq_dir = sys.argv[2]
with open(manifest_path, 'w') as f:
    f.write('sample-id\tabsolute-filepath\n')
    for fq in sorted(glob.glob(os.path.join(fastq_dir, '*_1.fastq*'))):
        basename = os.path.basename(fq)
        sample = basename.rsplit('_', 1)[0]
        f.write(f'{sample}\t{fq}\n')
" "${temp_file_path}/${dataset_ID}_manifest.tsv" "$fastq_path"

                                # Re-import as SE
                                sequence_type="single"
                                export sequence_type
                                rm -f "${dataset_path}/tmp/step_03_qza_import/${dataset_ID}.qza"
                                Amplicon_Common_ImportFastqToQiime2

                                # Re-run DADA2 as SE (denoise-single)
                                Amplicon_Illumina_DenosingDada2

                                # Check SE retention
                                denoised_table_se="${dataset_path}/tmp/step_05_denoise/${dataset_ID}-table-denoising.qza"
                                if [[ -f "$denoised_table_se" ]]; then
                                    final_reads_se=$(Count_Feature_Table_Reads "$denoised_table_se")
                                    retention_se=$(python3 -c "print(f'{$final_reads_se / $raw_pairs * 100:.1f}')")
                                    echo "  SE retention: ${retention_se}% (${final_reads_se}/${raw_pairs})"

                                    if python3 -c "import sys; sys.exit(0 if $final_reads_se / $raw_pairs < 0.5 else 1)"; then
                                        echo "  ⚠️ WARNING: SE retention still < 50%. This dataset may have low-quality data."
                                        echo "  Skipping cleanup to preserve intermediate files for debugging."
                                        _skip_cleanup=true
                                        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - LOW_QUALITY - PE: ${retention_pct}%, SE: ${retention_se}%" >> "$low_quality_log"
                                    fi
                                fi
                            fi
                        fi
                    fi
                fi

                if [[ "$_skip_cleanup" == true ]]; then
                    # Copy final outputs but preserve intermediate files for debugging
                    local _denoise="${dataset_path}/tmp/step_05_denoise"
                    if [[ -f "${_denoise}/${dataset_ID}-table-denoising.qza" ]]; then
                        cp "${_denoise}/${dataset_ID}-rep-seqs-denoising.qza" "${dataset_path}/${dataset_ID}-final-rep-seqs.qza"
                        cp "${_denoise}/${dataset_ID}-table-denoising.qza" "${dataset_path}/${dataset_ID}-final-table.qza"
                    fi
                    echo "  [LOW_QUALITY] Intermediate files preserved in: ${dataset_path}/tmp/"
                else
                    Amplicon_Common_FinalFilesCleaning
                fi
            fi

        elif [[ "$platform" == "LS454" ]]; then
            fastp_path="${dataset_path}/tmp/step_02_fastp"
            n_srr=$(wc -l < "${sra_file_name}" | tr -d ' ')

            # ── Resume checkpoint: check if fastp data is intact ──
            n_fastp_fq=0
            if [[ -d "$fastp_path" ]]; then
                n_fastp_fq=$(find "$fastp_path" -type f -name '*.fastq*' | wc -l)
            fi

            if [[ "$n_fastp_fq" -gt 0 ]] && [[ "$n_fastp_fq" -eq "$n_srr" ]]; then
                # 454 is always SE, so n_fastp_fq == n_srr
                echo ">>> Resuming: found $n_fastp_fq fastp files for $n_srr SRR accessions"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Resuming from checkpoint" >&3

                # Clean up downstream directories
                rm -rf "${dataset_path}/tmp/step_02b_adaptive_trim"
                rm -rf "${dataset_path}/tmp/step_03_qza_import"
                rm -rf "${dataset_path}/tmp/step_04_qza_import_QualityFilter"
                rm -rf "${dataset_path}/tmp/step_05_dedupicate"
                rm -rf "${dataset_path}/tmp/step_06_ChimerasRemoval"
                rm -rf "${dataset_path}/tmp/step_07_cluster"
                rm -rf "${dataset_path}/tmp/temp_file"
            else
                if [[ -d "${dataset_path}/tmp" ]]; then
                    echo ">>> No valid fastp checkpoint ($n_fastp_fq files, expected $n_srr). Cleaning and re-running..."
                    rm -rf "${dataset_path}/tmp"
                    rm -rf "${dataset_path}/ori_fastq"
                fi

                echo ">>> Downloading SRA data..."
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Downloading..." >&3
                if ! Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}" -b "$dataset_ID"; then
                    echo "Error: Download failed for dataset $dataset_ID" >&2
                    exit 1
                fi

                Common_CountRawReads "$dataset_path" "$sra_file_name"

                adapter_removed_path="${dataset_path}/tmp/step_01_adapter_removed"
                mkdir -p "$adapter_removed_path"

                for fq in "${ori_fastq_path}/"*.fastq*; do
                    [[ -f "$fq" ]] || continue
                    fastp -i "$fq" \
                          -o "${adapter_removed_path}/$(basename "$fq")" \
                          --disable_quality_filtering \
                          --disable_length_filtering \
                          -w "$cpu" \
                          -j "${adapter_removed_path}/fastp.json" \
                          -h "${adapter_removed_path}/fastp.html"
                done

                mkdir -p "$fastp_path"
                python3 "${SCRIPTS}/entropy_primer_detect.py" \
                    -i "$adapter_removed_path" \
                    -o "$fastp_path" || {
                    echo "  ✗ Entropy primer detection failed"
                    exit 1
                }

                rm -rf "$ori_fastq_path"
                rm -rf "$adapter_removed_path"
            fi

            sequence_type="single"
            export sequence_type
            original_sequence_type="$sequence_type"

            _now=$(date +%s); _prep_elapsed=$(( _now - _ds_start ))
            echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [3/3] Processing... (prep: $(( _prep_elapsed / 60 ))m$(( _prep_elapsed % 60 ))s)" >&3
            _proc_start=$_now

            # ── Step D: Adaptive tail trimming (data-driven N removal) ──
            # Analyses per-position N frequency at 3' end, trims elevated-N
            # tail, then computes P95 of remaining N counts for QC threshold.
            adaptive_trim_path="${dataset_path}/tmp/step_02b_adaptive_trim"
            mkdir -p "$adaptive_trim_path"

            echo ">>> Adaptive tail trimming..."
            trim_result=$(python3 "${SCRIPTS}/py_16s.py" adaptive_tail_trim \
                --input_dir "$fastp_path" \
                --output_dir "$adaptive_trim_path" \
                --max_sample_reads 10000)

            trim_length=$(echo "$trim_result" | grep "^TRIM_LENGTH=" | cut -d= -f2)
            max_ambiguous=$(echo "$trim_result" | grep "^MAX_AMBIGUOUS=" | cut -d= -f2)
            export max_ambiguous

            echo "  Trim length: ${trim_length} bp"
            echo "  Max ambiguous (P95): ${max_ambiguous}"

            # Clean up pre-trim FASTQ
            rm -rf "$fastp_path"

            # ── Step E: QIIME2 Import → Length filter → Dedup → Chimera → OTU 97% → Filter low-freq OTUs ──
            fastq_path="$adaptive_trim_path"
            export fastq_path
            Common_SanitizeFastq
            Amplicon_Common_MakeManifestFileForQiime2
            Amplicon_Common_ImportFastqToQiime2
            Amplicon_LS454_QualityControlForQZA
            Amplicon_LS454_Deduplication
            Amplicon_LS454_ChimerasRemoval
            Amplicon_LS454_ClusterDenovo
            Amplicon_LS454_FilterLowFreqOTUs
            Amplicon_Common_FinalFilesCleaning

        elif [[ "$platform" == "ION_TORRENT" ]]; then
            fastp_path="${dataset_path}/tmp/step_02_fastp"
            n_srr=$(wc -l < "${sra_file_name}" | tr -d ' ')

            # ── Resume checkpoint: check if fastp data is intact ──
            n_fastp_fq=0
            if [[ -d "$fastp_path" ]]; then
                n_fastp_fq=$(find "$fastp_path" -type f -name '*.fastq*' | wc -l)
            fi

            if [[ "$n_fastp_fq" -gt 0 ]] && [[ "$n_fastp_fq" -eq "$n_srr" ]]; then
                echo ">>> Resuming: found $n_fastp_fq fastp files for $n_srr SRR accessions"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Resuming from checkpoint" >&3
                rm -rf "${dataset_path}/tmp/step_03_qza_import"
                rm -rf "${dataset_path}/tmp/step_04_qza_import_QualityFilter"
                rm -rf "${dataset_path}/tmp/step_05_denoise"
                rm -rf "${dataset_path}/tmp/temp_file"
            else
                if [[ -d "${dataset_path}/tmp" ]]; then
                    echo ">>> No valid fastp checkpoint ($n_fastp_fq files, expected $n_srr). Cleaning and re-running..."
                    rm -rf "${dataset_path}/tmp"
                    rm -rf "${dataset_path}/ori_fastq"
                fi

                echo ">>> Downloading SRA data..."
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Downloading..." >&3
                if ! Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}" -b "$dataset_ID"; then
                    echo "Error: Download failed for dataset $dataset_ID" >&2
                    exit 1
                fi

                Common_CountRawReads "$dataset_path" "$sra_file_name"

                adapter_removed_path="${dataset_path}/tmp/step_01_adapter_removed"
                mkdir -p "$adapter_removed_path"

                for fq in "${ori_fastq_path}/"*.fastq*; do
                    [[ -f "$fq" ]] || continue
                    fastp -i "$fq" \
                          -o "${adapter_removed_path}/$(basename "$fq")" \
                          --disable_quality_filtering \
                          --disable_length_filtering \
                          -w "$cpu" \
                          -j "${adapter_removed_path}/fastp.json" \
                          -h "${adapter_removed_path}/fastp.html"
                done

                mkdir -p "$fastp_path"
                python3 "${SCRIPTS}/entropy_primer_detect.py" \
                    -i "$adapter_removed_path" \
                    -o "$fastp_path" || {
                    echo "  ✗ Entropy primer detection failed"
                    exit 1
                }

                rm -rf "$ori_fastq_path"
                rm -rf "$adapter_removed_path"
            fi

            sequence_type="single"
            export sequence_type
            original_sequence_type="$sequence_type"

            _now=$(date +%s); _prep_elapsed=$(( _now - _ds_start ))
            echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [3/3] Processing... (prep: $(( _prep_elapsed / 60 ))m$(( _prep_elapsed % 60 ))s)" >&3
            _proc_start=$_now

            # ── Step D: QIIME2 Import → Quality filter → DADA2 denoise-pyro ──
            # Ion Torrent signal instability in the first ~10bp is handled by
            # DADA2 denoise-pyro --p-trim-left 10 (passed via -s 10).
            # trunc-len is computed automatically from QC visualization.
            fastq_path="$fastp_path"
            export fastq_path
            Common_SanitizeFastq
            Amplicon_Common_MakeManifestFileForQiime2
            Amplicon_Common_ImportFastqToQiime2
            Amplicon_IonTorrent_QualityControlForQZA
            Amplicon_Illumina_DenosingDada2 -s 10
            Amplicon_Common_FinalFilesCleaning

        elif [[ "$platform" == "OXFORD_NANOPORE" ]]; then
            echo "⚠️ OXFORD_NANOPORE platform is not supported for amplicon analysis."
            echo "   Skipping dataset: $dataset_ID"
            echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - SKIPPED - Platform: $platform (not supported)" >> "$failed_log"
            continue

        elif [[ "$platform" == "PACBIO_SMRT" ]]; then
            adapter_removed_path="${dataset_path}/tmp/step_01_adapter_removed"
            n_srr=$(wc -l < "${sra_file_name}" | tr -d ' ')

            # ── Resume checkpoint: check if adapter-removed data is intact ──
            # PacBio skips primer detection (handled by DADA2 denoise-ccs)
            n_adapter_fq=0
            if [[ -d "$adapter_removed_path" ]]; then
                n_adapter_fq=$(find "$adapter_removed_path" -type f -name '*.fastq*' ! -name 'fastp.*' | wc -l)
            fi

            if [[ "$n_adapter_fq" -gt 0 ]] && [[ "$n_adapter_fq" -eq "$n_srr" ]]; then
                echo ">>> Resuming: found $n_adapter_fq adapter-removed files for $n_srr SRR accessions"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Resuming from checkpoint" >&3
                rm -rf "${dataset_path}/tmp/step_03_qza_import"
                rm -rf "${dataset_path}/tmp/step_04_qza_import_QualityFilter"
                rm -rf "${dataset_path}/tmp/step_05_denoise"
                rm -rf "${dataset_path}/tmp/temp_file"
            else
                if [[ -d "${dataset_path}/tmp" ]]; then
                    echo ">>> No valid checkpoint ($n_adapter_fq files, expected $n_srr). Cleaning and re-running..."
                    rm -rf "${dataset_path}/tmp"
                    rm -rf "${dataset_path}/ori_fastq"
                fi

                echo ">>> Downloading SRA data..."
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Downloading..." >&3
                if ! Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}" -b "$dataset_ID"; then
                    echo "Error: Download failed for dataset $dataset_ID" >&2
                    exit 1
                fi

                Common_CountRawReads "$dataset_path" "$sra_file_name"

                mkdir -p "$adapter_removed_path"
                for fq in "${ori_fastq_path}/"*.fastq*; do
                    [[ -f "$fq" ]] || continue
                    fastp -i "$fq" \
                          -o "${adapter_removed_path}/$(basename "$fq")" \
                          --disable_quality_filtering \
                          --disable_length_filtering \
                          -w "$cpu" \
                          -j "${adapter_removed_path}/fastp.json" \
                          -h "${adapter_removed_path}/fastp.html"
                done

                rm -rf "$ori_fastq_path"
            fi

            sequence_type="single"
            export sequence_type
            original_sequence_type="$sequence_type"

            _now=$(date +%s); _prep_elapsed=$(( _now - _ds_start ))
            echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [3/3] Processing... (prep: $(( _prep_elapsed / 60 ))m$(( _prep_elapsed % 60 ))s)" >&3
            _proc_start=$_now

            # ── Step B2: Read length check on first sample ──
            # Sample the first 1000 reads from the first FASTQ file to
            # determine whether these are near-full-length 16S CCS reads.
            first_fq=$(ls "${adapter_removed_path}/"*.fastq* 2>/dev/null | head -n 1)
            if [[ -z "$first_fq" ]]; then
                echo "❌ ERROR: No FASTQ files found after adapter removal"
                exit 1
            fi

            echo ">>> Checking read lengths from first sample: $(basename "$first_fq")"
            long_read_ratio=$(python3 -c "
import sys, gzip, os

fq_path = '${first_fq}'
open_fn = gzip.open if fq_path.endswith('.gz') else open
count = 0
long_count = 0
with open_fn(fq_path, 'rt') as fh:
    while count < 1000:
        header = fh.readline()
        if not header:
            break
        seq = fh.readline().strip()
        fh.readline()  # +
        fh.readline()  # qual
        count += 1
        if len(seq) > 1400:
            long_count += 1
if count == 0:
    print('0.0')
else:
    print(f'{long_count / count:.4f}')
")
            echo "  Reads > 1400 bp ratio: ${long_read_ratio} (from first 1000 reads)"

            # ── Sub-condition A: Full-length 16S CCS reads (majority > 1400bp) ──
            if python3 -c "sys_exit = __import__('sys').exit; sys_exit(0 if float('${long_read_ratio}') > 0.5 else 1)"; then
                echo ">>> Full-length 16S detected (>50% reads > 1400bp). Using 27F/1492R primers."

                # Read primer sequences from reference FASTA files
                DOCS_DIR="${SCRIPT_DIR}/docs"
                primer_front=$(python3 -c "
with open('${DOCS_DIR}/27F.fas') as f:
    lines = f.read().strip().split('\n')
    print(lines[1].strip())
")
                primer_adapter=$(python3 -c "
with open('${DOCS_DIR}/1492R.fas') as f:
    lines = f.read().strip().split('\n')
    print(lines[1].strip())
")
                # Import adapter-removed reads directly into QIIME2
                fastq_path="$adapter_removed_path"
                export fastq_path
                Common_SanitizeFastq
                Amplicon_Common_MakeManifestFileForQiime2
                Amplicon_Common_ImportFastqToQiime2
                Amplicon_Pacbio_QualityControlForQZA

                # DADA2 denoise-ccs with known 27F/1492R primers
                export primer_front
                export primer_adapter
                Amplicon_Pacbio_DenosingDada2
                Amplicon_Pacbio_ExtractReads
                Amplicon_Common_FinalFilesCleaning

            else
                echo "❌ ERROR: PacBio reads are mostly < 1400bp."
                echo "   This pipeline currently only supports full-length 16S PacBio CCS reads."
                echo "   Skipping dataset: $dataset_ID"
                echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - SKIPPED - PacBio reads too short (ratio > 1400bp: ${long_read_ratio})" >> "$failed_log"
                continue
            fi

        else
            echo "❌ Unknown platform: $platform"
            exit 1
        fi

        # Generate summary CSV entry
        echo ">>> Generating summary for $dataset_ID..."
        python "${SCRIPTS}/py_16s.py" append_summary \
            --dataset_id "$dataset_ID" \
            --sra_file "${dataset_path}/${sra_file_name}" \
            --raw_counts "${dataset_path}/${dataset_ID}_raw_read_counts.tsv" \
            --final_table "${dataset_path}/${dataset_ID}-final-table.qza" \
            --output_csv "$summary_csv" \
            --sequence_type "$original_sequence_type"

        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - SUCCESS - Platform: $platform" >> "$success_log"

    )
    local _rc=$?
    local _ds_end=$(date +%s)
    local _total=$(( _ds_end - _ds_start ))
    local _total_fmt="$(( _total / 60 ))m$(( _total % 60 ))s"
    if [[ $_rc -ne 0 ]]; then
        echo "❌ Pipeline failed for $dataset_ID — skipping to next dataset"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - FAILED" >> "$failed_log"
        echo "[$(date '+%H:%M:%S')] [${dataset_ID}] FAILED (${_total_fmt}) - see logs/${dataset_ID}.log" >&3
    else
        echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SUCCESS (${_total_fmt})" >&3
    fi
    }

    set +e
    if [[ "$MAX_PARALLEL" -gt 1 ]]; then
        _process_one_dataset >> "$log_file" 2>&1 &
        running_jobs=$((running_jobs + 1))
    else
        _process_one_dataset 2>&1 | tee "$log_file"
    fi
    set -e
done

# Wait for all remaining background jobs (parallel mode)
if [[ "$MAX_PARALLEL" -gt 1 ]]; then
    echo ">>> Waiting for remaining background datasets to finish..."
    wait
fi

################################################################################
#                          FINAL SUMMARY                                       #
################################################################################

n_success=$(wc -l < "$success_log" 2>/dev/null || echo 0)
n_failed=$(wc -l < "$failed_log" 2>/dev/null || echo 0)
n_skipped=$(wc -l < "$skipped_log" 2>/dev/null || echo 0)
n_low_quality=$(wc -l < "$low_quality_log" 2>/dev/null || echo 0)
_pipeline_end=$(date +%s)
_pipeline_elapsed=$(( _pipeline_end - _pipeline_start ))
_pipeline_min=$(( _pipeline_elapsed / 60 ))
_pipeline_sec=$(( _pipeline_elapsed % 60 ))

echo "========================================="
echo "ALL DONE  (total: ${_pipeline_min}m${_pipeline_sec}s)"
echo "========================================="
echo "  Success:      $n_success"
echo "  Failed:       $n_failed"
echo "  Skipped:      $n_skipped"
echo "  Low quality:  $n_low_quality"
echo "========================================="

if [[ "$n_failed" -gt 0 ]]; then
    echo ""
    echo "Failed datasets:"
    cat "$failed_log"
    echo ""
fi

if [[ "$n_low_quality" -gt 0 ]]; then
    echo ""
    echo "Low quality datasets (retention < 50% after PE+SE fallback):"
    cat "$low_quality_log"
    echo ""
fi

echo "Logs: ${OUTPUT}/logs/"