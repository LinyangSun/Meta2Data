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
MODE=""
LOCAL_MODE=0
LOCAL_PLATFORM=""
PRIMER_FWD=""
PRIMER_REV=""
declare -A LOCAL_SRC    # dataset_id -> source FASTQ directory (local mode)

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
        --mode) MODE="$2"; shift 2 ;;
        --local) LOCAL_MODE=1; shift ;;
        --platform) LOCAL_PLATFORM="$2"; shift 2 ;;
        --primer-fwd) PRIMER_FWD="$2"; shift 2 ;;
        --primer-rev) PRIMER_REV="$2"; shift 2 ;;
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

# Validate denoising mode (mandatory; normally enforced by the bin wrapper,
# re-checked here so run.sh is safe to call directly). asv = DADA2 only,
# otu = vsearch only — no default, no auto-mix.
case "$MODE" in
    asv|otu) ;;
    *) echo "Error: --mode must be 'asv' or 'otu', got '$MODE'"; exit 1 ;;
esac
export MODE

# --local: read FASTQ straight from a folder (no download, no NCBI detection).
# Requires an explicit --platform; --primer-fwd/--primer-rev are optional (without
# them the entropy auto-detector is used, exactly as in download mode).
export LOCAL_MODE LOCAL_PLATFORM PRIMER_FWD PRIMER_REV
if [[ "$LOCAL_MODE" == "1" ]]; then
    case "$LOCAL_PLATFORM" in
        ILLUMINA|LS454|ION_TORRENT|PACBIO_SMRT|OXFORD_NANOPORE) ;;
        *) echo "Error: --local requires --platform one of ILLUMINA|LS454|ION_TORRENT|PACBIO_SMRT|OXFORD_NANOPORE (got '$LOCAL_PLATFORM')"; exit 1 ;;
    esac
fi

# Strip trailing slashes from paths
METADATA="${METADATA%/}"
OUTPUT="${OUTPUT%/}"

# Validate input
if [[ -z "$METADATA" ]]; then
    echo "Error: --metadata is required"
    exit 1
fi

if [[ "$LOCAL_MODE" == "1" ]]; then
    if [[ ! -d "$METADATA" ]]; then
        echo "Error: --local input folder not found: '$METADATA'. Please check the path."
        exit 1
    fi
elif [[ ! -f "$METADATA" ]]; then
    echo "Error: Metadata file not found: '$METADATA'. Please check the path."
    exit 1
fi

if [[ -z "$OUTPUT" ]]; then
    if [[ "$LOCAL_MODE" == "1" ]]; then
        echo "Error: --local requires an explicit -o/--output directory (the input folder must not double as the output)"; exit 1
    fi
    OUTPUT=$(dirname "$METADATA")
fi

# Compute per-dataset thread count: total threads ÷ max parallel datasets
THREADS_PER_DATASET=$(( THREADS / MAX_PARALLEL ))
if [[ "$THREADS_PER_DATASET" -lt 1 ]]; then
    THREADS_PER_DATASET=1
fi
export THREADS_PER_DATASET
export cpu=$THREADS_PER_DATASET

if [[ -f "${SCRIPTS}/AmpliconFunction.sh" ]]; then
    source "${SCRIPTS}/AmpliconFunction.sh"
else
    echo "Error: AmpliconFunction.sh not found at ${SCRIPTS}/AmpliconFunction.sh"
    exit 1
fi

_fastp_se_adapter_remove() {
    local ori_fastq_path="$1"
    local adapter_removed_path="$2"
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
}

# Emit the "[3/3] Processing..." milestone with prep-phase elapsed time.
# Args: $1 = dataset start epoch (_ds_start), $2 = dataset id (dataset_ID).
# Writes to fd 3 (the saved console stdout), inherited from the caller.
_emit_prep_done() {
    local _ds_start="$1"
    local dataset_ID="$2"
    local _now=$(date +%s); local _prep_elapsed=$(( _now - _ds_start ))
    echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [3/3] Processing... (prep: $(( _prep_elapsed / 60 ))m$(( _prep_elapsed % 60 ))s)" >&3
}

cd "$OUTPUT" || exit 1
mkdir -p "${OUTPUT}/logs"

_stage_local_reads() {
    # _stage_local_reads <dataset_path> <dataset_id> — symlink local FASTQ into
    # ori_fastq (read-only on the user's originals; the pipeline only reads them
    # and removes the symlinks afterwards). Source dir comes from LOCAL_SRC[id].
    local dest="${1%/}/ori_fastq"
    local src="${LOCAL_SRC[$2]}"
    [[ -n "$src" && -d "$src" ]] || { echo "Error: no local source dir for $2" >&2; return 1; }
    mkdir -p "$dest"
    local n=0 fq dname
    shopt -s nullglob
    for fq in "$src"/*.fastq "$src"/*.fastq.gz "$src"/*.fq "$src"/*.fq.gz; do
        dname=$(basename "$fq")
        # Normalize .fq -> .fastq so the rest of the pipeline (which only matches
        # *.fastq*) processes the file.
        case "$dname" in
            *.fq.gz) dname="${dname%.fq.gz}.fastq.gz" ;;
            *.fq)    dname="${dname%.fq}.fastq" ;;
        esac
        ln -sf "$(readlink -f "$fq")" "$dest/$dname"
        n=$((n + 1))
    done
    shopt -u nullglob
    [[ "$n" -gt 0 ]] || { echo "Error: no FASTQ files in $src" >&2; return 1; }
    echo "  Staged $n local FASTQ file(s) from $src"
}

_obtain_reads() {
    # _obtain_reads <dataset_path> <sra_file_name> <dataset_id>
    # Download (normal mode) or symlink local files (--local).
    if [[ "${LOCAL_MODE:-0}" == "1" ]]; then
        _stage_local_reads "$1" "$3"
    else
        Common_SRADownloadToFastq_MultiSource -d "$1" -a "$2" -b "$3"
    fi
}

_trim_primers() {
    # _trim_primers <in_dir> <out_dir>
    # Explicit primers (cutadapt) when --primer-fwd is given, else the entropy
    # auto-detector (identical to download-mode behaviour).
    local in_dir="$1" out_dir="$2"
    mkdir -p "$out_dir"
    if [[ -n "${PRIMER_FWD:-}" ]]; then
        echo "  Trimming primers with cutadapt (fwd=${PRIMER_FWD}${PRIMER_REV:+, rev=$PRIMER_REV})..."
        local fq base r2
        shopt -s nullglob
        for fq in "$in_dir"/*.fastq "$in_dir"/*.fastq.gz "$in_dir"/*.fq "$in_dir"/*.fq.gz; do
            base=$(basename "$fq")
            [[ "$base" =~ _2\.fastq || "$base" =~ _R2 ]] && continue   # handled via its R1
            if [[ "$base" =~ _1\.fastq || "$base" =~ _R1 ]]; then
                if [[ "$base" =~ _1\.fastq ]]; then r2="${fq/_1./_2.}"; else r2="${fq/_R1/_R2}"; fi
                if [[ -n "${PRIMER_REV:-}" && -f "$r2" ]]; then
                    cutadapt -g "$PRIMER_FWD" -G "$PRIMER_REV" \
                        -o "$out_dir/$base" -p "$out_dir/$(basename "$r2")" "$fq" "$r2" >/dev/null
                else
                    cutadapt -g "$PRIMER_FWD" -o "$out_dir/$base" "$fq" >/dev/null
                fi
            else
                cutadapt -g "$PRIMER_FWD" -o "$out_dir/$base" "$fq" >/dev/null
            fi
        done
        shopt -u nullglob
    else
        python3 "${SCRIPTS}/entropy_primer_detect.py" -i "$in_dir" -o "$out_dir"
    fi
}

_local_register_dataset() {
    # _local_register_dataset <id> <src_dir> — register one local dataset:
    # create its dir, map its source, and write a synthetic <id>_sra.txt
    # (Run<TAB>SampleName per unique sample prefix) so Common_CountRawReads and
    # append_summary work exactly as in download mode.
    local id="$1" src="$2"
    local dpath="${OUTPUT}/${id}"
    mkdir -p "$dpath"
    LOCAL_SRC["$id"]="$src"
    Dataset_ID_sets+=("$id")
    echo "$id" >> "${OUTPUT}/datasets_ID.txt"
    local sra="${dpath}/${id}_sra.txt"; : > "$sra"
    local fq base p
    declare -A _seen=()
    shopt -s nullglob
    for fq in "$src"/*.fastq "$src"/*.fastq.gz "$src"/*.fq "$src"/*.fq.gz; do
        base=$(basename "$fq")
        p="${base%.gz}"; p="${p%.fastq}"; p="${p%.fq}"
        p="${p%_R1}"; p="${p%_R2}"; p="${p%_1}"; p="${p%_2}"
        [[ -n "${_seen[$p]:-}" ]] && continue
        _seen["$p"]=1
        printf '%s\t%s\n' "$p" "$p" >> "$sra"
    done
    shopt -u nullglob
    # Warn on prefix collisions: Common_CountRawReads globs "<prefix>*.fastq*", so a
    # sample name that is a prefix of another would double-count raw reads.
    local -a _pf; mapfile -t _pf < <(cut -f1 "$sra")
    local _a _b
    for _a in "${_pf[@]}"; do
        for _b in "${_pf[@]}"; do
            [[ "$_a" != "$_b" && "$_b" == "$_a"* ]] && \
                echo "  ⚠️  Warning: sample '$_a' is a prefix of '$_b' — raw-read counts may be inflated; rename to avoid overlap (e.g. sample01/sample10)." >&2
        done
    done
    echo "  Dataset '$id': $(wc -l < "$sra") sample(s) from $src"
}

_local_discover_datasets() {
    # --local takes ONE folder = ONE dataset; every FASTQ directly inside it is a
    # sample of that dataset (dataset id = folder name). No sub-folder recursion:
    # since --platform is a single value, one --local run handles exactly one
    # dataset / one platform. Mixed-platform data must be run separately, one
    # folder per run.
    local input="${METADATA%/}"
    Dataset_ID_sets=()
    : > "${OUTPUT}/datasets_ID.txt"
    if compgen -G "$input/"*.fastq* >/dev/null 2>&1 || compgen -G "$input/"*.fq* >/dev/null 2>&1; then
        echo "  Local mode: single dataset '$(basename "$input")' from $input"
        _local_register_dataset "$(basename "$input")" "$input"
    else
        echo "Error: no FASTQ files found directly in '$input'"; exit 1
    fi
    [[ ${#Dataset_ID_sets[@]} -gt 0 ]] || { echo "Error: no local dataset discovered"; exit 1; }
}

################################################################################
#                          PHASE 1: DATASET PREPARATION                        #
################################################################################

echo "========================================="
echo "PHASE 1: Dataset Preparation"
echo "Started: $(date)"
echo "========================================="

if [[ "$LOCAL_MODE" == "1" ]]; then
    _local_discover_datasets
else
    if ! python "${SCRIPTS}/py_16s.py" GenerateDatasetsIDsFile --FilePath "$METADATA" --Bioproject "$COL_BIOPROJECT" --OutputDir "$OUTPUT"; then
        echo "❌ ERROR: Failed to generate dataset IDs, please check your metadata file and column names for BioProject."
        exit 1
    fi

    mapfile -t Dataset_ID_sets < <(awk '{print $1}' "${OUTPUT}/datasets_ID.txt")

    if [ ${#Dataset_ID_sets[@]} -eq 0 ]; then
        echo "❌ ERROR: No datasets found, please check your metadata file and column names for BioProject."
        exit 1
    fi

    if ! python "${SCRIPTS}/py_16s.py" GenerateSRAsFile --FilePath "$METADATA" --Bioproject "$COL_BIOPROJECT" --SRA_Number "$COL_SRA" --OutputDir "$OUTPUT"; then
        echo "❌ ERROR: Failed to generate SRA file lists, please check your metadata file and column names for SRA."
        exit 1
    fi
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
if [[ "$LOCAL_MODE" != "1" ]]; then
_pairs_file="${OUTPUT}/.platform_query_pairs.txt"
: > "$PLATFORM_CACHE_FILE"
: > "$_pairs_file"

# Collect dataset_id<TAB>first_srr[<TAB>bioproject_id] for datasets that still
# lack a platform (i.e. not already processed).
for _ds_id in "${Dataset_ID_sets[@]}"; do
    _ds_path="${OUTPUT}/${_ds_id}"

    # Skip already processed (mode-specific: an asv run does not block a later otu run)
    [[ -f "${_ds_path}/${_ds_id}-${MODE}-final-rep-seqs.qza" ]] && continue

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

# Single Python call for the remainder: batch Entrez for NCBI, serial for CNCB
if [[ -s "$_pairs_file" ]]; then
    _n_queries=$(wc -l < "$_pairs_file" | tr -d ' ')
    echo "  Querying NCBI/CNCB for ${_n_queries} dataset(s)..."
    python "${SCRIPTS}/py_16s.py" batch_get_sequencing_platforms --pairs_file "$_pairs_file" >> "$PLATFORM_CACHE_FILE"
    echo "  Platform cache now has $(wc -l < "$PLATFORM_CACHE_FILE" | tr -d ' ') entries"
fi

rm -f "$_pairs_file"
fi   # end Phase 1.5 (skipped in --local mode)
echo ""

################################################################################
#                   PHASE 2: INDIVIDUAL DATASET PROCESSING                     #
################################################################################

echo "========================================="
echo "PHASE 2: Individual Dataset Processing"
echo "Started: $(date)"
echo "========================================="

# Unified status log (append-only audit trail; NEVER truncated). One line per
# event:  <timestamp>\t<STATUS>\t<dataset_id>\t<detail>
# STATUS in SUCCESS | FAILED | SKIPPED | LOW_QUALITY
RUN_LOG="${OUTPUT}/datasets.log"
summary_csv="${OUTPUT}/summary.csv"

_log_status() {
    # _log_status <STATUS> <dataset_id> <detail>
    printf '%s\t%s\t%s\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$1" "$2" "$3" >> "$RUN_LOG"
}

# Dated run header, then record the current line count so the end-of-run tally
# counts only THIS run's events (the log accumulates across runs).
printf '# === RUN %s | mode=%s | metadata=%s ===\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$MODE" "$METADATA" >> "$RUN_LOG"
_log_start=$(wc -l < "$RUN_LOG" 2>/dev/null || echo 0)

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

    # Check if already processed (mode-specific name)
    if [ -f "${dataset_path}/${dataset_ID}-${MODE}-final-rep-seqs.qza" ]; then
        echo "✓ Already processed. Skipping."
        _log_status SKIPPED "$dataset_ID" "ALREADY_DONE"
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

        # 1. Platform Detection — --local uses --platform; otherwise the
        #    pre-detected cache, with an API fallback on cache miss.
        if [[ "$LOCAL_MODE" == "1" ]]; then
            platform="$LOCAL_PLATFORM"
        else
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
        fi

        echo "Detected platform: $platform"
        # Persist platform per dataset so the Phase 3 summary stays correct on
        # re-runs (the shared .platform_cache.txt is rebuilt each run and skips
        # already-completed datasets).
        echo "$platform" > "${dataset_path}/${dataset_ID}_platform.txt"
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
                    # Migrate pre-A3b cache value ("degraded" -> "degraded_binned").
                    [[ "$quality_status" == "degraded" ]] && quality_status="degraded_binned"
                    # Re-test if the cached token is stale/invalid (the only valid
                    # tokens are normal | degraded_binned). Prevents a stale value
                    # from mis-routing the asv skip / otu preprocess branches.
                    if [[ "$quality_status" != "normal" && "$quality_status" != "degraded_binned" ]]; then
                        echo ">>> Quality cache stale/invalid ('$quality_status'); re-testing..."
                        quality_result=$(python3 "${SCRIPTS}/py_16s.py" check_quality_diversity \
                            --input_dir "$fastp_path" --n_samples 3 --n_reads 1000)
                        quality_status=$(echo "$quality_result" | grep "^QUALITY_STATUS=" | cut -d= -f2)
                    fi
                    echo "$quality_status" > "$quality_cache"
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
                if [[ "$MODE" == "asv" && "$LOCAL_MODE" != "1" ]]; then
                    # ── Early quality probe (asv only) ──
                    # asv skips degraded/binned-quality data. Quality (Q-score
                    # binning) is a dataset-wide property, so probe the first 2
                    # samples first; if binned, skip BEFORE downloading the rest
                    # (avoids fetching a large dataset just to drop it). Probe
                    # goes to a separate dir; on proceed it is folded into
                    # ori_fastq AFTER the rest downloads (so a NCBI-fallback wipe
                    # in the rest download can't clobber it).
                    probe_base="${dataset_path}/tmp/quality_probe"
                    rm -rf "$probe_base"; mkdir -p "$probe_base"
                    head -n 2 "${dataset_path}/${sra_file_name}" > "${probe_base}/probe_sra.txt"
                    if ! Common_SRADownloadToFastq_MultiSource -d "$probe_base" -a "probe_sra.txt" -b "$dataset_ID"; then
                        echo "Error: probe download failed for dataset $dataset_ID" >&2
                        exit 1
                    fi
                    echo ">>> Checking quality score diversity (early probe, first 2 samples)..."
                    quality_result=$(python3 "${SCRIPTS}/py_16s.py" check_quality_diversity \
                        --input_dir "${probe_base}/ori_fastq" --n_samples 2 --n_reads 1000)
                    quality_status=$(echo "$quality_result" | grep "^QUALITY_STATUS=" | cut -d= -f2)
                    echo "$quality_status" > "$quality_cache"
                    echo "Quality status (early probe): $quality_status"

                    if [[ "$quality_status" == "degraded_binned" ]]; then
                        echo ">>> SKIP: degraded/binned quality is incompatible with --asv (DADA2)."
                        echo ">>>       Skipped before downloading the rest of the dataset."
                        _log_status SKIPPED "$dataset_ID" "asv mode, quality=degraded_binned (early probe)"
                        echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SKIPPED (asv + degraded_binned, early)" >&3
                        rm -rf "$probe_base"
                        exit 0
                    fi

                    # Normal quality → download the remaining samples (lines 3..N)
                    tail -n +3 "${dataset_path}/${sra_file_name}" > "${dataset_path}/.rest_sra.txt"
                    if [[ -s "${dataset_path}/.rest_sra.txt" ]]; then
                        if ! Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a ".rest_sra.txt" -b "$dataset_ID"; then
                            echo "Error: Download failed for dataset $dataset_ID" >&2
                            exit 1
                        fi
                    fi
                    # Fold the probe samples into ori_fastq (after the rest download)
                    mkdir -p "$ori_fastq_path"
                    if compgen -G "${probe_base}/ori_fastq/*" > /dev/null 2>&1; then
                        mv "${probe_base}/ori_fastq/"* "$ori_fastq_path/"
                    fi
                    rm -rf "$probe_base" "${dataset_path}/.rest_sra.txt"
                else
                    # otu mode: no early skip (quality only selects maxee vs
                    # truncation later), so download everything up front.
                    if ! _obtain_reads "$dataset_path" "${sra_file_name}" "$dataset_ID"; then
                        echo "Error: Download failed for dataset $dataset_ID" >&2
                        exit 1
                    fi
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
                # asv already determined quality from the early probe above; only
                # otu needs it here (to pick maxee vs truncation preprocess).
                if [[ "$MODE" != "asv" || "$LOCAL_MODE" == "1" ]]; then
                    echo ">>> Checking quality score diversity..."
                    quality_result=$(python3 "${SCRIPTS}/py_16s.py" check_quality_diversity \
                        --input_dir "$ori_fastq_path" --n_samples 3 --n_reads 1000)
                    quality_status=$(echo "$quality_result" | grep "^QUALITY_STATUS=" | cut -d= -f2)
                    echo "$quality_status" > "$quality_cache"
                    echo "Quality status: $quality_status"
                fi

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
                    _fastp_se_adapter_remove "$ori_fastq_path" "$adapter_removed_path"
                fi

                # ── Step C: Entropy-based primer detection & trimming ──
                mkdir -p "$fastp_path"

                _trim_primers "$adapter_removed_path" "$fastp_path" || {
                    echo "  ✗ Entropy primer detection failed"
                    exit 1
                }

                # Delete original and intermediate fastq files to save space
                rm -rf "$ori_fastq_path"
                rm -rf "$adapter_removed_path"
            fi

            # ── From here: same flow regardless of resume or fresh run ──
            _emit_prep_done "$_ds_start" "$dataset_ID"

            if [[ "$MODE" == "otu" ]]; then
                # ── OTU mode: merge+maxee (normal) / forward-only truncation
                #    (binned) preprocess, then the shared pooled vsearch chain. ──
                Amplicon_Illumina_OTU_Preprocess
                OTU_STRAND="plus"; export OTU_STRAND   # Illumina short reads: plus only
                Amplicon_OTU_RunPooledChain
            elif [[ "$quality_status" == "degraded_binned" ]]; then
                # ── asv mode: degraded/binned quality is incompatible with DADA2 ──
                # DADA2's error model needs reliable per-base quality scores; binned
                # quality (NovaSeq/HiSeq Q-score compression, re-uploaded data) breaks
                # it. Per the OTU/ASV split contract, --asv NEVER reroutes to vsearch —
                # such datasets are skipped here and belong to --otu instead.
                # (The previous vsearch orchestration for this case is preserved in
                #  AmpliconPIP_OTU_ASV_changelist.md Appendix A for the OTU back-end.)
                echo ">>> SKIP: degraded/binned quality is incompatible with --asv (DADA2)."
                echo ">>>       Use --otu for this dataset (vsearch handles binned quality)."
                _log_status SKIPPED "$dataset_ID" "asv mode, quality=degraded_binned"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SKIPPED (asv + degraded_binned)" >&3
                exit 0
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
                                        _log_status LOW_QUALITY "$dataset_ID" "PE: ${retention_pct}%, SE: ${retention_se}%"
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
                        cp "${_denoise}/${dataset_ID}-rep-seqs-denoising.qza" "${dataset_path}/${dataset_ID}-${MODE}-final-rep-seqs.qza"
                        cp "${_denoise}/${dataset_ID}-table-denoising.qza" "${dataset_path}/${dataset_ID}-${MODE}-final-table.qza"
                    fi
                    echo "  [LOW_QUALITY] Intermediate files preserved in: ${dataset_path}/tmp/"
                else
                    Amplicon_Common_FinalFilesCleaning
                fi
            fi

        elif [[ "$platform" == "LS454" ]]; then
            if [[ "$MODE" == "asv" ]]; then
                # 454 has no DADA2 method → cannot produce ASVs. Skip (belongs to --otu).
                echo ">>> SKIP: LS454 (454) has no DADA2 method → not supported in --asv. Use --otu."
                _log_status SKIPPED "$dataset_ID" "asv mode, platform=LS454"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SKIPPED (asv: 454 unsupported)" >&3
                exit 0
            fi
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
                if ! _obtain_reads "$dataset_path" "${sra_file_name}" "$dataset_ID"; then
                    echo "Error: Download failed for dataset $dataset_ID" >&2
                    exit 1
                fi

                Common_CountRawReads "$dataset_path" "$sra_file_name"

                adapter_removed_path="${dataset_path}/tmp/step_01_adapter_removed"
                mkdir -p "$adapter_removed_path"

                _fastp_se_adapter_remove "$ori_fastq_path" "$adapter_removed_path"

                mkdir -p "$fastp_path"
                _trim_primers "$adapter_removed_path" "$fastp_path" || {
                    echo "  ✗ Entropy primer detection failed"
                    exit 1
                }

                rm -rf "$ori_fastq_path"
                rm -rf "$adapter_removed_path"
            fi

            sequence_type="single"
            export sequence_type
            original_sequence_type="$sequence_type"

            _emit_prep_done "$_ds_start" "$dataset_ID"

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

            # ── OTU back-end (unified pooled vsearch chain) ──
            # 454 is OTU-only (no DADA2 method). §4.x (user-accepted) behaviour
            # change: the old QIIME2 dedup→chimera→cluster path is replaced by the
            # shared UNOISE3 → cluster_fast 97% → map-back chain. Abundance now
            # comes from read map-back (not cluster size), with added UNOISE3
            # denoising. NOTE for review: the old LS454_QualityControlForQZA
            # q-score/length filter is dropped; adaptive_tail_trim handles the
            # 3' N-tail but there is no explicit quality (maxee) filter for 454.
            fastq_path="$adaptive_trim_path"
            export fastq_path
            sequence_type="single"; export sequence_type
            OTU_STRAND="plus"; export OTU_STRAND
            Amplicon_OTU_RunPooledChain

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
                if ! _obtain_reads "$dataset_path" "${sra_file_name}" "$dataset_ID"; then
                    echo "Error: Download failed for dataset $dataset_ID" >&2
                    exit 1
                fi

                Common_CountRawReads "$dataset_path" "$sra_file_name"

                adapter_removed_path="${dataset_path}/tmp/step_01_adapter_removed"
                mkdir -p "$adapter_removed_path"

                _fastp_se_adapter_remove "$ori_fastq_path" "$adapter_removed_path"

                mkdir -p "$fastp_path"
                _trim_primers "$adapter_removed_path" "$fastp_path" || {
                    echo "  ✗ Entropy primer detection failed"
                    exit 1
                }

                rm -rf "$ori_fastq_path"
                rm -rf "$adapter_removed_path"
            fi

            sequence_type="single"
            export sequence_type
            original_sequence_type="$sequence_type"

            _emit_prep_done "$_ds_start" "$dataset_ID"

            if [[ "$MODE" == "otu" ]]; then
                # ── OTU: strip 5' 10bp + maxee → shared pooled vsearch chain ──
                Amplicon_IonTorrent_OTU_Preprocess
                OTU_STRAND="plus"; export OTU_STRAND
                Amplicon_OTU_RunPooledChain
            else
                # ── ASV: QIIME2 Import → Quality filter → DADA2 denoise-pyro ──
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
            fi

        elif [[ "$platform" == "OXFORD_NANOPORE" ]]; then
            if [[ "$MODE" == "asv" ]]; then
                # Single-base ASV resolution is conceptually invalid for ONT
                # (~5-10% error + indels): true sequences explode into a cloud of
                # spurious variants. Skip in --asv; ONT belongs to --otu.
                echo ">>> SKIP: ONT single-base ASV is invalid (5-10% error+indels) → not supported in --asv. Use --otu."
                _log_status SKIPPED "$dataset_ID" "asv mode, platform=OXFORD_NANOPORE"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SKIPPED (asv: ONT unsupported)" >&3
                exit 0
            fi
            # ── Oxford Nanopore: long-read amplicon (faithful ONT-AmpSeq port) ──
            # Lazy dependency check: only ONT datasets need these long-read
            # tools, so non-ONT runs are never gated on them.
            ont_missing=()
            for _t in chopper minimap2 racon vsearch; do
                command -v "$_t" >/dev/null 2>&1 || ont_missing+=("$_t")
            done
            if [[ ${#ont_missing[@]} -gt 0 ]]; then
                # Report to stderr only; the outer subshell handler owns the
                # single canonical FAILED log entry (matches other branches).
                echo "❌ ERROR: ONT processing requires missing tools: ${ont_missing[*]}" >&2
                echo "   Install via: conda install -c bioconda ${ont_missing[*]}  (or module load)" >&2
                exit 1
            fi

            primer_trim_path="${dataset_path}/tmp/step_02_primer"
            chopper_path="${dataset_path}/tmp/step_03_chopper"

            # ── Resume checkpoint: chopper output complete? ──
            if [[ -f "${chopper_path}/.chopper_done" ]]; then
                echo ">>> Resuming: found completed chopper-filtered reads"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Resuming from checkpoint" >&3
                # Clean only downstream directories
                rm -rf "${dataset_path}/tmp/step_06_ont"
                rm -rf "${dataset_path}/tmp/step_07_cluster"
                rm -rf "${dataset_path}/tmp/temp_file"
            else
                if [[ -d "${dataset_path}/tmp" ]]; then
                    echo ">>> No valid chopper checkpoint. Cleaning and re-running..."
                    rm -rf "${dataset_path}/tmp"
                    rm -rf "${dataset_path}/ori_fastq"
                fi

                # ── Step A: Download ──
                echo ">>> Downloading SRA data..."
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Downloading..." >&3
                if ! _obtain_reads "$dataset_path" "${sra_file_name}" "$dataset_ID"; then
                    echo "Error: Download failed for dataset $dataset_ID" >&2
                    exit 1
                fi

                Common_CountRawReads "$dataset_path" "$sra_file_name"

                # ── Step B: Remove sequencing adapters with fastp (SE) ──
                adapter_removed_path="${dataset_path}/tmp/step_01_adapter_removed"
                mkdir -p "$adapter_removed_path"
                _fastp_se_adapter_remove "$ori_fastq_path" "$adapter_removed_path"

                # ── Step C: Entropy-based primer detection & trimming ──
                mkdir -p "$primer_trim_path"
                _trim_primers "$adapter_removed_path" "$primer_trim_path" || {
                    echo "  ✗ Entropy primer detection failed"
                    exit 1
                }

                rm -rf "$ori_fastq_path"
                rm -rf "$adapter_removed_path"

                # ── Step D: chopper length-window (auto peak) + quality filter ──
                fastq_path="$primer_trim_path"
                export fastq_path
                Amplicon_ONT_ChopperFilter

                rm -rf "$primer_trim_path"
            fi

            sequence_type="single"
            export sequence_type
            original_sequence_type="$sequence_type"
            export ONT_FASTQ_DIR="$chopper_path"

            _emit_prep_done "$_ds_start" "$dataset_ID"

            # ── Manifest from chopper-filtered reads (sample-id = run accession) ──
            fastq_path="$chopper_path"
            export fastq_path
            Amplicon_Common_MakeManifestFileForQiime2

            # ── ONT-AmpSeq core: per-sample UNOISE3 -> racon polish -> cluster 97% ──
            Amplicon_ONT_ClusterPerSample
            Amplicon_ONT_PolishRacon
            Amplicon_ONT_RelabelMerge
            Amplicon_ONT_ClusterID
            # ── Abundance via read mapping (DegradedQ-style) + QIIME2 import ──
            Amplicon_ONT_MapReadsToOTUs
            # ── Filter low-freq OTUs + finalize (reused) ──
            Amplicon_LS454_FilterLowFreqOTUs
            Amplicon_Common_FinalFilesCleaning

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
                if ! _obtain_reads "$dataset_path" "${sra_file_name}" "$dataset_ID"; then
                    echo "Error: Download failed for dataset $dataset_ID" >&2
                    exit 1
                fi

                Common_CountRawReads "$dataset_path" "$sra_file_name"

                mkdir -p "$adapter_removed_path"
                _fastp_se_adapter_remove "$ori_fastq_path" "$adapter_removed_path"

                rm -rf "$ori_fastq_path"
            fi

            sequence_type="single"
            export sequence_type
            original_sequence_type="$sequence_type"

            _emit_prep_done "$_ds_start" "$dataset_ID"

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
                echo ">>> Full-length 16S detected (>50% reads > 1400bp)."

                if [[ "$MODE" == "otu" ]]; then
                # ── OTU: length-window + maxee_rate → shared pooled chain (strand both) ──
                Amplicon_Pacbio_OTU_Preprocess
                OTU_STRAND="both"; export OTU_STRAND
                Amplicon_OTU_RunPooledChain
                else
                # ── ASV: DADA2 denoise-ccs with known 27F/1492R primers ──
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
                fi

            else
                echo ">>> SKIP: PacBio reads are mostly < 1400bp (full-length 16S CCS only)."
                _log_status SKIPPED "$dataset_ID" "PacBio reads too short (ratio >1400bp: ${long_read_ratio})"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SKIPPED (PacBio reads too short)" >&3
                exit 0
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
            --final_table "${dataset_path}/${dataset_ID}-${MODE}-final-table.qza" \
            --output_csv "$summary_csv" \
            --sequence_type "$original_sequence_type"

        _log_status SUCCESS "$dataset_ID" "Platform: $platform"

    )
    local _rc=$?
    local _ds_end=$(date +%s)
    local _total=$(( _ds_end - _ds_start ))
    local _total_fmt="$(( _total / 60 ))m$(( _total % 60 ))s"
    if [[ $_rc -eq 99 ]]; then
        # Untrustworthy data — already logged (LOW_QUALITY) inside subshell
        echo "⚠ Skipped $dataset_ID — untrustworthy data (see ${dataset_ID}-UNTRUSTABLE.txt)"
        echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SKIPPED-UNTRUSTABLE (${_total_fmt})" >&3
    elif [[ $_rc -ne 0 ]]; then
        echo "❌ Pipeline failed for $dataset_ID — skipping to next dataset"
        _log_status FAILED "$dataset_ID" "see logs/${dataset_ID}.log"
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
#                   PHASE 3: PER-DATASET SUMMARY                               #
################################################################################
# One row per successful dataset: platform + quality status + amplified V-region
# (rep-seqs aligned to the E. coli 16S reference). Upserted, so re-runs never drop
# previously-summarised datasets. Non-fatal: failure here never fails the pipeline.

echo "========================================="
echo "PHASE 3: Per-dataset summary (platform + quality + amplified region)"
echo "Started: $(date)"
echo "========================================="
python "${SCRIPTS}/py_16s.py" build_per_dataset_summary \
    --output_dir "$OUTPUT" \
    --mode "$MODE" \
    --ecoli_ref "${SCRIPT_DIR}/docs/ecoli_16S_J01859.fasta" \
    --output_csv "$summary_csv" \
    --threads "$THREADS" || echo "  Warning: per-dataset summary generation failed (non-fatal)"
echo ""

################################################################################
#                          FINAL SUMMARY                                       #
################################################################################

n_success=$(awk -F'\t' -v s="$_log_start" 'NR>s && $2=="SUCCESS"' "$RUN_LOG" 2>/dev/null | wc -l)
n_failed=$(awk -F'\t' -v s="$_log_start" 'NR>s && $2=="FAILED"' "$RUN_LOG" 2>/dev/null | wc -l)
n_skipped=$(awk -F'\t' -v s="$_log_start" 'NR>s && $2=="SKIPPED"' "$RUN_LOG" 2>/dev/null | wc -l)
n_low_quality=$(awk -F'\t' -v s="$_log_start" 'NR>s && $2=="LOW_QUALITY"' "$RUN_LOG" 2>/dev/null | wc -l)
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
    awk -F'\t' -v s="$_log_start" 'NR>s && $2=="FAILED"' "$RUN_LOG"
    echo ""
fi

if [[ "$n_low_quality" -gt 0 ]]; then
    echo ""
    echo "Low quality datasets (low retention or untrustworthy single-sample data):"
    awk -F'\t' -v s="$_log_start" 'NR>s && $2=="LOW_QUALITY"' "$RUN_LOG"
    echo ""
fi

echo "Logs: ${OUTPUT}/logs/"
echo "Status log: $RUN_LOG"