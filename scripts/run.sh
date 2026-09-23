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
# re-checked here so run.sh is safe to call directly). dada2 = DADA2 only,
# vsearch = vsearch only — no default, no auto-mix.
case "$MODE" in
    dada2|vsearch) ;;
    *) echo "Error: --mode must be 'dada2' or 'vsearch', got '$MODE'"; exit 1 ;;
esac
export MODE
if [[ -z "${PRIMER_WINDOW:-}" ]]; then
    _defaults=$(python3 "${SCRIPTS}/parameters.py" shell) || exit 2
    eval "$_defaults"
fi

# --local: read FASTQ straight from a folder (no download, no NCBI detection).
# Requires an explicit --platform; --primer-fwd/--primer-rev are optional (without
# them the entropy auto-detector is used, exactly as in download mode).
export LOCAL_MODE LOCAL_PLATFORM PRIMER_FWD PRIMER_REV
if [[ "$LOCAL_MODE" == "1" ]]; then
    case "$LOCAL_PLATFORM" in
        ILLUMINA|LS454|ION_TORRENT|PACBIO_SMRT|OXFORD_NANOPORE) ;;
        *) echo "Error: --local requires --platform one of ILLUMINA|LS454|ION_TORRENT|PACBIO_SMRT|OXFORD_NANOPORE (got '$LOCAL_PLATFORM')"; exit 1 ;;
    esac
    # --local always resolves to exactly one dataset / one platform (see
    # _local_discover_datasets), so there is nothing to parallelize ACROSS
    # datasets. Force --max-parallel 1 so the single dataset gets ALL --threads;
    # otherwise THREADS_PER_DATASET = THREADS / MAX_PARALLEL (default 2) would
    # leave (MAX_PARALLEL-1)/MAX_PARALLEL of the threads idle.
    if [[ "$MAX_PARALLEL" -ne 1 ]]; then
        echo "Note: --local processes a single dataset; forcing --max-parallel 1 (was $MAX_PARALLEL) so it uses all $THREADS thread(s)."
        MAX_PARALLEL=1
    fi
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
    METADATA=$(cd "$METADATA" && pwd -P)
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

# A worker needs at least one CPU; never launch more workers than the budget.
if [[ "$MAX_PARALLEL" -gt "$THREADS" ]]; then
    echo "Note: limiting --max-parallel $MAX_PARALLEL to the $THREADS available CPU(s)."
    MAX_PARALLEL="$THREADS"
fi
export M2D_PROFILE_CPUS="$THREADS"

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

_fastp_checked_sample() {
    local sample="$1" first="$2" second="$3" target="$4" report_stem="$5"
    local arguments=(--sample-id "$sample" --in1 "$first"
        --out1 "${target}/$(basename "$first")"
        --report-json "${READ_COUNTS_REPORTS}/fastp/${report_stem}.json"
        --report-html "${READ_COUNTS_REPORTS}/fastp/${report_stem}.html"
        --audit-dir "${READ_COUNTS_REPORTS}/fastp_guard/${sample}"
        --work-dir "${dataset_path}/tmp/fastp_guard_work"
        --db-manifest "$M2D_ADAPTER_GUARD_DB_MANIFEST"
        --cache-dir "$M2D_ADAPTER_GUARD_CACHE" --threads "$cpu"
        --min-length "$M2D_ADAPTER_GUARD_MIN_LENGTH"
        --min-identity "$M2D_ADAPTER_GUARD_MIN_IDENTITY"
        --min-coverage "$M2D_ADAPTER_GUARD_MIN_COVERAGE")
    if [[ -n "$second" ]]; then
        arguments+=(--in2 "$second" --out2 "${target}/$(basename "$second")")
    fi
    python3 "${SCRIPTS}/fastp_checked.py" "${arguments[@]}" || return $?
}

_fastp_se_adapter_remove() {
    local ori_fastq_path="$1"
    local adapter_removed_path="$2"
    mkdir -p "${READ_COUNTS_REPORTS}/fastp"
    for fq in "${ori_fastq_path}/"*.fastq*; do
        [[ -f "$fq" ]] || continue
        if [[ "${M2D_ADAPTER_GUARD_ENABLED:-0}" == 1 ]]; then
            local sample_name
            sample_name=$(basename "$fq")
            sample_name="${sample_name%.gz}"; sample_name="${sample_name%.fastq}"
            _fastp_checked_sample "$sample_name" "$fq" "" "$adapter_removed_path" "$(basename "$fq")" || return $?
            continue
        fi
        fastp -i "$fq" \
              -o "${adapter_removed_path}/$(basename "$fq")" \
              --disable_quality_filtering \
              --disable_length_filtering \
              -w "$cpu" \
              -j "${READ_COUNTS_REPORTS}/fastp/$(basename "$fq").json" \
              -h "${READ_COUNTS_REPORTS}/fastp/$(basename "$fq").html"
    done
    Audit_Fastq fastp_reads "$adapter_removed_path" RawReads
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
    local dest="${1%/}/ori_fastq"
    rm -rf "$dest"
    python3 "${SCRIPTS}/read_layout.py" stage \
        --input "${1%/}/read_layout.json" --output "$dest"
}

_obtain_reads() {
    # _obtain_reads <dataset_path> <sra_file_name> <dataset_id>
    # Download (normal mode) or symlink local files (--local).
    if [[ "${LOCAL_MODE:-0}" == "1" ]]; then
        _stage_local_reads "$1" "$3"
    else
        # Restore the owned download cache before a retry; normalized files may
        # be symlinks into it. Download verification can then reuse intact reads.
        if [[ -d "${1%/}/downloaded_fastq" ]]; then
            rm -rf "${1%/}/ori_fastq"
            mv "${1%/}/downloaded_fastq" "${1%/}/ori_fastq"
        fi
        Common_SRADownloadToFastq_MultiSource -d "$1" -a "$2" -b "$3"
    fi
}

_validate_platform_layout() {
    local actual_layout
    actual_layout=$(python3 "${SCRIPTS}/read_layout.py" layout --input "$1")
    if [[ "$platform" != "ILLUMINA" && "$actual_layout" == "paired" ]]; then
        echo "SKIP: paired-end data is unsupported for platform $platform; check --platform or separate the input dataset." >&2
        _log_status SKIPPED "$dataset_ID" "Unsupported paired-end layout for $platform"
        exit 98
    fi
}

_trim_primers() {
    local in_dir="$1" out_dir="$2"
    shift 2
    mkdir -p "$out_dir"
    local status=0
    if [[ -n "${PRIMER_FWD:-}" ]]; then
        python3 "${SCRIPTS}/explicit_primers.py" --input "$in_dir" --output "$out_dir" \
            --forward "$PRIMER_FWD" --reverse "${PRIMER_REV:-}" "$@" || status=$?
    else
        python3 "${SCRIPTS}/entropy_primer_detect.py" -i "$in_dir" -o "$out_dir" "$@" || status=$?
    fi
    if [[ -f "${out_dir}/primer_info.json" ]]; then
        cp "${out_dir}/primer_info.json" "${dataset_path}/${dataset_ID}-${MODE}-primer_info.json"
        cp "${out_dir}/primer_info.json" "${READ_COUNTS_REPORTS}/primer_info.json"
    fi
    if [[ "$status" -eq 98 ]]; then
        _log_status SKIPPED "$dataset_ID" "UNKNOWN_PRIMER; see ${dataset_ID}-${MODE}-primer_info.json"
        exit 98
    elif [[ "$status" -ne 0 ]]; then
        echo "Error: primer detection/trimming failed for $dataset_ID" >&2
        exit "$status"
    fi
    # detect-only (PacBio DADA2) does not produce trimmed reads.
    if [[ " $* " != *" --detect-only "* ]]; then
        Audit_Fastq primer_trimmed_reads "$out_dir" fastp_reads
    fi
    touch "${out_dir}/.primer_done"
}

_local_register_dataset() {
    # _local_register_dataset <id> <src_dir> — register one local dataset:
    # create its dir, map its source, and write a synthetic <id>_sra.txt
    # (Run<TAB>SampleName per unique sample prefix) so Common_CountRawReads and
    # append_summary work exactly as in download mode.
    local id="$1" src="$2"
    local dpath="${OUTPUT}/${id}"
    python3 "${SCRIPTS}/read_layout.py" validate-local --input "$src" --output "$dpath"
    mkdir -p "$dpath"
    LOCAL_SRC["$id"]="$src"
    Dataset_ID_sets+=("$id")
    echo "$id" >> "${OUTPUT}/datasets_ID.txt"
    python3 "${SCRIPTS}/read_layout.py" register --input "$src" \
        --output "${dpath}/read_layout.json" --samples "${dpath}/${id}_sra.txt"

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
    echo "  Local mode: single dataset '$(basename "$input")' from $input"
    _local_register_dataset "$(basename "$input")" "$input"
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
        echo "[ERROR] Failed to generate dataset IDs, please check your metadata file and column names for BioProject."
        exit 1
    fi

    mapfile -t Dataset_ID_sets < <(awk '{print $1}' "${OUTPUT}/datasets_ID.txt")

    if [ ${#Dataset_ID_sets[@]} -eq 0 ]; then
        echo "[ERROR] No datasets found, please check your metadata file and column names for BioProject."
        exit 1
    fi

    if ! python "${SCRIPTS}/py_16s.py" GenerateSRAsFile --FilePath "$METADATA" --Bioproject "$COL_BIOPROJECT" --SRA_Number "$COL_SRA" --OutputDir "$OUTPUT"; then
        echo "[ERROR] Failed to generate SRA file lists, please check your metadata file and column names for SRA."
        exit 1
    fi
fi

# Preserve checkpoints only for the same inputs and effective settings.
for _id in "${Dataset_ID_sets[@]}"; do
    python3 "${SCRIPTS}/pip_state.py" --dataset "${OUTPUT}/${_id}" --method "$MODE" \
        --local-source "${LOCAL_SRC[$_id]:-}" --platform "$LOCAL_PLATFORM" \
        --forward "$PRIMER_FWD" --reverse "$PRIMER_REV"
done

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

    # Skip already processed (mode-specific: an dada2 run does not block a later vsearch run)
    [[ -s "${_ds_path}/${_ds_id}-${MODE}-final-rep-seqs.qza" && -s "${_ds_path}/${_ds_id}-${MODE}-final-table.qza" ]] && continue

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

# Each profiled dataset must be a real child process so wait4 accounts for the
# complete worker, including shell orchestration and commands with no PATH shim.
# Serialize only the worker's declared context and our shell helpers. Biological
# functions are loaded from the same versioned source as the main runner.
_launch_dataset() (
    export M2D_PROFILE_PROJECT="$dataset_ID"
    export M2D_PROFILE_CPUS="$THREADS_PER_DATASET"
    if [[ -z "${M2D_PROFILE_DIR:-}" ]]; then
        _process_one_dataset
        exit $?
    fi
    local worker_file worker_rc=0
    worker_file=$(mktemp "${OUTPUT}/logs/.resource-worker-${dataset_ID}.XXXXXX") || exit 1
    {
        printf '#!/bin/bash\n'
        declare -p SCRIPT_DIR SCRIPTS OUTPUT THREADS THREADS_PER_DATASET cpu MODE \
            LOCAL_MODE LOCAL_PLATFORM PRIMER_FWD PRIMER_REV LOCAL_SRC \
            dataset_ID dataset_path sra_file_name platform PLATFORM_CACHE_FILE RUN_LOG
        printf 'source "${SCRIPTS}/AmpliconFunction.sh"\nset +e\n'
        declare -f _fastp_checked_sample _fastp_se_adapter_remove _emit_prep_done _stage_local_reads \
            _obtain_reads _validate_platform_layout _trim_primers _log_status _process_one_dataset
        printf '_process_one_dataset\n'
    } > "$worker_file"
    # M2D_PROFILE_PYTHON is resolved before installing PATH shims, avoiding an
    # accidental extra profiler around the dataset profiler itself.
    "${M2D_PROFILE_PYTHON:-python3}" "${SCRIPTS}/resource_profile.py" run \
        --scope dataset --stage pip_dataset_total -- /bin/bash "$worker_file" || worker_rc=$?
    rm -f "$worker_file"
    exit "$worker_rc"
)

worker_pids=()
_worker_failed=0
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
    if [[ -s "${dataset_path}/${dataset_ID}-${MODE}-final-rep-seqs.qza" && -s "${dataset_path}/${dataset_ID}-${MODE}-final-table.qza" ]]; then
        python3 "${SCRIPTS}/read_counts.py" recover \
            --dataset "$dataset_path" --output "$OUTPUT" --mode "$MODE" \
            --input "${dataset_path}/${dataset_ID}-${MODE}-final-table.qza"
        echo "[OK] Already processed. Skipping."
        _log_status SKIPPED "$dataset_ID" "ALREADY_DONE"
        continue
    fi

    # Wait for a slot if running at max parallel capacity
    if [[ "$MAX_PARALLEL" -gt 1 ]]; then
        while [[ "$(jobs -pr | wc -l)" -ge "$MAX_PARALLEL" ]]; do
            # The exact exit status is collected by PID below, including jobs
            # that finished before wait -n was called (which can return 127).
            wait -n 2>/dev/null || true
        done
    fi

    _process_one_dataset() {
    local _ds_start=$(date +%s)
    (
        set -e
        cd "$dataset_path"
        READ_COUNTS_STATE=$(python3 "${SCRIPTS}/read_counts.py" begin \
            --dataset "$dataset_path" --output "$OUTPUT" --mode "$MODE")
        export READ_COUNTS_STATE
        READ_COUNTS_REPORTS="$(dirname "$READ_COUNTS_STATE")/reports"
        mkdir -p "$READ_COUNTS_REPORTS"
        export READ_COUNTS_REPORTS
        trap 'rc=$?; Audit_Exit "$rc" || exit 1' EXIT

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

            if [[ -f "${fastp_path}/.primer_done" && "$n_fastp_fq" -gt 0 ]] && \
               [[ "$n_fastp_fq" -eq "$n_srr" || "$n_fastp_fq" -eq $(( n_srr * 2 )) ]]; then
                # fastp data is intact — resume from here
                echo ">>> Resuming: found $n_fastp_fq fastp files for $n_srr SRR accessions"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Resuming from checkpoint" >&3

                Audit_Counts inherit --input primer

                # Clean up downstream intermediate directories (keep step_02c as checkpoint)
                rm -rf "${dataset_path}/tmp/step_03_qza_import"
                rm -rf "${dataset_path}/tmp/step_04_qza_import_QualityFilter"
                rm -rf "${dataset_path}/tmp/step_05_dedupicate"
                rm -rf "${dataset_path}/tmp/step_05_denoise"
                rm -rf "${dataset_path}/tmp/step_06_vsearch_cli"
                rm -rf "${dataset_path}/tmp/step_06_ChimerasRemoval"
                rm -rf "${dataset_path}/tmp/step_07_cluster"
                rm -rf "${dataset_path}/tmp/temp_file"

                sequence_type=$(python3 "${SCRIPTS}/read_layout.py" layout --input "$fastp_path")
                export sequence_type
                original_sequence_type="$sequence_type"

                # Read cached quality status or re-check
                if [[ -f "$quality_cache" ]]; then
                    quality_status=$(cat "$quality_cache")
                    # Migrate pre-A3b cache value ("degraded" -> "degraded_binned").
                    [[ "$quality_status" == "degraded" ]] && quality_status="degraded_binned"
                    # Re-test if the cached token is stale/invalid (the only valid
                    # tokens are normal | degraded_binned). Prevents a stale value
                    # from mis-routing the dada2 skip / vsearch preprocess branches.
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
                    # Keep downloaded raw files: the downloader validates and
                    # reuses intact mates; local inputs are re-staged separately.
                fi

                # ── Step A: Download ──
                echo ">>> Downloading SRA data..."
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Downloading..." >&3
                if [[ "$MODE" == "dada2" && "$LOCAL_MODE" != "1" ]]; then
                    # ── Early quality probe (dada2 only) ──
                    # dada2 skips degraded/binned-quality data. Quality (Q-score
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
                        echo ">>> SKIP: degraded/binned quality is incompatible with --dada2 (DADA2)."
                        echo ">>>       Skipped before downloading the rest of the dataset."
                        _log_status SKIPPED "$dataset_ID" "dada2 mode, quality=degraded_binned (early probe)"
                        echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SKIPPED (dada2 + degraded_binned, early)" >&3
                        rm -rf "$probe_base"
                        exit 98
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
                    # vsearch mode: no early skip (quality only selects maxee vs
                    # truncation later), so download everything up front.
                    if ! _obtain_reads "$dataset_path" "${sra_file_name}" "$dataset_ID"; then
                        echo "Error: Download failed for dataset $dataset_ID" >&2
                        exit 1
                    fi
                fi

                # Count raw reads before any processing
                if [[ "$LOCAL_MODE" != "1" ]]; then
                    python3 "${SCRIPTS}/read_layout.py" normalize --input "$ori_fastq_path" \
                        --output "${dataset_path}/read_layout.json"
                fi
                _validate_platform_layout "$ori_fastq_path"
                Common_CountRawReads "$dataset_path" "$sra_file_name"

                sequence_type=$(python3 "${SCRIPTS}/read_layout.py" layout --input "$ori_fastq_path")
                echo "Sequence type: ${sequence_type^^}"
                export sequence_type
                original_sequence_type="$sequence_type"

                # ── Quality Score Diversity Check (first 3 samples) ──
                # dada2 already determined quality from the early probe above; only
                # vsearch needs it here (to pick maxee vs truncation preprocess).
                if [[ "$MODE" != "dada2" || "$LOCAL_MODE" == "1" ]]; then
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
                    mkdir -p "${READ_COUNTS_REPORTS}/fastp"
                    pair_rows=$(python3 "${SCRIPTS}/read_layout.py" pairs --input "$ori_fastq_path")
                    while IFS=$'\t' read -r sample r1 r2; do
                        if [[ "${M2D_ADAPTER_GUARD_ENABLED:-0}" == 1 ]]; then
                            _fastp_checked_sample "$sample" "$r1" "$r2" "$adapter_removed_path" "$sample" || exit $?
                            continue
                        fi
                        fastp -i "$r1" -I "$r2" \
                            -o "${adapter_removed_path}/$(basename "$r1")" \
                            -O "${adapter_removed_path}/$(basename "$r2")" \
                            --detect_adapter_for_pe \
                            --disable_quality_filtering --disable_length_filtering \
                            -w "$cpu" -j "${READ_COUNTS_REPORTS}/fastp/${sample}.json" \
                            -h "${READ_COUNTS_REPORTS}/fastp/${sample}.html"
                    done <<< "$pair_rows"
                    Audit_Fastq fastp_reads "$adapter_removed_path" RawReads
                else
                    _fastp_se_adapter_remove "$ori_fastq_path" "$adapter_removed_path" || exit $?
                touch "${adapter_removed_path}/.adapters_done"
                fi

                # ── Step C: Entropy-based primer detection & trimming ──
                mkdir -p "$fastp_path"

                _trim_primers "$adapter_removed_path" "$fastp_path"

                # Delete original and intermediate fastq files to save space
                rm -rf "$ori_fastq_path"
                rm -rf "$adapter_removed_path"
            fi

            # ── From here: same flow regardless of resume or fresh run ──
            _emit_prep_done "$_ds_start" "$dataset_ID"

            if [[ "$MODE" == "vsearch" ]]; then
                # ── vsearch mode: merge+maxee (normal) / forward-only truncation
                #    (binned) preprocess, then the shared pooled vsearch chain. ──
                Amplicon_Illumina_Vsearch_Preprocess
                VSEARCH_STRAND="plus"; export VSEARCH_STRAND   # Illumina short reads: plus only
                Amplicon_Vsearch_RunPooledChain
            elif [[ "$quality_status" == "degraded_binned" ]]; then
                # ── dada2 mode: degraded/binned quality is incompatible with DADA2 ──
                # DADA2's error model needs reliable per-base quality scores; binned
                # quality (NovaSeq/HiSeq Q-score compression, re-uploaded data) breaks
                # it. Per the method selection contract, --dada2 NEVER reroutes to vsearch —
                # such datasets are skipped here and belong to --vsearch instead.
                # (The previous vsearch orchestration for this case is preserved in
                #  AmpliconPIP_OTU_ASV_changelist.md Appendix A for the OTU back-end.)
                echo ">>> SKIP: degraded/binned quality is incompatible with --dada2 (DADA2)."
                echo ">>>       Use --vsearch for this dataset (vsearch handles binned quality)."
                _log_status SKIPPED "$dataset_ID" "dada2 mode, quality=degraded_binned"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SKIPPED (dada2 + degraded_binned)" >&3
                exit 98
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
                                echo "  [WARNING] PE retention < 50%. Falling back to SE (forward reads only)..."

                                # Clean up PE denoise outputs
                                rm -rf "${dataset_path}/tmp/step_05_denoise/"
                                rm -rf "${dataset_path}/tmp/temp_file/QualityFilter_vis/"

                                # Create SE manifest from forward reads, preserving PE sample names
                                temp_file_path="${dataset_path}/tmp/temp_file"
                                mkdir -p "$temp_file_path"
                                python3 -c "
import sys
sys.path.insert(0, sys.argv[3])
from read_layout import discover, manifest
rows = discover(sys.argv[2])
manifest([row['r1'] for row in rows], sys.argv[1], False, sample_ids=[row['sample'] for row in rows])
" "${temp_file_path}/${dataset_ID}_manifest.tsv" "$fastq_path" "$SCRIPTS"

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
                                        echo "  [WARNING] SE retention still < 50%. This dataset may have low-quality data."
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
            if [[ "$MODE" == "dada2" ]]; then
                # 454 has no DADA2 method → cannot produce ASVs. Skip (belongs to --vsearch).
                echo ">>> SKIP: LS454 (454) has no DADA2 method → not supported in --dada2. Use --vsearch."
                _log_status SKIPPED "$dataset_ID" "dada2 mode, platform=LS454"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SKIPPED (dada2: 454 unsupported)" >&3
                exit 98
            fi
            fastp_path="${dataset_path}/tmp/step_02_fastp"
            n_srr=$(wc -l < "${sra_file_name}" | tr -d ' ')

            # ── Resume checkpoint: check if fastp data is intact ──
            n_fastp_fq=0
            if [[ -d "$fastp_path" ]]; then
                n_fastp_fq=$(find "$fastp_path" -type f -name '*.fastq*' | wc -l)
            fi

            if [[ -f "${fastp_path}/.primer_done" && "$n_fastp_fq" -gt 0 ]] && [[ "$n_fastp_fq" -eq "$n_srr" ]]; then
                # 454 is always SE, so n_fastp_fq == n_srr
                echo ">>> Resuming: found $n_fastp_fq fastp files for $n_srr SRR accessions"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Resuming from checkpoint" >&3

                Audit_Counts inherit --input primer

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
                    # Keep downloaded raw files: the downloader validates and
                    # reuses intact mates; local inputs are re-staged separately.
                fi

                echo ">>> Downloading SRA data..."
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Downloading..." >&3
                if ! _obtain_reads "$dataset_path" "${sra_file_name}" "$dataset_ID"; then
                    echo "Error: Download failed for dataset $dataset_ID" >&2
                    exit 1
                fi

                if [[ "$LOCAL_MODE" != "1" ]]; then
                    python3 "${SCRIPTS}/read_layout.py" normalize --input "$ori_fastq_path" \
                        --output "${dataset_path}/read_layout.json"
                fi
                _validate_platform_layout "$ori_fastq_path"
                Common_CountRawReads "$dataset_path" "$sra_file_name"

                adapter_removed_path="${dataset_path}/tmp/step_01_adapter_removed"
                mkdir -p "$adapter_removed_path"

                _fastp_se_adapter_remove "$ori_fastq_path" "$adapter_removed_path" || exit $?
                touch "${adapter_removed_path}/.adapters_done"

                mkdir -p "$fastp_path"
                _trim_primers "$adapter_removed_path" "$fastp_path"

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

            Audit_Fastq vsearch_adaptive_trimmed_reads "$adaptive_trim_path" primer_trimmed_reads

            # Clean up pre-trim FASTQ
            rm -rf "$fastp_path"

            # ── OTU back-end (unified pooled vsearch chain) ──
            # 454 is vsearch-only (no DADA2 method). §4.x (user-accepted) behaviour
            # change: the old QIIME2 dedup→chimera→cluster path is replaced by the
            # shared UNOISE3 → cluster_fast 97% → map-back chain. Abundance now
            # comes from read map-back (not cluster size), with added UNOISE3
            # denoising. NOTE for review: the old LS454_QualityControlForQZA
            # q-score/length filter is dropped; adaptive_tail_trim handles the
            # 3' N-tail but there is no explicit quality (maxee) filter for 454.
            fastq_path="$adaptive_trim_path"
            export fastq_path
            sequence_type="single"; export sequence_type
            VSEARCH_STRAND="plus"; export VSEARCH_STRAND
            Amplicon_Vsearch_RunPooledChain

        elif [[ "$platform" == "ION_TORRENT" ]]; then
            fastp_path="${dataset_path}/tmp/step_02_fastp"
            n_srr=$(wc -l < "${sra_file_name}" | tr -d ' ')

            # ── Resume checkpoint: check if fastp data is intact ──
            n_fastp_fq=0
            if [[ -d "$fastp_path" ]]; then
                n_fastp_fq=$(find "$fastp_path" -type f -name '*.fastq*' | wc -l)
            fi

            if [[ -f "${fastp_path}/.primer_done" && "$n_fastp_fq" -gt 0 ]] && [[ "$n_fastp_fq" -eq "$n_srr" ]]; then
                echo ">>> Resuming: found $n_fastp_fq fastp files for $n_srr SRR accessions"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Resuming from checkpoint" >&3
                Audit_Counts inherit --input primer
                rm -rf "${dataset_path}/tmp/step_03_qza_import"
                rm -rf "${dataset_path}/tmp/step_04_qza_import_QualityFilter"
                rm -rf "${dataset_path}/tmp/step_05_denoise"
                rm -rf "${dataset_path}/tmp/temp_file"
            else
                if [[ -d "${dataset_path}/tmp" ]]; then
                    echo ">>> No valid fastp checkpoint ($n_fastp_fq files, expected $n_srr). Cleaning and re-running..."
                    rm -rf "${dataset_path}/tmp"
                    # Keep downloaded raw files: the downloader validates and
                    # reuses intact mates; local inputs are re-staged separately.
                fi

                echo ">>> Downloading SRA data..."
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Downloading..." >&3
                if ! _obtain_reads "$dataset_path" "${sra_file_name}" "$dataset_ID"; then
                    echo "Error: Download failed for dataset $dataset_ID" >&2
                    exit 1
                fi

                if [[ "$LOCAL_MODE" != "1" ]]; then
                    python3 "${SCRIPTS}/read_layout.py" normalize --input "$ori_fastq_path" \
                        --output "${dataset_path}/read_layout.json"
                fi
                _validate_platform_layout "$ori_fastq_path"
                Common_CountRawReads "$dataset_path" "$sra_file_name"

                adapter_removed_path="${dataset_path}/tmp/step_01_adapter_removed"
                mkdir -p "$adapter_removed_path"

                _fastp_se_adapter_remove "$ori_fastq_path" "$adapter_removed_path" || exit $?
                touch "${adapter_removed_path}/.adapters_done"

                mkdir -p "$fastp_path"
                _trim_primers "$adapter_removed_path" "$fastp_path"

                rm -rf "$ori_fastq_path"
                rm -rf "$adapter_removed_path"
            fi

            sequence_type="single"
            export sequence_type
            original_sequence_type="$sequence_type"

            _emit_prep_done "$_ds_start" "$dataset_ID"

            if [[ "$MODE" == "vsearch" ]]; then
                # ── vsearch: configurable extra 5' trimming + maxee → shared pooled vsearch chain ──
                Amplicon_IonTorrent_Vsearch_Preprocess
                VSEARCH_STRAND="plus"; export VSEARCH_STRAND
                Amplicon_Vsearch_RunPooledChain
            else
                # ── DADA2: QIIME2 Import → Quality filter → DADA2 denoise-pyro ──
                # Additional Ion Torrent trimming is explicit and configurable.
                # The default is zero after primer removal.
                # trunc-len is computed automatically from QC visualization.
                fastq_path="$fastp_path"
                export fastq_path
                Common_SanitizeFastq
                Amplicon_Common_MakeManifestFileForQiime2
                Amplicon_Common_ImportFastqToQiime2
                Amplicon_IonTorrent_QualityControlForQZA
                Amplicon_Illumina_DenosingDada2 -s "${DADA2_ION_TRIM_LEFT}"
                Amplicon_Common_FinalFilesCleaning
            fi

        elif [[ "$platform" == "OXFORD_NANOPORE" ]]; then
            if [[ "$MODE" == "dada2" ]]; then
                # Single-base ASV resolution is conceptually invalid for ONT
                # (~5-10% error + indels): true sequences explode into a cloud of
                # spurious variants. Skip in --dada2; ONT belongs to --vsearch.
                echo ">>> SKIP: ONT single-base ASV is invalid (5-10% error+indels) → not supported in --dada2. Use --vsearch."
                _log_status SKIPPED "$dataset_ID" "dada2 mode, platform=OXFORD_NANOPORE"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SKIPPED (dada2: ONT unsupported)" >&3
                exit 98
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
                echo "[ERROR] ONT processing requires missing tools: ${ont_missing[*]}" >&2
                echo "   Install via: conda install -c bioconda ${ont_missing[*]}  (or module load)" >&2
                exit 1
            fi

            primer_trim_path="${dataset_path}/tmp/step_02_primer"
            chopper_path="${dataset_path}/tmp/step_03_chopper"

            # ── Resume checkpoint: chopper output complete? ──
            if [[ -f "${chopper_path}/.chopper_done" ]]; then
                echo ">>> Resuming: found completed chopper-filtered reads"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Resuming from checkpoint" >&3
                Audit_Counts inherit --input chopper
                # Clean only downstream directories
                rm -rf "${dataset_path}/tmp/step_06_ont"
                rm -rf "${dataset_path}/tmp/step_07_cluster"
                rm -rf "${dataset_path}/tmp/temp_file"
            else
                if [[ -d "${dataset_path}/tmp" ]]; then
                    echo ">>> No valid chopper checkpoint. Cleaning and re-running..."
                    rm -rf "${dataset_path}/tmp"
                    # Keep downloaded raw files: the downloader validates and
                    # reuses intact mates; local inputs are re-staged separately.
                fi

                # ── Step A: Download ──
                echo ">>> Downloading SRA data..."
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Downloading..." >&3
                if ! _obtain_reads "$dataset_path" "${sra_file_name}" "$dataset_ID"; then
                    echo "Error: Download failed for dataset $dataset_ID" >&2
                    exit 1
                fi

                if [[ "$LOCAL_MODE" != "1" ]]; then
                    python3 "${SCRIPTS}/read_layout.py" normalize --input "$ori_fastq_path" \
                        --output "${dataset_path}/read_layout.json"
                fi
                _validate_platform_layout "$ori_fastq_path"
                Common_CountRawReads "$dataset_path" "$sra_file_name"

                # ── Step B: Remove sequencing adapters with fastp (SE) ──
                adapter_removed_path="${dataset_path}/tmp/step_01_adapter_removed"
                mkdir -p "$adapter_removed_path"
                _fastp_se_adapter_remove "$ori_fastq_path" "$adapter_removed_path" || exit $?
                touch "${adapter_removed_path}/.adapters_done"

                # ── Step C: Entropy-based primer detection & trimming ──
                mkdir -p "$primer_trim_path"
                _trim_primers "$adapter_removed_path" "$primer_trim_path"

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
            # CCS detection is repeated below; DADA2 itself performs the trimming.
            n_adapter_fq=0
            if [[ -d "$adapter_removed_path" ]]; then
                n_adapter_fq=$(find "$adapter_removed_path" -type f -name '*.fastq*' ! -name 'fastp.*' | wc -l)
            fi

            if [[ -f "${adapter_removed_path}/.adapters_done" && "$n_adapter_fq" -gt 0 ]] && [[ "$n_adapter_fq" -eq "$n_srr" ]]; then
                echo ">>> Resuming: found $n_adapter_fq adapter-removed files for $n_srr SRR accessions"
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Resuming from checkpoint" >&3
                Audit_Counts inherit --input fastp
                rm -rf "${dataset_path}/tmp/step_03_qza_import"
                rm -rf "${dataset_path}/tmp/step_04_qza_import_QualityFilter"
                rm -rf "${dataset_path}/tmp/step_05_denoise"
                rm -rf "${dataset_path}/tmp/temp_file"
            else
                if [[ -d "${dataset_path}/tmp" ]]; then
                    echo ">>> No valid checkpoint ($n_adapter_fq files, expected $n_srr). Cleaning and re-running..."
                    rm -rf "${dataset_path}/tmp"
                    # Keep downloaded raw files: the downloader validates and
                    # reuses intact mates; local inputs are re-staged separately.
                fi

                echo ">>> Downloading SRA data..."
                echo "[$(date '+%H:%M:%S')] [${dataset_ID}] [2/3] Downloading..." >&3
                if ! _obtain_reads "$dataset_path" "${sra_file_name}" "$dataset_ID"; then
                    echo "Error: Download failed for dataset $dataset_ID" >&2
                    exit 1
                fi

                if [[ "$LOCAL_MODE" != "1" ]]; then
                    python3 "${SCRIPTS}/read_layout.py" normalize --input "$ori_fastq_path" \
                        --output "${dataset_path}/read_layout.json"
                fi
                _validate_platform_layout "$ori_fastq_path"
                Common_CountRawReads "$dataset_path" "$sra_file_name"

                mkdir -p "$adapter_removed_path"
                _fastp_se_adapter_remove "$ori_fastq_path" "$adapter_removed_path" || exit $?
                touch "${adapter_removed_path}/.adapters_done"

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
                echo "[ERROR] No FASTQ files found after adapter removal"
                exit 1
            fi

            echo ">>> Checking read lengths from first sample: $(basename "$first_fq")"
            long_read_ratio=$(python3 -c "
import sys, gzip, os

fq_path = sys.argv[1]
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
" "$first_fq")
            echo "  Reads > 1400 bp ratio: ${long_read_ratio} (from first 1000 reads)"

            # ── Sub-condition A: Full-length 16S CCS reads (majority > 1400bp) ──
            if python3 -c "sys_exit = __import__('sys').exit; sys_exit(0 if float('${long_read_ratio}') > 0.5 else 1)"; then
                echo ">>> Full-length 16S detected (>50% reads > 1400bp)."

                if [[ "$MODE" == "vsearch" ]]; then
                # ── OTU: length-window + maxee_rate → shared pooled chain (strand both) ──
                pacbio_primer_path="${dataset_path}/tmp/step_02_primer"
                _trim_primers "$adapter_removed_path" "$pacbio_primer_path" --mixed-orientation
                adapter_removed_path="$pacbio_primer_path"
                Amplicon_Pacbio_Vsearch_Preprocess
                VSEARCH_STRAND="both"; export VSEARCH_STRAND
                Amplicon_Vsearch_RunPooledChain
                else
                # Detect without trimming: denoise-ccs needs the forward primer for orientation.
                pacbio_primer_path="${dataset_path}/tmp/step_02_primer"
                _trim_primers "$adapter_removed_path" "$pacbio_primer_path" --detect-only
                primer_spec=$(python3 "${SCRIPTS}/pacbio_primers.py" \
                    "${pacbio_primer_path}/primer_info.json") || {
                    _log_status SKIPPED "$dataset_ID" "DADA2 CCS requires a known or explicit forward primer for orientation"
                    exit 98
                }
                IFS=$'\t' read -r primer_front primer_adapter <<< "$primer_spec"
                # Import adapter-removed reads directly into QIIME2
                fastq_path="$adapter_removed_path"
                export fastq_path
                Common_SanitizeFastq
                Amplicon_Common_MakeManifestFileForQiime2
                Amplicon_Common_ImportFastqToQiime2
                Amplicon_Pacbio_QualityControlForQZA

                # DADA2 denoise-ccs handles orientation and primer removal once
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
                exit 98
            fi

        else
            echo "SKIP: unsupported platform: $platform"
            _log_status SKIPPED "$dataset_ID" "Unsupported platform: $platform"
            exit 98
        fi

        # Preserve the final per-sample abundances in the existing summary.csv.
        Audit_Table "${MODE}_final_reads" "${dataset_path}/${dataset_ID}-${MODE}-final-table.qza"

        _log_status SUCCESS "$dataset_ID" "Platform: $platform"

    )
    local _rc=$?
    local _ds_end=$(date +%s)
    local _total=$(( _ds_end - _ds_start ))
    local _total_fmt="$(( _total / 60 ))m$(( _total % 60 ))s"
    if [[ $_rc -eq 98 ]]; then
        echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SKIPPED (${_total_fmt})" >&3
    elif [[ $_rc -eq 99 ]]; then
        _log_status SKIPPED "$dataset_ID" "Untrustworthy single-sample abundance; see ${dataset_ID}-UNTRUSTABLE.txt"
        # Untrustworthy data details are saved inside the dataset.
        echo "[WARNING] Skipped $dataset_ID — untrustworthy data (see ${dataset_ID}-UNTRUSTABLE.txt)"
        echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SKIPPED-UNTRUSTABLE (${_total_fmt})" >&3
    elif [[ $_rc -ne 0 ]]; then
        echo "[ERROR] Pipeline failed for $dataset_ID — skipping to next dataset"
        _log_status FAILED "$dataset_ID" "see logs/${dataset_ID}.log"
        echo "[$(date '+%H:%M:%S')] [${dataset_ID}] FAILED (${_total_fmt}) - see logs/${dataset_ID}.log" >&3
    else
        echo "[$(date '+%H:%M:%S')] [${dataset_ID}] SUCCESS (${_total_fmt})" >&3
    fi
    # Intentional biological skips remain non-fatal. Actual worker failures
    # propagate after every parallel dataset has had a chance to finish.
    if [[ "$_rc" -ne 0 && "$_rc" -ne 98 && "$_rc" -ne 99 ]]; then
        return "$_rc"
    fi
    return 0
    }

    set +e
    if [[ "$MAX_PARALLEL" -gt 1 ]]; then
        _launch_dataset >> "$log_file" 2>&1 &
        worker_pids+=("$!")
    else
        _launch_dataset 2>&1 | tee "$log_file"
        [[ "${PIPESTATUS[0]}" -eq 0 ]] || _worker_failed=1
    fi
    set -e
done

# Wait for all remaining background jobs (parallel mode)
if [[ "$MAX_PARALLEL" -gt 1 ]]; then
    echo ">>> Waiting for remaining background datasets to finish..."
    for worker_pid in "${worker_pids[@]}"; do
        wait "$worker_pid" || _worker_failed=1
    done
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

if [[ "$n_failed" -gt 0 || "$_worker_failed" -ne 0 ]]; then
    exit 1
fi
