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
    -t, --threads INT          Number of CPU threads (default: 4)
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

# Strip trailing slashes from paths
METADATA="${METADATA%/}"
OUTPUT="${OUTPUT%/}"

# Validate input
if [[ -z "$METADATA" ]]; then
    echo "Error: --metadata is required"
    exit 1
fi

if [[ "$METADATA" != /* ]]; then
    echo "Error: --metadata must be a full absolute path (starting with '/')"
    exit 1
fi

if [[ ! -f "$METADATA" ]]; then
    echo "Error: Metadata file '$METADATA' does not exist"
    exit 1
fi

if [[ -z "$OUTPUT" ]]; then
    OUTPUT=$(dirname "$METADATA")
fi

export cpu=$THREADS

# Source function library (AmpliconFunction.sh contains all processing functions)
if [[ -f "${SCRIPTS}/AmpliconFunction.sh" ]]; then
    source "${SCRIPTS}/AmpliconFunction.sh"
else
    echo "Error: AmpliconFunction.sh not found at ${SCRIPTS}/AmpliconFunction.sh"
    exit 1
fi

cd "$OUTPUT" || exit 1
mkdir -p "${OUTPUT}/analysis/"

################################################################################
#                     PRE-FLIGHT: QIIME2 SANITY CHECK                          #
################################################################################

echo "========================================="
echo "Pre-flight: Checking QIIME2 availability"
echo "========================================="

# Use "qiime tools import --help" instead of "qiime info" for the pre-flight
# check. "qiime info" tries to write a vendored parallel_config file, which
# fails on read-only HPC container filesystems (e.g., Tykky on Puhti) with
# OSError: [Errno 30] Read-only file system. "tools import --help" exercises
# the same plugin loading path without triggering that config write.
if ! qiime tools import --help > /dev/null 2>&1; then
    echo ""
    echo "❌ QIIME2 PRE-FLIGHT FAILED"
    echo ""
    # Capture the actual error for diagnostics
    qiime_err=$(qiime tools import --help 2>&1 || true)
    if echo "$qiime_err" | grep -q "numpy"; then
        echo "ROOT CAUSE: NumPy version incompatibility."
        echo "  The installed NumPy $(python3 -c 'import numpy; print(numpy.__version__)' 2>/dev/null || echo '(unknown)') is incompatible"
        echo "  with compiled extensions (skbio, etc.) that require NumPy 1.x."
        echo ""
        echo "FIX: Downgrade NumPy in your conda environment:"
        echo "  conda install 'numpy<2' or pip install 'numpy<2'"
    elif echo "$qiime_err" | grep -q "Read-only file system"; then
        echo "ROOT CAUSE: Read-only filesystem prevents QIIME2 config writes."
        echo "  QIIME2 is trying to write to: $(echo "$qiime_err" | grep -oP "(?<=Read-only file system: ').*(?=')" || echo '(unknown path)')"
        echo ""
        echo "FIX: Set a writable QIIME2 config path before running the pipeline:"
        echo "  export QIIME2_CONFIG=\$HOME/.qiime2_config.toml"
    else
        echo "QIIME2 failed to initialize. Error output:"
        echo "$qiime_err" | head -20
    fi
    echo ""
    echo "Aborting pipeline — all datasets would fail at the QIIME2 import step."
    exit 1
fi
echo "✓ QIIME2 is functional"

################################################################################
#                          PHASE 1: DATASET PREPARATION                        #
################################################################################

echo "========================================="
echo "PHASE 1: Dataset Preparation"
echo "Started: $(date)"
echo "========================================="

# Step 1: Generate dataset ID list
if ! python "${SCRIPTS}/py_16s.py" GenerateDatasetsIDsFile --FilePath "$METADATA" --Bioproject "$COL_BIOPROJECT" --OutputDir "$OUTPUT"; then
    echo "❌ ERROR: Failed to generate dataset IDs"
    exit 1
fi

mapfile -t Dataset_ID_sets < <(awk '{print $1}' "${OUTPUT}/datasets_ID.txt")

if [ ${#Dataset_ID_sets[@]} -eq 0 ]; then
    echo "❌ ERROR: No datasets found"
    exit 1
fi

# Step 2: Generate SRA file list
if ! python "${SCRIPTS}/py_16s.py" GenerateSRAsFile --FilePath "$METADATA" --Bioproject "$COL_BIOPROJECT" --SRA_Number "$COL_SRA" --OutputDir "$OUTPUT"; then
    echo "❌ ERROR: Failed to generate SRA file lists"
    exit 1
fi

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
summary_csv="${OUTPUT}/summary.csv"

: > "$failed_log"
: > "$success_log"
: > "$skipped_log"

for i in "${!Dataset_ID_sets[@]}"; do
    dataset_ID="${Dataset_ID_sets[$i]}"
    dataset_path="${OUTPUT}/${dataset_ID}"
    sra_file_name="${dataset_ID}_sra.txt"
    platform="Unknown" 

    echo "----------------------------------------"
    echo "Dataset $((i+1))/${#Dataset_ID_sets[@]}: $dataset_ID"
    
    # Check if already processed
    if [ -f "${dataset_path}/${dataset_ID}-final-rep-seqs.qza" ]; then
        echo "✓ Already processed. Skipping."
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - ALREADY_DONE" >> "$skipped_log"
        continue
    fi

    # Capture stderr from the subshell so we can log meaningful error messages
    dataset_errlog="${dataset_path}/${dataset_ID}_pipeline.err"

    set +e
    (
        set -e
        cd "$dataset_path"

        # 1. Dynamic Platform Detection (BEFORE downloading)
        echo ">>> Detecting sequencing platform..."
        first_srr=$(awk 'NR==1 {print $1}' "${sra_file_name}")

        # Query API for platform (pass BioProject ID for CNCB/CRR accessions)
        if [[ "$first_srr" =~ ^CRR ]]; then
            # CRR accession - need to pass BioProject ID to CNCB API
            echo "  CNCB accession detected, using BioProject: $dataset_ID"
            platform=$(python "${SCRIPTS}/py_16s.py" get_sequencing_platform --srr_id "$first_srr" --bioproject_id "$dataset_ID")
        else
            # NCBI accession (SRR/ERR/DRR)
            platform=$(python "${SCRIPTS}/py_16s.py" get_sequencing_platform --srr_id "$first_srr")
        fi

        echo "Detected platform: $platform"

        # 2. Platform-specific pipeline (download + processing)
        export dataset_path
        export dataset_name="$dataset_ID"
        ori_fastq_path="${dataset_path}/ori_fastq"

        if [[ "$platform" == "ILLUMINA" ]]; then
            # ── Step A: Download ──
            echo ">>> Downloading SRA data..."
            if ! Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}"; then
                echo "Error: Download failed for dataset $dataset_ID" >&2
                exit 1
            fi

            # Count raw reads before any processing
            Common_CountRawReads "$dataset_path" "$sra_file_name"

            # Detect PE/SE after download (check for _1/_2 paired files)
            r1_count=$(find "$ori_fastq_path" -type f -name '*_1.fastq*' 2>/dev/null | wc -l)
            r2_count=$(find "$ori_fastq_path" -type f -name '*_2.fastq*' 2>/dev/null | wc -l)
            if [ "$r1_count" -gt 0 ] && [ "$r1_count" -eq "$r2_count" ]; then
                sequence_type="paired"
            else
                sequence_type="single"
            fi
            echo "Sequence type: ${sequence_type^^}"
            export sequence_type

            # ── Step B: Remove sequencing adapters with fastp ──
            adapter_removed_path="${dataset_path}/tmp/step_01_adapter_removed"
            mkdir -p "$adapter_removed_path"

            if [[ "$sequence_type" == "paired" ]]; then
                # PE: find R1/R2 pairs and run fastp in PE mode
                # Output filenames are normalised to _1/_2 so that all
                # downstream tools (entropy_primer_detect, mk_manifest_PE,
                # QIIME2 import) see a consistent naming convention.
                pe_done=false
                for r1 in "${ori_fastq_path}/"*_R1*.fastq*; do
                    [[ -f "$r1" ]] || continue
                    r2="${r1/_R1/_R2}"
                    [[ -f "$r2" ]] || continue
                    # Normalise _R1 → _1, _R2 → _2
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
                    # Fallback: files already use _1/_2 pattern
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
                # SE: run fastp on each file
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

            # ── Step B: Entropy-based primer detection & trimming ──
            fastp_path="${dataset_path}/tmp/step_02_fastp"
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

            # Run Illumina pipeline (dada2 denoise-paired for PE, denoise-single for SE)
            fastq_path="$fastp_path"
            export fastq_path
            Amplicon_Common_MakeManifestFileForQiime2
            Amplicon_Common_ImportFastqToQiime2
            Amplicon_Illumina_QualityControlForQZA
            Amplicon_Illumina_DenosingDada2
            Amplicon_Common_FinalFilesCleaning

        elif [[ "$platform" == "LS454" ]]; then
            # ── Step A: Download with normal prefetch + fasterq-dump (same as Illumina) ──
            echo ">>> Downloading SRA data (454)..."
            if ! Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}"; then
                echo "Error: Download failed for dataset $dataset_ID" >&2
                exit 1
            fi

            # Count raw reads before any processing
            Common_CountRawReads "$dataset_path" "$sra_file_name"

            sequence_type="single"
            export sequence_type

            # ── Step B: Remove sequencing adapters with fastp (SE mode) ──
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

            # ── Step C: Entropy-based primer detection & trimming (same as Illumina) ──
            fastp_path="${dataset_path}/tmp/step_02_fastp"
            mkdir -p "$fastp_path"

            python3 "${SCRIPTS}/entropy_primer_detect.py" \
                -i "$adapter_removed_path" \
                -o "$fastp_path" || {
                echo "  ✗ Entropy primer detection failed"
                exit 1
            }

            rm -rf "$ori_fastq_path"
            rm -rf "$adapter_removed_path"

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
            Amplicon_Common_MakeManifestFileForQiime2
            Amplicon_Common_ImportFastqToQiime2
            Amplicon_LS454_QualityControlForQZA
            Amplicon_LS454_Deduplication
            Amplicon_LS454_ChimerasRemoval
            Amplicon_LS454_ClusterDenovo
            Amplicon_LS454_FilterLowFreqOTUs
            Amplicon_Common_FinalFilesCleaning

        elif [[ "$platform" == "ION_TORRENT" ]]; then
            # ── Step A: Download with normal prefetch + fasterq-dump (same as 454) ──
            echo ">>> Downloading SRA data (Ion Torrent)..."
            if ! Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}"; then
                echo "Error: Download failed for dataset $dataset_ID" >&2
                exit 1
            fi

            # Count raw reads before any processing
            Common_CountRawReads "$dataset_path" "$sra_file_name"

            sequence_type="single"
            export sequence_type

            # ── Step B: Remove sequencing adapters with fastp (SE mode) ──
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

            # ── Step C: Entropy-based primer detection & trimming (same as Illumina) ──
            fastp_path="${dataset_path}/tmp/step_02_fastp"
            mkdir -p "$fastp_path"

            python3 "${SCRIPTS}/entropy_primer_detect.py" \
                -i "$adapter_removed_path" \
                -o "$fastp_path" || {
                echo "  ✗ Entropy primer detection failed"
                exit 1
            }

            rm -rf "$ori_fastq_path"
            rm -rf "$adapter_removed_path"

            # ── Step D: QIIME2 Import → Quality filter → DADA2 denoise-pyro ──
            # Ion Torrent signal instability in the first ~10bp is handled by
            # DADA2 denoise-pyro --p-trim-left 10 (passed via -s 10).
            # trunc-len is computed automatically from QC visualization.
            fastq_path="$fastp_path"
            export fastq_path
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
            # ── Step A: Download ──
            echo ">>> Downloading SRA data..."
            if ! Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}"; then
                echo "Error: Download failed for dataset $dataset_ID" >&2
                exit 1
            fi

            # Count raw reads before any processing
            Common_CountRawReads "$dataset_path" "$sra_file_name"

            sequence_type="single"
            export sequence_type

            # ── Step B: Remove sequencing adapters with fastp ──
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

            # Clean up original FASTQ to save space
            rm -rf "$ori_fastq_path"

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
            --output_csv "$summary_csv"

        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - SUCCESS - Platform: $platform" >> "$success_log"

    ) 2> >(tee "$dataset_errlog" >&2)
    rc=$?
    set -e
    if [[ $rc -ne 0 ]]; then
        # Extract a meaningful error reason from the captured stderr
        err_reason=""
        if [[ -f "$dataset_errlog" ]]; then
            if grep -q "numpy" "$dataset_errlog" 2>/dev/null; then
                err_reason="NumPy version incompatibility (QIIME2 environment broken)"
            elif grep -q "Download failed\|All download sources failed" "$dataset_errlog" 2>/dev/null; then
                err_reason="Data download failed"
            elif grep -q "No FASTQ files" "$dataset_errlog" 2>/dev/null; then
                err_reason="No FASTQ files produced"
            else
                # Use last non-empty stderr line as reason
                err_reason=$(grep -v '^$' "$dataset_errlog" 2>/dev/null | tail -1 | head -c 200)
            fi
        fi
        echo "❌ Pipeline failed for $dataset_ID — skipping to next dataset"
        [[ -n "$err_reason" ]] && echo "   Reason: $err_reason"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - FAILED - ${err_reason:-unknown}" >> "$failed_log"
    else
        # Clean up error log on success
        rm -f "$dataset_errlog"
    fi
done

################################################################################
#                          FINAL SUMMARY                                       #
################################################################################

n_success=$(wc -l < "$success_log" 2>/dev/null || echo 0)
n_failed=$(wc -l < "$failed_log" 2>/dev/null || echo 0)
n_skipped=$(wc -l < "$skipped_log" 2>/dev/null || echo 0)

echo "========================================="
echo "ALL DONE"
echo "========================================="
echo "  Success:  $n_success"
echo "  Failed:   $n_failed"
echo "  Skipped:  $n_skipped"
echo "========================================="

if [[ "$n_failed" -gt 0 ]]; then
    echo ""
    echo "Failed datasets:"
    cat "$failed_log"
    echo ""
fi