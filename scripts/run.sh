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
#                          FIND REFERENCE FILE                                 #
################################################################################

# Find the 16S reference file in multiple possible locations
find_reference_file() {
    local ref_file="J01859.1.fna"
    local search_paths=(
        "${SCRIPT_DIR}/docs/${ref_file}"                    # Development/git repo
        "/usr/local/share/Meta2Data/docs/${ref_file}"       # System install
        "${HOME}/.local/share/Meta2Data/docs/${ref_file}"   # User install
        "${CONDA_PREFIX}/share/Meta2Data/docs/${ref_file}"  # Conda install
        "${PREFIX}/share/Meta2Data/docs/${ref_file}"        # Alternative prefix
    )

    for path in "${search_paths[@]}"; do
        if [[ -f "$path" ]]; then
            echo "$path"
            return 0
        fi
    done

    echo "ERROR: 16S reference file '${ref_file}' not found in any location:" >&2
    for path in "${search_paths[@]}"; do
        echo "  - $path" >&2
    done
    return 1
}

# Find and export reference file path
REF_16S_PATH=$(find_reference_file)
if [[ $? -ne 0 ]]; then
    echo "❌ ERROR: Cannot find 16S reference sequence file"
    exit 1
fi
echo "Using 16S reference: $REF_16S_PATH"

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
    dataset_path="${OUTPUT}/${dataset_ID}/"
    sra_file_name="${dataset_ID}_sra.txt"
    platform="Unknown" 

    echo "----------------------------------------"
    echo "Dataset $((i+1))/${#Dataset_ID_sets[@]}: $dataset_ID"
    
    # Check if already processed
    if [ -f "${dataset_path}${dataset_ID}-final-rep-seqs.qza" ]; then
        echo "✓ Already processed. Skipping."
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - ALREADY_DONE" >> "$skipped_log"
        continue
    fi

    {
        cd "$dataset_path" || exit 1

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
        ori_fastq_path="${dataset_path}ori_fastq/"

        if [[ "$platform" == "ILLUMINA" ]]; then
            # ── Step A: Download ──
            echo ">>> Downloading SRA data..."
            Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}"

            # Count raw reads before any processing
            Common_CountRawReads "$dataset_path" "$sra_file_name"

            # Detect PE/SE after download
            line_count=$(wc -l < "${dataset_path}${sra_file_name}")
            file_count=$(find "$ori_fastq_path" -type f 2>/dev/null | wc -l)
            if [ $((line_count * 2)) -eq $file_count ]; then
                sequence_type="paired"
            else
                sequence_type="single"
            fi
            echo "Sequence type: ${sequence_type^^}"
            export sequence_type

            # ── Step B: Remove sequencing adapters with fastp ──
            echo ">>> Removing sequencing adapters with fastp..."
            adapter_removed_path="${dataset_path}temp/step_01_adapter_removed/"
            mkdir -p "$adapter_removed_path"

            if [[ "$sequence_type" == "paired" ]]; then
                # PE: find R1/R2 pairs and run fastp in PE mode
                # Output filenames are normalised to _1/_2 so that all
                # downstream tools (entropy_primer_detect, mk_manifest_PE,
                # QIIME2 import) see a consistent naming convention.
                pe_done=false
                for r1 in "$ori_fastq_path"*_R1*.fastq*; do
                    [[ -f "$r1" ]] || continue
                    r2="${r1/_R1/_R2}"
                    [[ -f "$r2" ]] || continue
                    # Normalise _R1 → _1, _R2 → _2
                    r1_out=$(basename "$r1"); r1_out="${r1_out/_R1/_1}"
                    r2_out=$(basename "$r2"); r2_out="${r2_out/_R2/_2}"
                    fastp -i "$r1" -I "$r2" \
                          -o "${adapter_removed_path}${r1_out}" \
                          -O "${adapter_removed_path}${r2_out}" \
                          --detect_adapter_for_pe \
                          --disable_quality_filtering \
                          --disable_length_filtering \
                          -w "$cpu" \
                          -j "${adapter_removed_path}fastp.json" \
                          -h "${adapter_removed_path}fastp.html"
                    echo "  ✓ Adapter removal: $(basename "$r1") + $(basename "$r2") → ${r1_out} + ${r2_out}"
                    pe_done=true
                done
                if [[ "$pe_done" == false ]]; then
                    # Fallback: files already use _1/_2 pattern
                    for r1 in "$ori_fastq_path"*_1.fastq*; do
                        [[ -f "$r1" ]] || continue
                        r2="${r1/_1.fastq/_2.fastq}"
                        [[ -f "$r2" ]] || continue
                        fastp -i "$r1" -I "$r2" \
                              -o "${adapter_removed_path}$(basename "$r1")" \
                              -O "${adapter_removed_path}$(basename "$r2")" \
                              --detect_adapter_for_pe \
                              --disable_quality_filtering \
                              --disable_length_filtering \
                              -w "$cpu" \
                              -j "${adapter_removed_path}fastp.json" \
                              -h "${adapter_removed_path}fastp.html"
                        echo "  ✓ Adapter removal: $(basename "$r1") + $(basename "$r2")"
                    done
                fi
            else
                # SE: run fastp on each file
                for fq in "$ori_fastq_path"*.fastq*; do
                    [[ -f "$fq" ]] || continue
                    fastp -i "$fq" \
                          -o "${adapter_removed_path}$(basename "$fq")" \
                          --disable_quality_filtering \
                          --disable_length_filtering \
                          -w "$cpu" \
                          -j "${adapter_removed_path}fastp.json" \
                          -h "${adapter_removed_path}fastp.html"
                    echo "  ✓ Adapter removal: $(basename "$fq")"
                done
            fi
            echo "✓ Adapter removal completed"

            # ── Step B: Entropy-based primer detection & trimming ──
            echo ">>> Detecting and removing primers (entropy method)..."
            fastp_path="${dataset_path}temp/step_02_fastp/"
            mkdir -p "$fastp_path"

            python3 "${SCRIPTS}/entropy_primer_detect.py" \
                -i "$adapter_removed_path" \
                -o "$fastp_path" || {
                echo "  ✗ Entropy primer detection failed"
                exit 1
            }

            echo "✓ Entropy primer detection/removal completed"

            # Delete original and intermediate fastq files to save space
            echo ">>> Cleaning up intermediate files..."
            rm -rf "$ori_fastq_path"
            rm -rf "$adapter_removed_path"
            echo "✓ Removed ori_fastq/ and adapter_removed/ to save space"

            # Run Illumina pipeline (dada2 denoise-paired for PE, denoise-pyro for SE)
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
            Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}"

            # Count raw reads before any processing
            Common_CountRawReads "$dataset_path" "$sra_file_name"

            # 454 is always single-end
            sequence_type="single"
            echo "Sequence type: ${sequence_type^^}"
            export sequence_type

            # ── Step B: Remove sequencing adapters with fastp (SE mode) ──
            echo ">>> Removing sequencing adapters with fastp..."
            adapter_removed_path="${dataset_path}temp/step_01_adapter_removed/"
            mkdir -p "$adapter_removed_path"

            for fq in "$ori_fastq_path"*.fastq*; do
                [[ -f "$fq" ]] || continue
                fastp -i "$fq" \
                      -o "${adapter_removed_path}$(basename "$fq")" \
                      --disable_quality_filtering \
                      --disable_length_filtering \
                      -w "$cpu" \
                      -j "${adapter_removed_path}fastp.json" \
                      -h "${adapter_removed_path}fastp.html"
                echo "  ✓ Adapter removal: $(basename "$fq")"
            done
            echo "✓ Adapter removal completed"

            # ── Step C: Entropy-based primer detection & trimming (same as Illumina) ──
            echo ">>> Detecting and removing primers (entropy method)..."
            fastp_path="${dataset_path}temp/step_02_fastp/"
            mkdir -p "$fastp_path"

            python3 "${SCRIPTS}/entropy_primer_detect.py" \
                -i "$adapter_removed_path" \
                -o "$fastp_path" || {
                echo "  ✗ Entropy primer detection failed"
                exit 1
            }

            echo "✓ Entropy primer detection/removal completed"

            echo ">>> Cleaning up intermediate files..."
            rm -rf "$ori_fastq_path"
            rm -rf "$adapter_removed_path"
            echo "✓ Removed ori_fastq/ and adapter_removed/ to save space"

            # ── Step D: QIIME2 Import → Length filter → Dedup → Chimera → OTU 97% ──
            fastq_path="$fastp_path"
            export fastq_path
            Amplicon_Common_MakeManifestFileForQiime2
            Amplicon_Common_ImportFastqToQiime2
            Amplicon_LS454_QualityControlForQZA
            Amplicon_LS454_Deduplication
            Amplicon_LS454_ChimerasRemoval
            Amplicon_LS454_ClusterDenovo
            Amplicon_Common_FinalFilesCleaning

        elif [[ "$platform" == "ION_TORRENT" ]]; then
            # ── Step A: Download with normal prefetch + fasterq-dump (same as 454) ──
            echo ">>> Downloading SRA data (Ion Torrent)..."
            Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}"

            # Count raw reads before any processing
            Common_CountRawReads "$dataset_path" "$sra_file_name"

            # Ion Torrent treated as single-end (same as 454)
            sequence_type="single"
            echo "Sequence type: ${sequence_type^^}"
            export sequence_type

            # ── Step B: Remove sequencing adapters with fastp (SE mode) ──
            echo ">>> Removing sequencing adapters with fastp..."
            adapter_removed_path="${dataset_path}temp/step_01_adapter_removed/"
            mkdir -p "$adapter_removed_path"

            for fq in "$ori_fastq_path"*.fastq*; do
                [[ -f "$fq" ]] || continue
                fastp -i "$fq" \
                      -o "${adapter_removed_path}$(basename "$fq")" \
                      --disable_quality_filtering \
                      --disable_length_filtering \
                      -w "$cpu" \
                      -j "${adapter_removed_path}fastp.json" \
                      -h "${adapter_removed_path}fastp.html"
                echo "  ✓ Adapter removal: $(basename "$fq")"
            done
            echo "✓ Adapter removal completed"

            # ── Step C: Entropy-based primer detection & trimming (same as Illumina) ──
            echo ">>> Detecting and removing primers (entropy method)..."
            fastp_path="${dataset_path}temp/step_02_fastp/"
            mkdir -p "$fastp_path"

            python3 "${SCRIPTS}/entropy_primer_detect.py" \
                -i "$adapter_removed_path" \
                -o "$fastp_path" || {
                echo "  ✗ Entropy primer detection failed"
                exit 1
            }

            echo "✓ Entropy primer detection/removal completed"

            # ── Step C2: Trim first 15bp of biological sequence (post-primer) ──
            # Ion Torrent reads have low-quality bases in the first ~15bp after
            # the primer due to signal instability at sequencing start.
            echo ">>> Trimming first 15bp of biological sequence (post-primer)..."
            trimmed_path="${dataset_path}temp/step_02b_trimmed/"
            mkdir -p "$trimmed_path"

            for fq in "$fastp_path"*.fastq*; do
                [[ -f "$fq" ]] || continue
                fastp -i "$fq" \
                      -o "${trimmed_path}$(basename "$fq")" \
                      --trim_front1 15 \
                      --disable_adapter_trimming \
                      --disable_quality_filtering \
                      --disable_length_filtering \
                      -w "$cpu" \
                      -j "${trimmed_path}fastp_trim.json" \
                      -h "${trimmed_path}fastp_trim.html"
                echo "  ✓ Front 15bp trimmed: $(basename "$fq")"
            done
            echo "✓ Post-primer 15bp trim completed"

            echo ">>> Cleaning up intermediate files..."
            rm -rf "$ori_fastq_path"
            rm -rf "$adapter_removed_path"
            rm -rf "$fastp_path"
            echo "✓ Removed ori_fastq/, adapter_removed/, and fastp/ to save space"

            # ── Step D: QIIME2 Import → Length filter → Dedup → Chimera → OTU 97% ──
            fastq_path="$trimmed_path"
            export fastq_path
            Amplicon_Common_MakeManifestFileForQiime2
            Amplicon_Common_ImportFastqToQiime2
            Amplicon_IonTorrent_QualityControlForQZA
            Amplicon_LS454_Deduplication
            Amplicon_LS454_ChimerasRemoval
            Amplicon_LS454_ClusterDenovo
            Amplicon_Common_FinalFilesCleaning

        elif [[ "$platform" == "OXFORD_NANOPORE" ]]; then
            echo "⚠️ OXFORD_NANOPORE platform is not supported for amplicon analysis."
            echo "   Skipping dataset: $dataset_ID"
            echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - SKIPPED - Platform: $platform (not supported)" >> "$failed_log"
            continue

        elif [[ "$platform" == "PACBIO_SMRT" ]]; then
            # ── Step A: Download ──
            echo ">>> Downloading SRA data..."
            Common_SRADownloadToFastq_MultiSource -d "$dataset_path" -a "${sra_file_name}"

            # Count raw reads before any processing
            Common_CountRawReads "$dataset_path" "$sra_file_name"

            # PacBio is typically SE
            sequence_type="single"
            echo "Sequence type: ${sequence_type^^}"
            export sequence_type

            # ── Step B: Remove sequencing adapters with fastp ──
            echo ">>> Removing sequencing adapters with fastp..."
            adapter_removed_path="${dataset_path}temp/step_01_adapter_removed/"
            mkdir -p "$adapter_removed_path"

            # PacBio is typically SE
            for fq in "$ori_fastq_path"*.fastq*; do
                [[ -f "$fq" ]] || continue
                fastp -i "$fq" \
                      -o "${adapter_removed_path}$(basename "$fq")" \
                      --disable_quality_filtering \
                      --disable_length_filtering \
                      -w "$cpu" \
                      -j "${adapter_removed_path}fastp.json" \
                      -h "${adapter_removed_path}fastp.html"
                echo "  ✓ Adapter removal: $(basename "$fq")"
            done
            echo "✓ Adapter removal completed"

            # ── Step B2: Entropy-based primer DETECTION (no trimming) ──
            # For PacBio CCS, DADA2 denoise-ccs needs the primer sequence
            # (--p-front) to re-orient reads (CCS reads are in random
            # forward/reverse-complement orientations). So we only detect
            # the primer here and let DADA2 handle both orientation and
            # primer removal.
            echo ">>> Detecting primers (entropy method, no trimming)..."
            primer_detect_path="${dataset_path}temp/step_02_primer_detect/"
            mkdir -p "$primer_detect_path"

            python3 "${SCRIPTS}/entropy_primer_detect.py" \
                -i "$adapter_removed_path" \
                -o "$primer_detect_path" || {
                echo "  ✗ Entropy primer detection failed"
                exit 1
            }

            echo "✓ Entropy primer detection completed"

            # Read detected primer consensus from JSON
            primer_info_json="${primer_detect_path}primer_info.json"
            if [[ -f "$primer_info_json" ]]; then
                detected_primer=$(python3 -c "
import json, sys
with open('${primer_info_json}') as f:
    info = json.load(f)
fp = info.get('forward_primer', {})
if fp.get('detected', False) and fp.get('consensus', ''):
    print(fp['consensus'])
else:
    print('')
")
                echo "  Detected forward primer: ${detected_primer:-none}"
            else
                echo "  ⚠ primer_info.json not found"
                detected_primer=""
            fi

            # Clean up detection output (we only needed the JSON)
            echo ">>> Cleaning up intermediate files..."
            rm -rf "$ori_fastq_path"
            rm -rf "$primer_detect_path"
            echo "✓ Removed ori_fastq/ and primer_detect/ to save space"

            # ── Step C: QIIME2 Import → QC → DADA2 denoise-ccs ──
            # Import adapter-removed reads (still containing primers for DADA2)
            fastq_path="$adapter_removed_path"
            export fastq_path
            export detected_primer
            Amplicon_Common_MakeManifestFileForQiime2
            Amplicon_Common_ImportFastqToQiime2
            Amplicon_Pacbio_QualityControlForQZA
            Amplicon_Pacbio_DenosingDada2
            Amplicon_Common_FinalFilesCleaning

        else
            echo "❌ Unknown platform: $platform"
            exit 1
        fi

        # Generate summary CSV entry
        echo ">>> Generating summary for $dataset_ID..."
        python "${SCRIPTS}/py_16s.py" append_summary \
            --dataset_id "$dataset_ID" \
            --sra_file "${dataset_path}${sra_file_name}" \
            --raw_counts "${dataset_path}${dataset_ID}_raw_read_counts.tsv" \
            --final_table "${dataset_path}${dataset_ID}-final-table.qza" \
            --output_csv "$summary_csv"

        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - SUCCESS - Platform: $platform" >> "$success_log"

    } || {
        echo "❌ Pipeline failed for $dataset_ID"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - FAILED - Platform: $platform" >> "$failed_log"
    }
done

echo "========================================="
echo "ALL DONE. Check logs for summary."
echo "========================================="