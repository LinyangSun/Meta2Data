#!/bin/bash
#
# AmpliconFunction.sh - Unified Amplicon Processing Functions
#
# This file consolidates all amplicon sequencing processing functions
# for different platforms (454, Illumina, PacBio, ONT) and common utilities.
#
# Functions are organized by platform and processing stage:
#   - Common Functions: Shared utilities across all platforms
#   - 454 Functions: Roche 454 pyrosequencing pipeline
#   - Illumina Functions: Illumina short-read pipeline
#   - PacBio Functions: PacBio HiFi/CCS long-read pipeline
#   - ONT Functions: Oxford Nanopore long-read pipeline (TODO)
#
# Usage: Source this file in your pipeline script
#   source scripts/AmpliconFunction.sh
#
# Author: Meta2Data Project
# Date: $(date +%Y-%m-%d)
#

set -e


################################################################################
#                         COMMON FUNCTIONS                                     #
################################################################################
# Functions shared across all sequencing platforms

Verify_Fastq_Integrity() {
    # Verify that a FASTQ file is valid (line count divisible by 4, has content).
    # Args: $1 = path to FASTQ file
    # Returns: 0 if valid, 1 if invalid
    local fq_file="$1"

    if [[ ! -f "$fq_file" ]]; then
        echo "  Integrity check: file not found: $fq_file" >&2
        return 1
    fi

    local file_size
    file_size=$(stat -c%s "$fq_file" 2>/dev/null || stat -f%z "$fq_file" 2>/dev/null)
    if [[ "$file_size" -eq 0 ]]; then
        echo "  Integrity check: empty file: $(basename "$fq_file")" >&2
        return 1
    fi

    local line_count
    if [[ "$fq_file" == *.gz ]]; then
        line_count=$(zcat "$fq_file" 2>/dev/null | wc -l)
    else
        line_count=$(wc -l < "$fq_file")
    fi

    if [[ "$line_count" -eq 0 ]]; then
        echo "  Integrity check: no lines in file: $(basename "$fq_file")" >&2
        return 1
    fi

    if (( line_count % 4 != 0 )); then
        echo "  Integrity check: line count ($line_count) not divisible by 4: $(basename "$fq_file")" >&2
        return 1
    fi

    return 0
}

Download_From_ENA() {
    # Download FASTQ file(s) directly from ENA (European Nucleotide Archive).
    # ENA provides pre-converted FASTQ files, serving as a reliable fallback
    # when NCBI prefetch fails.
    #
    # Args:
    #   $1 = SRR/ERR/DRR accession
    #   $2 = target directory for FASTQ files
    #   $3 = rename prefix for output files
    #
    # Returns: 0 on success, 1 on failure
    local srr="$1"
    local target_dir="$2"
    local rename_prefix="$3"
    local max_retries=5

    echo "  [ENA] Attempting ENA download for $srr..."

    # ENA URL structure: first 6 chars of accession, then full accession
    # e.g., SRR123456 -> /vol1/fastq/SRR123/006/SRR1234566/
    local acc_prefix="${srr:0:6}"
    local acc_len=${#srr}

    # Build ENA directory path based on accession length
    local ena_dir
    if (( acc_len <= 9 )); then
        ena_dir="https://ftp.sra.ebi.ac.uk/vol1/fastq/${acc_prefix}/${srr}"
    elif (( acc_len == 10 )); then
        ena_dir="https://ftp.sra.ebi.ac.uk/vol1/fastq/${acc_prefix}/00${srr: -1}/${srr}"
    elif (( acc_len == 11 )); then
        ena_dir="https://ftp.sra.ebi.ac.uk/vol1/fastq/${acc_prefix}/0${srr: -2}/${srr}"
    else
        ena_dir="https://ftp.sra.ebi.ac.uk/vol1/fastq/${acc_prefix}/${srr: -3}/${srr}"
    fi

    echo "  [ENA] Base URL: $ena_dir"

    # Try PE first (_1.fastq.gz and _2.fastq.gz), then SE (.fastq.gz)
    local -a try_files=("${srr}_1.fastq.gz" "${srr}_2.fastq.gz" "${srr}.fastq.gz")
    local downloaded_any=false

    for fname in "${try_files[@]}"; do
        local url="${ena_dir}/${fname}"
        local out_file="${target_dir}/${fname}"
        local success=false

        echo "  [ENA] Trying: $fname"
        for ((attempt=1; attempt<=max_retries; attempt++)); do
            local wget_err=""
            if wget_err=$(wget --timeout=60 --tries=1 "$url" -O "$out_file" 2>&1); then
                # Verify the file is a valid gzip
                if gzip -t "$out_file" 2>/dev/null; then
                    local dl_size
                    dl_size=$(stat -c%s "$out_file" 2>/dev/null || stat -f%z "$out_file" 2>/dev/null || echo "?")
                    echo "  [ENA] ✓ Downloaded $fname (${dl_size} bytes)"
                    success=true
                    break
                else
                    echo "  [ENA] Corrupt download for $fname (attempt $attempt/$max_retries)"
                    rm -f "$out_file"
                fi
            else
                # Extract HTTP status code from wget output
                local http_code
                http_code=$(echo "$wget_err" | grep -oP 'HTTP request sent.*\K[0-9]{3}' | tail -1)
                if [[ -n "$http_code" ]]; then
                    echo "  [ENA] ✗ $fname: HTTP $http_code (attempt $attempt/$max_retries)"
                else
                    local wget_reason
                    wget_reason=$(echo "$wget_err" | tail -1)
                    echo "  [ENA] ✗ $fname: $wget_reason (attempt $attempt/$max_retries)"
                fi
                rm -f "$out_file"
            fi
            # Exponential backoff: 5, 10, 20, 40, 60 (capped)
            local wait=$(( 5 * (1 << (attempt - 1)) ))
            (( wait > 60 )) && wait=60
            [[ $attempt -lt $max_retries ]] && sleep "$wait"
        done

        if [[ "$success" == true ]]; then
            # Rename to match expected naming convention
            local base_filename="${fname/${srr}/}"
            mv "$out_file" "${target_dir}/${rename_prefix}${base_filename}"
            echo "  [ENA] → Renamed to: ${rename_prefix}${base_filename}"
            downloaded_any=true
        else
            rm -f "$out_file" 2>/dev/null
        fi
    done

    if [[ "$downloaded_any" == true ]]; then
        # If we got both _1 and _2, remove the SE file (.fastq.gz) if it also downloaded
        if [[ -f "${target_dir}/${rename_prefix}_1.fastq.gz" ]] && \
           [[ -f "${target_dir}/${rename_prefix}_2.fastq.gz" ]]; then
            rm -f "${target_dir}/${rename_prefix}.fastq.gz" 2>/dev/null
        fi
        return 0
    fi

    echo "  [ENA] ✗ ENA download failed for $srr — no FASTQ files available" >&2
    return 1
}

Download_CRR() {
    local crr=$1
    local target_dir=$2
    local rename_prefix=$3
    local max_retries=5

    echo "[CNCB] Processing $crr"

    # 1. Map CRR to parent CRA ID (with retry)
    local cra=""
    for ((attempt=1; attempt<=3; attempt++)); do
        cra=$(wget -qO- --timeout=30 --user-agent="Mozilla/5.0" "https://ngdc.cncb.ac.cn/gsa/search?searchTerm=${crr}" | grep -v "example" | grep -oe "CRA[0-9]\+" | uniq | head -n 1)
        [[ -n "$cra" ]] && break
        local wait=$(( 5 * attempt ))
        echo "  CRA lookup retry $attempt/3 (wait ${wait}s)..." >&2
        sleep "$wait"
    done

    if [[ -z "$cra" ]]; then
        echo "Error: Could not map $crr to a CRA project." >&2
        return 1
    fi

    # 2. Fetch project MD5 list and path prefix (with retry)
    local md5_file="${cra}_md5.txt"
    local path_prefix=""
    for ((attempt=1; attempt<=3; attempt++)); do
        path_prefix=$(wget -qO- --timeout=30 --post-data="searchTerm=${cra}&totalDatas=1&downLoadCount=1" --user-agent="Mozilla/5.0" "https://ngdc.cncb.ac.cn/gsa/search/getRunInfoByCra" | grep -oe "gsa[0-9]\+/${cra}\|gsa/${cra}" | uniq | head -n 1)
        [[ -n "$path_prefix" ]] && break
        sleep $(( 5 * attempt ))
    done

    if [[ ! -s "$md5_file" ]]; then
        wget -q --timeout=30 "https://download.cncb.ac.cn/${path_prefix}/md5sum.txt" -O "$md5_file"
    fi

    # 3. Identify files belonging to the CRR accession
    local -a files
    while IFS= read -r filepath; do
        files+=("$(basename "$filepath")")
    done < <(grep "$crr" "$md5_file" | awk '{print $2}')

    for filename in "${files[@]}"; do
        local expected_md5=$(grep "$filename" "$md5_file" | awk '{print tolower($1)}')
        local success=false

        for ((i=1; i<=max_retries; i++)); do
            local wget_url="https://download.cncb.ac.cn/${path_prefix}/${crr}/${filename}"
            if wget -q --timeout=120 --user-agent="Mozilla/5.0" "$wget_url" -O "$filename"; then
                if [[ -f "$filename" ]]; then
                    local current_md5=$(md5sum "$filename" | awk '{print $1}')
                    if [[ "$current_md5" == "$expected_md5" ]]; then
                        # Standardize filename: _r1/_R1 → _1, _r2/_R2 → _2
                        local new_filename=$(echo "$filename" | sed -E 's/_[rR]1([._])/_1\1/; s/_[rR]2([._])/_2\1/')
                        # Strip CRR ID from filename (keep the separator)
                        local base_filename=$(echo "$new_filename" | sed "s/${crr}//")
                        mv "$filename" "${target_dir}/${rename_prefix}${base_filename}"
                        success=true
                        echo "  ✓ Downloaded: $filename → ${rename_prefix}${base_filename}"
                        break
                    else
                        echo "  MD5 mismatch (attempt $i/$max_retries), retrying..."
                        rm -f "$filename"
                    fi
                fi
            else
                echo "  Download failed (attempt $i/$max_retries)..." >&2
                rm -f "$filename"
            fi
            # Exponential backoff: 5, 10, 20, 40, 60 (capped)
            local wait=$(( 5 * (1 << (i - 1)) ))
            (( wait > 60 )) && wait=60
            [[ $i -lt $max_retries ]] && sleep "$wait"
        done

        if [[ "$success" = false ]]; then
            echo "Error: Failed to download $filename after $max_retries attempts." >&2
            return 1
        fi
    done
    return 0
}


Common_SRADownloadToFastq_MultiSource() {
    local dir_path="" acc_file=""
    OPTIND=1
    while getopts ":d:a:" opt; do
        case $opt in
            d) dir_path=$OPTARG ;;
            a) acc_file=$OPTARG ;;
            ?) echo "Unknown Parameter: -$OPTARG" >&2; return 1 ;;
        esac
    done

    if [[ -z "$dir_path" || -z "$acc_file" ]]; then
        echo "Usage: Common_SRADownloadToFastq_MultiSource -d <dir> -a <accession_tsv>" >&2
        return 1
    fi

    # Pre-scan accession file to determine required tools
    local has_ncbi_accessions=false
    local has_cncb_accessions=false
    local base_dir="${dir_path%/}"

    while IFS=$'\t' read -r srr _; do
        [[ -z "$srr" ]] && continue
        if [[ "$srr" =~ ^CRR ]]; then
            has_cncb_accessions=true
        elif [[ "$srr" =~ ^[EDS]RR ]]; then
            has_ncbi_accessions=true
        fi
    done < "${base_dir}/${acc_file}"

    # Check dependencies based on what we'll download
    if [[ "$has_ncbi_accessions" == true ]]; then
        command -v prefetch >/dev/null 2>&1 || { echo "Error: 'prefetch' not found (required for NCBI SRR/ERR/DRR downloads)." >&2; return 1; }
        command -v fasterq-dump >/dev/null 2>&1 || { echo "Error: 'fasterq-dump' not found (required for NCBI SRR/ERR/DRR downloads)." >&2; return 1; }
    fi

    if [[ "$has_cncb_accessions" == true ]]; then
        command -v wget >/dev/null 2>&1 || { echo "Error: 'wget' not found (required for CNCB CRR downloads)." >&2; return 1; }
    fi

    # Quick connectivity check — fail fast instead of waiting through retries
    if [[ "$has_ncbi_accessions" == true ]]; then
        echo "  Checking NCBI connectivity..."
        if ! wget -q --spider --timeout=10 "https://ftp.ncbi.nlm.nih.gov/" 2>/dev/null; then
            echo "  ⚠️  NCBI unreachable — will rely on ENA fallback"
            echo "  Checking ENA connectivity..."
            if ! wget -q --spider --timeout=10 "https://ftp.sra.ebi.ac.uk/" 2>/dev/null; then
                echo "  ✗ Both NCBI and ENA are unreachable. Check your network connection." >&2
                return 1
            fi
            echo "  ✓ ENA reachable"
        else
            echo "  ✓ NCBI reachable"
        fi
    fi

    # base_dir already declared above
    local fastq_path="${base_dir}/ori_fastq"
    local temp_dl_path="${base_dir}/temp1"
    mkdir -p "$fastq_path" "$temp_dl_path"

    # Phase 1: Data Acquisition
    while IFS=$'\t' read -r srr rename _; do
        [[ -z "$srr" || -z "$rename" ]] && continue

        # Route by ID type (CRR vs SRR/ERR/DRR)
        if [[ "$srr" =~ ^CRR ]]; then
            # CRR files go directly to ori_fastq/ with rename prefix
            Download_CRR "$srr" "$fastq_path" "$rename"

        elif [[ "$srr" =~ ^[EDS]RR ]]; then
            echo "[NCBI] Processing $srr"
            local max_attempts=3
            local attempt=0
            local dl_ok=false
            local prefetch_err=""
            while (( attempt < max_attempts )); do
                (( ++attempt ))
                echo "  [prefetch] Attempt $attempt/$max_attempts for $srr..."
                prefetch_err=$(prefetch -O "$temp_dl_path" "$srr" 2>&1) && { echo "  [prefetch] ✓ Success"; dl_ok=true; break; }
                echo "  [prefetch] ✗ Failed (attempt $attempt/$max_attempts)"
                echo "  [prefetch] Error: ${prefetch_err##*$'\n'}"
                if (( attempt < max_attempts )); then
                    # Exponential backoff: 10, 20, 40, 60, 60, ... (capped at 60s)
                    local wait=$(( 10 * (1 << (attempt - 1)) ))
                    (( wait > 60 )) && wait=60
                    echo "  [prefetch] Waiting ${wait}s before retry..."
                    sleep "$wait"
                fi
            done

            if [[ "$dl_ok" == true ]]; then
                # Extract downloaded SRA files to temp
                local sra_files=()
                while IFS= read -r -d '' file; do
                    sra_files+=("$file")
                done < <(find "$temp_dl_path/$srr" -type f \( -name "*.sra" -o -name "*.fastq*" \) -print0 2>/dev/null)

                echo "  [prefetch] Found ${#sra_files[@]} file(s) to convert"

                # Process each file
                local convert_ok=true
                for file in "${sra_files[@]}"; do
                    local basename_file=$(basename "$file")
                    local ext="${basename_file##*.}"
                    local file_size
                    file_size=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null || echo "?")
                    echo "  [convert] Processing: $basename_file (${file_size} bytes)"

                    if [[ "$ext" == "sra" ]]; then
                        # Convert SRA to FASTQ with rename prefix
                        local temp_outdir="${temp_dl_path}/fastq_${srr}"
                        mkdir -p "$temp_outdir"

                        echo "  [fasterq-dump] Converting $basename_file to FASTQ..."
                        local fqd_err=""
                        if fqd_err=$(fasterq-dump --split-3 "$file" --outdir "$temp_outdir" 2>&1); then
                            # List converted files
                            local fq_count
                            fq_count=$(find "$temp_outdir" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | wc -l)
                            echo "  [fasterq-dump] ✓ Converted to $fq_count FASTQ file(s)"
                            # Rename and move converted FASTQ files
                            find "$temp_outdir" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) -print0 |
                            while IFS= read -r -d '' fq_file; do
                                local fq_basename=$(basename "$fq_file")
                                # Standardize to _1/_2 format first
                                local normalized=$(echo "$fq_basename" | sed -E 's/_[rR]?1([._])/_1\1/; s/_[rR]?2([._])/_2\1/')
                                # Strip SRR/ERR/DRR ID from filename (keep separator)
                                local base_filename=$(echo "$normalized" | sed "s/${srr}//")
                                echo "  [fasterq-dump] → ${rename}${base_filename}"
                                mv "$fq_file" "${fastq_path}/${rename}${base_filename}"
                            done
                            rm -rf "$temp_outdir"
                        else
                            echo "  [fasterq-dump] ✗ Failed for $basename_file"
                            echo "  [fasterq-dump] Error: $fqd_err"
                            convert_ok=false
                        fi
                        rm -f "$file"

                    elif [[ "$ext" == "fastq" || "$basename_file" =~ \.fastq\.gz$ ]]; then
                        # Already FASTQ, standardize and strip SRR/ERR/DRR ID
                        # Standardize to _1/_2 format first
                        local normalized=$(echo "$basename_file" | sed -E 's/_[rR]?1([._])/_1\1/; s/_[rR]?2([._])/_2\1/')
                        # Strip SRR/ERR/DRR ID from filename (keep separator)
                        local base_filename=$(echo "$normalized" | sed "s/${srr}//")
                        echo "  [convert] Already FASTQ → ${rename}${base_filename}"
                        mv "$file" "${fastq_path}/${rename}${base_filename}"
                    fi
                done

                rm -rf "${temp_dl_path:?}/$srr"

                # If fasterq-dump conversion failed, fall through to ENA
                if [[ "$convert_ok" == false ]]; then
                    dl_ok=false
                fi
            fi

            # Fallback: download directly from ENA if prefetch/fasterq-dump failed
            if [[ "$dl_ok" != true ]]; then
                echo "  NCBI download failed for $srr after $max_attempts attempts — trying ENA fallback..."
                if Download_From_ENA "$srr" "$fastq_path" "$rename"; then
                    echo "  ✓ ENA fallback succeeded for $srr"
                else
                    echo "  ✗ All download sources failed for $srr (NCBI + ENA)" >&2
                    return 1
                fi
            fi

            # Verify FASTQ integrity for downloaded files
            echo "  [verify] Checking FASTQ integrity for $srr..."
            local verify_failed=false
            for fq in "${fastq_path}/${rename}"*.fastq*; do
                [[ -f "$fq" ]] || continue
                if ! Verify_Fastq_Integrity "$fq"; then
                    echo "  [verify] ✗ Invalid: $(basename "$fq") — removing"
                    rm -f "$fq"
                    verify_failed=true
                else
                    local fq_size
                    fq_size=$(stat -c%s "$fq" 2>/dev/null || stat -f%z "$fq" 2>/dev/null || echo "?")
                    echo "  [verify] ✓ OK: $(basename "$fq") (${fq_size} bytes)"
                fi
            done
            if [[ "$verify_failed" == true ]]; then
                # Check if any valid files remain
                local remaining_files
                remaining_files=$(find "$fastq_path" -name "${rename}*.fastq*" -type f 2>/dev/null | wc -l)
                if [[ "$remaining_files" -eq 0 ]]; then
                    echo "  [verify] ✗ No valid FASTQ files remaining for $srr" >&2
                    return 1
                fi
                echo "  [verify] $remaining_files valid file(s) remaining after cleanup"
            fi
        else
            echo "Warning: Unknown Accession format: $srr" >&2
        fi
    done < "${base_dir}/${acc_file}"
    rm -rf "$temp_dl_path"

    # Orphan file cleanup: fasterq-dump --split-3 may produce 3 files per sample
    # (prefix.fastq + prefix_1.fastq + prefix_2.fastq). The unpaired orphan
    # reads (prefix.fastq) must be removed when paired files exist, otherwise
    # downstream tools will misinterpret them as additional SE samples.
    local orphan_count=0
    for r1 in "${fastq_path}/"*_1.fastq*; do
        [[ -f "$r1" ]] || continue
        local prefix="${r1%_1.fastq*}"
        # Verify R2 also exists
        local has_r2=false
        for r2 in "${prefix}_2.fastq"*; do
            [[ -f "$r2" ]] && has_r2=true && break
        done
        $has_r2 || continue
        # Both R1 and R2 present — remove orphan (prefix.fastq / prefix.fastq.gz)
        for orphan in "${prefix}.fastq" "${prefix}.fastq.gz"; do
            if [[ -f "$orphan" ]]; then
                echo "  Removing orphan file: $(basename "$orphan")"
                rm -f "$orphan"
                orphan_count=$((orphan_count + 1))
            fi
        done
    done
    if [[ $orphan_count -gt 0 ]]; then
        echo "  Cleaned up $orphan_count orphan file(s)"
    fi
    return 0

}

Amplicon_Common_MakeManifestFileForQiime2() {
    cd "${dataset_path%/}" || { echo "ERROR: Cannot access dataset path: $dataset_path"; exit 1; }
    local temp_file_path="${dataset_path%/}/temp/temp_file"
    mkdir -p "$temp_file_path"
    local _dp="${dataset_path%/}"
    dataset_name="${_dp##*/}"
    find "$fastq_path" -type f -name "*.fastq*" > "${temp_file_path}/${dataset_name}-file.txt"
    if [ "$sequence_type" = "single" ]; then
        python "${SCRIPTS}/py_16s.py" mk_manifest_SE --FilePath "${temp_file_path}/${dataset_name}-file.txt"
    else
        python "${SCRIPTS}/py_16s.py" mk_manifest_PE --FilePath "${temp_file_path}/${dataset_name}-file.txt"
    fi
}
Amplicon_Common_ImportFastqToQiime2() {
    set -u
    cd "$dataset_path" || { echo "❌ Cannot access dataset path: $dataset_path"; exit 1; }

    local temp_path="${dataset_path%/}/temp/"
    local temp_file_path="${temp_path}temp_file/"
    local qza_path="${temp_path}step_03_qza_import/"
    mkdir -p "$temp_file_path" "$qza_path"

    local dataset_name="${dataset_path%/}"
    dataset_name="${dataset_name##*/}"

    local paired_manifest="${temp_file_path}${dataset_name}_manifest.tsv"

    # --- SINGLE-END ---
    if [ "${sequence_type:-paired}" = "single" ]; then
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path "$paired_manifest" \
            --output-path "${qza_path}${dataset_name}.qza" \
            --input-format SingleEndFastqManifestPhred33V2
        return
    fi

    # --- PAIRED-END ---
    # Import as paired-end; DADA2 denoise-paired handles merging internally
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$paired_manifest" \
        --output-path "${qza_path}${dataset_name}.qza" \
        --input-format PairedEndFastqManifestPhred33V2
}
Amplicon_Common_FinalFilesCleaning() {
    dataset_path="${dataset_path%/}/"
    local quality_filter_path="${dataset_path%/}/temp/step_04_qza_import_QualityFilter/"
    local denoising_path="${dataset_path%/}/temp/step_05_denoise/"
    local qf_vis_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/"
    local qf_view_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/qf_trim_pos/"
    local qc_vis="${dataset_path%/}/temp/temp_file/qc_vis/"
    local denoising_vis="${dataset_path%/}/temp/temp_file/denoising_vis/"
    
    cd "$dataset_path" || { echo "Error: dataset_path not found"; return 1; }
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    quality_filter_path="${dataset_path%/}/temp/step_04_qza_import_QualityFilter/"
    denoising_path="${dataset_path%/}/temp/step_05_denoise/"
    qf_trim_pos_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/qf_trim_pos/"
    qc_vis="${dataset_path%/}/temp/temp_file/qc_vis/"
    denoising_vis="${dataset_path%/}/temp/temp_file/denoising_vis/"
    mkdir -p "$qc_vis" "$denoising_vis"
    
    # Check for denoising output
    if [ -d "$denoising_path" ] && [ -f "${denoising_path%/}/${dataset_name}-table-denoising.qza" ]; then
        cp "${denoising_path%/}/${dataset_name}-rep-seqs-denoising.qza" "${dataset_path%/}/${dataset_name}-final-rep-seqs.qza"
        cp "${denoising_path%/}/${dataset_name}-table-denoising.qza" "${dataset_path%/}/${dataset_name}-final-table.qza"
        
        # Copy trim position files if they exist
        if [ -f "${qf_trim_pos_path%/}/Trim_position.txt" ]; then
            cp "${qf_trim_pos_path%/}/Trim_position.txt" "${dataset_path%/}/${dataset_name}-TrimPosition.txt"
        fi
        if [ -f "${qf_trim_pos_path%/}/forward-seven-number-summaries.tsv" ]; then
            cp "${qf_trim_pos_path%/}/forward-seven-number-summaries.tsv" "${dataset_path%/}/${dataset_name}-TrimPositionOriginalFile.tsv"
        fi
        
        # Extract stats
        if [ -f "${quality_filter_path%/}/${dataset_name}_filter-stats.qza" ]; then
            unzip -q "${quality_filter_path%/}/${dataset_name}_filter-stats.qza" -d "$qc_vis"
            find "$qc_vis" -type f -name 'stats.csv' -exec cp {} "${dataset_path%/}/${dataset_name}-QCStats.csv" \;
        fi
        
        if [ -f "${denoising_path%/}/${dataset_name}-denoising-stats.qza" ]; then
            unzip -q "${denoising_path%/}/${dataset_name}-denoising-stats.qza" -d "$denoising_vis"
            find "$denoising_vis" -type f -name 'stats.csv' -exec cp {} "${dataset_path%/}/${dataset_name}-DenoisingStats.csv" \;
        fi
        
        # Remove temporary directories but keep dataset_path itself
        find "$dataset_path" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
        
        # Remove unwanted logs
        rm -f "${dataset_path%/}/"{denoising.log,fastp.html,fastp.json}

        rm -rf "${dataset_path%/}/ori_fastq" 2>/dev/null || true
        rm -rf "${dataset_path%/}/working_fastq" 2>/dev/null || true
        
        return 0
        
    # Check for vsearch output (454 pipeline)
    elif [ -f "${dataset_path%/}/${dataset_name}-table-vsearch.qza" ]; then
        mv "${dataset_path%/}/${dataset_name}-table-vsearch.qza" "${dataset_path%/}/${dataset_name}-final-table.qza"
        mv "${dataset_path%/}/${dataset_name}-rep-seqs-vsearch.qza" "${dataset_path%/}/${dataset_name}-final-rep-seqs.qza"
        find "$dataset_path" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
        
        return 0
        
    else
        echo "❌ ERROR: The analysis failed! The final denoising output does not exist."
        echo "   Expected files:"
        echo "     - ${denoising_path%/}/${dataset_name}-table-denoising.qza "
        echo ""
        echo "   Check previous steps for errors."
        
        # Don't clean up temp files so we can debug
        return 1
    fi
}



Common_CountRawReads() {
    # Count raw reads per sample from ori_fastq and save to TSV
    # Args: $1 = dataset_path, $2 = sra_file_name
    local base_dir="${1%/}"
    local acc_file="$2"
    local fastq_path="${base_dir}/ori_fastq"
    local dataset_name="${base_dir##*/}"
    local raw_counts_file="${base_dir}/${dataset_name}_raw_read_counts.tsv"
    : > "$raw_counts_file"

    while IFS=$'\t' read -r srr rename _; do
        [[ -z "$srr" || -z "$rename" ]] && continue
        local total_lines=0
        for fq in "${fastq_path}/${rename}"*.fastq*; do
            [[ -f "$fq" ]] || continue
            local lines
            if [[ "$fq" == *.gz ]]; then
                lines=$(zcat "$fq" | wc -l)
            else
                lines=$(wc -l < "$fq")
            fi
            total_lines=$((total_lines + lines))
        done
        local total_reads=$((total_lines / 4))
        printf '%s\t%s\t%d\n' "$srr" "$rename" "$total_reads" >> "$raw_counts_file"
    done < "${base_dir}/${acc_file}"
}

################################################################################
#                       ILLUMINA PLATFORM FUNCTIONS                            #
################################################################################
# Functions for Illumina short-read sequencing data processing



Amplicon_Illumina_QualityControlForQZA(){
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local qza_path="${dataset_path%/}/temp/step_03_qza_import/"
    local quality_filter_path="${dataset_path%/}/temp/step_04_qza_import_QualityFilter/"
    local qf_vis_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/"
    mkdir -p "$quality_filter_path" "$qf_vis_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    qiime quality-filter q-score \
        --i-demux "${qza_path%/}/${dataset_name}.qza" \
        --p-min-quality 20 \
        --o-filtered-sequences "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --o-filter-stats "${quality_filter_path%/}/${dataset_name}_filter-stats.qza" \
        --verbose
}

Amplicon_Illumina_DenosingDada2() {
    # Parse flags: -s <start>, -e <end>
    local start_in="" end_in="" opt
    OPTIND=1
    while getopts ":s:e:" opt; do
        case "$opt" in
            s) start_in="$OPTARG" ;;
            e) end_in="$OPTARG" ;;
            \?) echo "Invalid option: -$OPTARG" >&2; return 2 ;;
            :)  echo "Option -$OPTARG requires an argument." >&2; return 2 ;;
        esac
    done
    shift $((OPTIND-1))
    set -euo pipefail
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    # Paths
    local base="${dataset_path%/}"
    local quality_filter_path="${base%/}/temp/step_04_qza_import_QualityFilter/"
    local denoising_path="${base%/}/temp/step_05_denoise/"
    local qf_vis_path="${base%/}/temp/temp_file/QualityFilter_vis/"
    local qf_view_path="${base%/}/temp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${base%/}/temp/temp_file/QualityFilter_vis/qf_trim_pos/"
    mkdir -p "$denoising_path" "$qf_view_path" "$qf_trim_pos_path" "$qf_vis_path"
    local dataset_name="${base##*/}"

    local qza_file="${quality_filter_path%/}/${dataset_name}_QualityFilter.qza"

    # Detect if data is paired-end or single-end
    local qza_type
    qza_type=$(qiime tools peek "$qza_file" | grep "Type:" | sed 's/.*Type:[[:space:]]*//')
    if [[ "$qza_type" == *"PairedEnd"* ]]; then
        # ── PAIRED-END: use dada2 denoise-paired ──
        local need_compute=true
        [[ -n "$start_in" && -n "$end_in" ]] && need_compute=false

        local start_f="" end_f="" start_r="" end_r=""

        if $need_compute; then
            qiime demux summarize \
                --i-data "$qza_file" \
                --o-visualization "${qf_vis_path%/}/${dataset_name}_import_cutadapt_QualityFilter.qzv"
            unzip -q -o "${qf_vis_path%/}/${dataset_name}_import_cutadapt_QualityFilter.qzv" -d "$qf_view_path"

            # Forward trim positions
            find "$qf_view_path" -type f -name 'forward-seven-number-summaries.tsv' -exec cp -f {} "${qf_trim_pos_path%/}/" \;
            local fwd_tsv="${qf_trim_pos_path%/}/forward-seven-number-summaries.tsv"
            if [[ ! -s "$fwd_tsv" ]]; then
                echo "ERROR: Forward seven-number-summaries.tsv not found" >&2
                return 1
            fi
            local fwd_result
            fwd_result="$(python "${SCRIPTS}/py_16s.py" trim_pos_deblur --FilePath "$fwd_tsv")"
            IFS=',' read -r start_f end_f <<< "$fwd_result"

            # Reverse trim positions
            find "$qf_view_path" -type f -name 'reverse-seven-number-summaries.tsv' -exec cp -f {} "${qf_trim_pos_path%/}/" \;
            local rev_tsv="${qf_trim_pos_path%/}/reverse-seven-number-summaries.tsv"
            if [[ -s "$rev_tsv" ]]; then
                local rev_result
                rev_result="$(python "${SCRIPTS}/py_16s.py" trim_pos_deblur --FilePath "$rev_tsv")"
                IFS=',' read -r start_r end_r <<< "$rev_result"
            else
                echo "⚠️ No reverse summary found, using forward positions for reverse"
                start_r="$start_f"
                end_r="$end_f"
            fi

            rm -rf "$qf_view_path"
        fi

        local start="${start_in:-$start_f}"
        local end="${end_in:-$end_f}"
        local start_rev="${start_in:-$start_r}"
        local end_rev="${end_in:-$end_r}"
        echo "$start $end $start_rev $end_rev" > "${qf_trim_pos_path%/}/Trim_position.txt"

        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$qza_file" \
            --p-trunc-len-f "$end" \
            --p-trunc-len-r "$end_rev" \
            --p-trim-left-f "$start" \
            --p-trim-left-r "$start_rev" \
            --o-representative-sequences "${denoising_path%/}/${dataset_name}-rep-seqs-denoising.qza" \
            --o-table "${denoising_path%/}/${dataset_name}-table-denoising.qza" \
            --o-denoising-stats "${denoising_path%/}/${dataset_name}-denoising-stats.qza"
    else
        # ── SINGLE-END: use dada2 denoise-pyro ──
        local need_compute=true
        [[ -n "$start_in" && -n "$end_in" ]] && need_compute=false

        local final_start="" final_end=""

        if $need_compute; then
            qiime demux summarize \
                --i-data "$qza_file" \
                --o-visualization "${qf_vis_path%/}/${dataset_name}_import_cutadapt_QualityFilter.qzv"
            unzip -q -o "${qf_vis_path%/}/${dataset_name}_import_cutadapt_QualityFilter.qzv" -d "$qf_view_path"
            find "$qf_view_path" -type f -name 'forward-seven-number-summaries.tsv' -exec cp -f {} "${qf_trim_pos_path%/}/" \;
            rm -rf "$qf_view_path"
            local tsv_path="${qf_trim_pos_path%/}/forward-seven-number-summaries.tsv"
            if [[ ! -s "$tsv_path" ]]; then
                echo "ERROR: Expected TSV not found: $tsv_path" >&2
                return 1
            fi
            local trim_pos_result
            trim_pos_result="$(python "${SCRIPTS}/py_16s.py" trim_pos_deblur --FilePath "$tsv_path")"
            IFS=',' read -r final_start final_end <<< "$trim_pos_result"
        fi

        local start="${start_in:-$final_start}"
        local end="${end_in:-$final_end}"
        echo "$start $end" > "${qf_trim_pos_path%/}/Trim_position.txt"

        qiime dada2 denoise-pyro \
            --i-demultiplexed-seqs "$qza_file" \
            --p-trunc-len "$end" \
            --p-trim-left "$start" \
            --o-representative-sequences "${denoising_path%/}/${dataset_name}-rep-seqs-denoising.qza" \
            --o-table "${denoising_path%/}/${dataset_name}-table-denoising.qza" \
            --o-denoising-stats "${denoising_path%/}/${dataset_name}-denoising-stats.qza"
    fi
}
################################################################################
#                        LS454 PLATFORM FUNCTIONS                              #
################################################################################
# Functions for 454 pyrosequencing data processing
# 454 pipeline: QC (length filter) → Deduplication → Chimera removal → OTU clustering (97%) → Filter low-freq OTUs
# Quality scores from fasterq-dump are unreliable for 454, so q-score filtering is disabled.

Amplicon_LS454_QualityControlForQZA() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local qza_path="${dataset_path%/}/temp/step_03_qza_import/"
    local quality_filter_path="${dataset_path%/}/temp/step_04_qza_import_QualityFilter/"
    mkdir -p "$quality_filter_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"

    # Use adaptive max_ambiguous if set, otherwise default to 0
    local p_max_ambiguous="${max_ambiguous:-0}"

    # 454 quality scores are unreliable (dummy values from fasterq-dump),
    # so disable q-score filtering (--p-min-quality 0) and apply:
    #   - length filtering (keep reads >= 85% of median length)
    #   - ambiguous base filtering (--p-max-ambiguous from adaptive_tail_trim)
    qiime quality-filter q-score \
        --i-demux "${qza_path%/}/${dataset_name}.qza" \
        --p-min-quality 0 \
        --p-min-length-fraction 0.85 \
        --p-max-ambiguous "$p_max_ambiguous" \
        --o-filtered-sequences "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --o-filter-stats "${quality_filter_path%/}/${dataset_name}_filter-stats.qza" \
        --verbose
}

Amplicon_LS454_Deduplication() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    local quality_filter_path="${dataset_path%/}/temp/step_04_qza_import_QualityFilter/"
    local dedupicate_path="${dataset_path%/}/temp/step_05_dedupicate/"
    mkdir -p "$dedupicate_path"

    # Dereplicate: collapse identical sequences
    # Note: Singleton removal is NOT done here. For 454 data with high error
    # rates, removing singletons before OTU clustering discards reads that
    # would otherwise cluster together. Low-frequency filtering is applied
    # after OTU clustering instead (see Amplicon_LS454_FilterLowFreqOTUs).
    qiime vsearch dereplicate-sequences \
        --i-sequences "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --o-dereplicated-table "${dedupicate_path%/}/${dataset_name}-table.qza" \
        --o-dereplicated-sequences "${dedupicate_path%/}/${dataset_name}-repseq.qza"
}

Amplicon_LS454_ChimerasRemoval() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    local dedupicate_path="${dataset_path%/}/temp/step_05_dedupicate/"
    local chimeras_path="${dataset_path%/}/temp/step_06_ChimerasRemoval/"
    mkdir -p "$chimeras_path"

    # De novo chimera detection (using dereplicated sequences directly,
    # without prior singleton removal)
    qiime vsearch uchime-denovo \
        --i-sequences "${dedupicate_path%/}/${dataset_name}-repseq.qza" \
        --i-table "${dedupicate_path%/}/${dataset_name}-table.qza" \
        --o-chimeras "${chimeras_path%/}/${dataset_name}-chimeras.qza" \
        --o-nonchimeras "${chimeras_path%/}/${dataset_name}-nonchimeras.qza" \
        --o-stats "${chimeras_path%/}/${dataset_name}-uchime-stats.qza"

    # Keep only non-chimeric features in the table
    qiime feature-table filter-features \
        --i-table "${dedupicate_path%/}/${dataset_name}-table.qza" \
        --m-metadata-file "${chimeras_path%/}/${dataset_name}-nonchimeras.qza" \
        --o-filtered-table "${chimeras_path%/}/${dataset_name}-table-nochimera.qza"

    # Filter sequences to non-chimeras
    qiime feature-table filter-seqs \
        --i-data "${dedupicate_path%/}/${dataset_name}-repseq.qza" \
        --m-metadata-file "${chimeras_path%/}/${dataset_name}-nonchimeras.qza" \
        --o-filtered-data "${chimeras_path%/}/${dataset_name}-repseq-nochimera.qza"
}

Amplicon_LS454_ClusterDenovo() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    local chimeras_path="${dataset_path%/}/temp/step_06_ChimerasRemoval/"
    local cluster_path="${dataset_path%/}/temp/step_07_cluster/"
    mkdir -p "$cluster_path"

    # OTU clustering at 97% identity
    qiime vsearch cluster-features-de-novo \
        --i-table "${chimeras_path%/}/${dataset_name}-table-nochimera.qza" \
        --i-sequences "${chimeras_path%/}/${dataset_name}-repseq-nochimera.qza" \
        --p-perc-identity 0.97 \
        --o-clustered-table "${cluster_path%/}/${dataset_name}-table-clustered.qza" \
        --o-clustered-sequences "${cluster_path%/}/${dataset_name}-repseq-clustered.qza"
}

Amplicon_LS454_FilterLowFreqOTUs() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    local cluster_path="${dataset_path%/}/temp/step_07_cluster/"

    # Remove singleton OTUs (total frequency < 2) AFTER clustering.
    # At the OTU level, singletons are truly rare/spurious sequences rather
    # than sequencing-error variants that failed to cluster.
    qiime feature-table filter-features \
        --i-table "${cluster_path%/}/${dataset_name}-table-clustered.qza" \
        --p-min-frequency 2 \
        --o-filtered-table "${dataset_path%/}/${dataset_name}-table-vsearch.qza"

    # Sync representative sequences with filtered OTU table
    qiime feature-table filter-seqs \
        --i-data "${cluster_path%/}/${dataset_name}-repseq-clustered.qza" \
        --i-table "${dataset_path%/}/${dataset_name}-table-vsearch.qza" \
        --o-filtered-data "${dataset_path%/}/${dataset_name}-rep-seqs-vsearch.qza"
}

################################################################################
#                      ION_TORRENT PLATFORM FUNCTIONS                          #
################################################################################
# Functions for Ion Torrent sequencing data processing
# Uses DADA2 denoise-pyro with --p-trim-left 10 to handle Ion Torrent
# signal instability in the first ~10bp.

Amplicon_IonTorrent_QualityControlForQZA() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local qza_path="${dataset_path%/}/temp/step_03_qza_import/"
    local quality_filter_path="${dataset_path%/}/temp/step_04_qza_import_QualityFilter/"
    mkdir -p "$quality_filter_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"

    # Length filter + N removal (no q-score filtering for Ion Torrent)
    # Note: primers are already removed by entropy_primer_detect.py before QIIME2
    # import. The first 10bp trim (Ion Torrent signal instability) is handled
    # downstream by DADA2 denoise-pyro --p-trim-left.
    qiime quality-filter q-score \
        --i-demux "${qza_path%/}/${dataset_name}.qza" \
        --p-min-quality 0 \
        --p-min-length-fraction 0.85 \
        --p-max-ambiguous 0 \
        --o-filtered-sequences "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --o-filter-stats "${quality_filter_path%/}/${dataset_name}_filter-stats.qza" \
        --verbose
}

################################################################################
#                        PACBIO PLATFORM FUNCTIONS                             #
################################################################################
# Functions for PacBio HiFi/CCS long-read sequencing data processing


Amplicon_Pacbio_QualityControlForQZA() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local qza_path="${dataset_path%/}/temp/step_03_qza_import/"
    local quality_filter_path="${dataset_path%/}/temp/step_04_qza_import_QualityFilter/"
    local qf_vis_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/"
    mkdir -p "$quality_filter_path" "$qf_vis_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    qiime quality-filter q-score \
        --i-demux "${qza_path%/}/${dataset_name}.qza" \
        --p-min-length-fraction 0.9 \
        --o-filtered-sequences "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --o-filter-stats "${quality_filter_path%/}/${dataset_name}_filter-stats.qza" \
        --verbose
}
Amplicon_Pacbio_DenosingDada2() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local base="${dataset_path%/}"
    local quality_filter_path="${base}/temp/step_04_qza_import_QualityFilter/"
    local denoising_path="${base}/temp/step_05_denoise/"
    mkdir -p "$denoising_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"

    # Validate that primer sequences are provided
    if [[ -z "$primer_front" ]]; then
        echo "❌ ERROR: No forward primer (primer_front) set for PacBio denoise-ccs."
        echo "   dada2 denoise-ccs requires --p-front to orient CCS reads."
        return 1
    fi
    # Build denoise-ccs command with required parameters
    local cmd=(
        qiime dada2 denoise-ccs
        --i-demultiplexed-seqs "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza"
        --p-front "$primer_front"
        --p-min-len 1000
        --p-max-len 1600
        --o-table "${denoising_path%/}/${dataset_name}-table-denoising.qza"
        --o-representative-sequences "${denoising_path%/}/${dataset_name}-rep-seqs-denoising.qza"
        --o-denoising-stats "${denoising_path%/}/${dataset_name}-denoising-stats.qza"
        --verbose
    )

    # Add reverse primer (--p-adapter) if provided
    if [[ -n "$primer_adapter" ]]; then
        cmd+=(--p-adapter "$primer_adapter")
    fi

    # DADA2 denoise-ccs: handles read orientation, primer removal, denoising,
    # and chimera removal in one step.
    # --p-front: forward primer (27F) — used to orient reads and trim 5' end
    # --p-adapter: reverse primer (1492R) — trims 3' end
    # Reads without primers are discarded; RC reads are re-oriented.
    "${cmd[@]}"
}

Amplicon_Pacbio_ExtractReads() {
    # Extract the V3-V4 region from full-length 16S PacBio rep-seqs so that
    # they are comparable with Illumina V3-V4 amplicon data.
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local base="${dataset_path%/}"
    local denoising_path="${base}/temp/step_05_denoise"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"

    local original="${denoising_path}/${dataset_name}-rep-seqs-denoising.qza"
    local ori_renamed="${denoising_path}/${dataset_name}-ori-rep-seqs-denoising.qza"

    if [[ ! -f "$original" ]]; then
        echo "❌ ERROR: Rep-seqs file not found: $original"
        return 1
    fi

    # Rename the full-length rep-seqs to ori-rep-seqs
    mv "$original" "$ori_renamed"
    echo ">>> Renamed full-length rep-seqs to: $(basename "$ori_renamed")"

    # Extract V3-V4 region using 341F / 785R primers
    echo ">>> Extracting V3-V4 region from full-length 16S rep-seqs..."
    qiime feature-classifier extract-reads \
        --i-sequences "$ori_renamed" \
        --p-f-primer CCTACGGGNGGCWGCAG \
        --p-r-primer GACTACHVGGGTATCTAATCC \
        --p-min-length 300 \
        --p-max-length 500 \
        --p-n-jobs "$cpu" \
        --o-reads "$original"

    echo ">>> V3-V4 extraction complete: $(basename "$original")"
}

################################################################################
#                         ONT PLATFORM FUNCTIONS                               #
################################################################################
# Functions for Oxford Nanopore long-read sequencing data processing
# TODO: These functions are placeholders and need to be implemented

# Amplicon_ONT_PrimerDetectionAndQualityControl() {
#     echo "TODO: Implement ONT primer detection and quality control"
#     echo "  - Adapter trimming for ONT reads"
#     echo "  - Quality filtering for ONT specific error profiles"
#     echo "  - Primer detection and trimming"
#     return 1
# }

# Amplicon_ONT_QualityControlForQZA() {
#     echo "TODO: Implement ONT quality control for QZA"
#     echo "  - Import ONT reads to QIIME2"
#     echo "  - Quality filtering based on ONT characteristics"
#     return 1
# }

# Amplicon_ONT_DenosingMethod() {
#     echo "TODO: Implement ONT denoising method"
#     echo "  - Evaluate denoising approaches for ONT (DADA2? Custom?)"
#     echo "  - Handle ONT-specific error patterns"
#     echo "  - Generate ASV/OTU table"
#     return 1
# }


################################################################################
#                              END OF FILE                                     #
################################################################################

