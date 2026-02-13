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
Download_CRR() {
    local crr=$1
    local target_dir=$2
    local rename_prefix=$3
    local max_retries=3

    echo "[CNCB] Processing $crr"

    # 1. Map CRR to parent CRA ID
    local cra=$(wget -qO- --user-agent="Mozilla/5.0" "https://ngdc.cncb.ac.cn/gsa/search?searchTerm=${crr}" | grep -v "example" | grep -oe "CRA[0-9]\+" | uniq | head -n 1)

    if [[ -z "$cra" ]]; then
        echo "Error: Could not map $crr to a CRA project." >&2
        return 1
    fi

    # 2. Fetch project MD5 list and path prefix
    local md5_file="${cra}_md5.txt"
    local path_prefix=""
    path_prefix=$(wget -qO- --post-data="searchTerm=${cra}&totalDatas=1&downLoadCount=1" --user-agent="Mozilla/5.0" "https://ngdc.cncb.ac.cn/gsa/search/getRunInfoByCra" | grep -oe "gsa[0-9]\+/${cra}\|gsa/${cra}" | uniq | head -n 1)
    if [[ ! -s "$md5_file" ]]; then
        wget -q "https://download.cncb.ac.cn/${path_prefix}/md5sum.txt" -O "$md5_file"
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
            if wget -q --user-agent="Mozilla/5.0" "$wget_url" -O "$filename"; then
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
            fi
            [[ $i -lt $max_retries ]] && sleep 2
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
    local base_dir="${dir_path%/}/"

    while IFS=$'\t' read -r srr _; do
        [[ -z "$srr" ]] && continue
        if [[ "$srr" =~ ^CRR ]]; then
            has_cncb_accessions=true
        elif [[ "$srr" =~ ^[EDS]RR ]]; then
            has_ncbi_accessions=true
        fi
    done < "${base_dir}${acc_file}"

    # Check dependencies based on what we'll download
    if [[ "$has_ncbi_accessions" == true ]]; then
        command -v prefetch >/dev/null 2>&1 || { echo "Error: 'prefetch' not found (required for NCBI SRR/ERR/DRR downloads)." >&2; return 1; }
        command -v fasterq-dump >/dev/null 2>&1 || { echo "Error: 'fasterq-dump' not found (required for NCBI SRR/ERR/DRR downloads)." >&2; return 1; }
    fi

    if [[ "$has_cncb_accessions" == true ]]; then
        command -v wget >/dev/null 2>&1 || { echo "Error: 'wget' not found (required for CNCB CRR downloads)." >&2; return 1; }
    fi

    # base_dir already declared above
    local fastq_path="${base_dir}ori_fastq/"
    local temp_dl_path="${base_dir}temp1/"
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
            local attempt=0
            # SRA Toolkit prefetch (Quiet mode)
            until prefetch -q -O "$temp_dl_path" "$srr" > /dev/null 2>&1; do
                (( attempt++ >= 10 )) && { echo "Error: Download failed for $srr" >&2; break; }
                sleep 5
            done

            # Extract downloaded SRA files to temp
            local sra_files=()
            while IFS= read -r -d '' file; do
                sra_files+=("$file")
            done < <(find "$temp_dl_path/$srr" -type f \( -name "*.sra" -o -name "*.fastq*" \) -print0 2>/dev/null)

            # Process each file
            for file in "${sra_files[@]}"; do
                local basename_file=$(basename "$file")
                local ext="${basename_file##*.}"

                if [[ "$ext" == "sra" ]]; then
                    # Convert SRA to FASTQ with rename prefix
                    local temp_outdir="${temp_dl_path}/fastq_${srr}"
                    mkdir -p "$temp_outdir"

                    if fasterq-dump --split-3 -q "$file" --outdir "$temp_outdir" > /dev/null 2>&1; then
                        # Rename and move converted FASTQ files
                        find "$temp_outdir" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) -print0 |
                        while IFS= read -r -d '' fq_file; do
                            local fq_basename=$(basename "$fq_file")
                            # Standardize to _1/_2 format first
                            local normalized=$(echo "$fq_basename" | sed -E 's/_[rR]?1([._])/_1\1/; s/_[rR]?2([._])/_2\1/')
                            # Strip SRR/ERR/DRR ID from filename (keep separator)
                            local base_filename=$(echo "$normalized" | sed "s/${srr}//")
                            mv "$fq_file" "${fastq_path}${rename}${base_filename}"
                        done
                        rm -rf "$temp_outdir"
                    else
                        echo "Error: fasterq-dump failed for $basename_file" >&2
                    fi
                    rm -f "$file"

                elif [[ "$ext" == "fastq" || "$basename_file" =~ \.fastq\.gz$ ]]; then
                    # Already FASTQ, standardize and strip SRR/ERR/DRR ID
                    # Standardize to _1/_2 format first
                    local normalized=$(echo "$basename_file" | sed -E 's/_[rR]?1([._])/_1\1/; s/_[rR]?2([._])/_2\1/')
                    # Strip SRR/ERR/DRR ID from filename (keep separator)
                    local base_filename=$(echo "$normalized" | sed "s/${srr}//")
                    mv "$file" "${fastq_path}${rename}${base_filename}"
                fi
            done

            rm -rf "${temp_dl_path:?}/$srr"
        else
            echo "Warning: Unknown Accession format: $srr" >&2
        fi
    done < "${base_dir}${acc_file}"
    rm -rf "$temp_dl_path"

    # [Note: Subsequent vsearch merging and reporting logic should follow here]

    echo "Process finished. Output directory: $fastq_path"
}


Common_SRADownloadViaSFF() {
    # Download SRA data via SFF format for LS454/ION_TORRENT platforms.
    # SFF preserves original flowgram signals and quality scores, which are
    # critical for dada2 denoise-pyro.
    #
    # Flow: prefetch → sff-dump → sff2fastq → FASTQ (with quality scores)
    # Fallback: prefetch → fasterq-dump → FASTQ (quality may be missing)
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
        echo "Usage: Common_SRADownloadViaSFF -d <dir> -a <accession_tsv>" >&2
        return 1
    fi

    local base_dir="${dir_path%/}/"
    local fastq_path="${base_dir}ori_fastq/"
    local temp_dl_path="${base_dir}temp1/"
    mkdir -p "$fastq_path" "$temp_dl_path"

    # Pre-scan for dependency check
    local has_ncbi=false has_cncb=false
    while IFS=$'\t' read -r srr _; do
        [[ -z "$srr" ]] && continue
        if [[ "$srr" =~ ^CRR ]]; then has_cncb=true
        elif [[ "$srr" =~ ^[EDS]RR ]]; then has_ncbi=true; fi
    done < "${base_dir}${acc_file}"

    if [[ "$has_ncbi" == true ]]; then
        command -v prefetch >/dev/null 2>&1 || { echo "Error: 'prefetch' not found" >&2; return 1; }
        if ! command -v sff-dump >/dev/null 2>&1; then
            echo "⚠️ 'sff-dump' not found; will fall back to fasterq-dump"
        fi
        if ! command -v sff2fastq >/dev/null 2>&1; then
            echo "⚠️ 'sff2fastq' not found; will fall back to fasterq-dump"
        fi
    fi
    if [[ "$has_cncb" == true ]]; then
        command -v wget >/dev/null 2>&1 || { echo "Error: 'wget' not found" >&2; return 1; }
    fi

    while IFS=$'\t' read -r srr rename _; do
        [[ -z "$srr" || -z "$rename" ]] && continue

        if [[ "$srr" =~ ^CRR ]]; then
            # CNCB: download via wget (already FASTQ)
            Download_CRR "$srr" "$fastq_path" "$rename"

        elif [[ "$srr" =~ ^[EDS]RR ]]; then
            echo "[NCBI/SFF] Processing $srr"

            # Step 1: prefetch .sra file
            local attempt=0
            until prefetch -q -O "$temp_dl_path" "$srr" > /dev/null 2>&1; do
                (( attempt++ >= 10 )) && { echo "Error: prefetch failed for $srr" >&2; break; }
                sleep 5
            done

            local sra_file
            sra_file=$(find "$temp_dl_path/$srr" -name "*.sra" -type f 2>/dev/null | head -1)

            # Step 2: Try sff-dump → sff2fastq (preserves quality)
            local used_sff=false
            if [[ -n "$sra_file" ]] && command -v sff-dump >/dev/null 2>&1 && command -v sff2fastq >/dev/null 2>&1; then
                local sff_outdir="${temp_dl_path}/sff_${srr}"
                mkdir -p "$sff_outdir"

                if sff-dump --outdir "$sff_outdir" "$sra_file" 2>/dev/null; then
                    for sff_file in "$sff_outdir"/*.sff; do
                        [[ -f "$sff_file" ]] || continue
                        sff2fastq "$sff_file" -o "${fastq_path}${rename}.fastq"
                        echo "  ✓ SFF→FASTQ: ${srr} → ${rename}.fastq (quality scores preserved)"
                        used_sff=true
                    done
                else
                    echo "  ⚠️ sff-dump failed (data may not be SFF format)"
                fi
                rm -rf "$sff_outdir"
            fi

            # Step 3: Fallback to fasterq-dump
            if [[ "$used_sff" == false ]]; then
                echo "  ↩️ Falling back to fasterq-dump..."
                if [[ -n "$sra_file" ]]; then
                    local temp_fq="${temp_dl_path}/fastq_${srr}"
                    mkdir -p "$temp_fq"

                    if fasterq-dump --split-3 -q "$sra_file" --outdir "$temp_fq" 2>/dev/null; then
                        find "$temp_fq" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) -print0 |
                        while IFS= read -r -d '' fq_file; do
                            local fq_basename=$(basename "$fq_file")
                            local normalized=$(echo "$fq_basename" | sed -E 's/_[rR]?1([._])/_1\1/; s/_[rR]?2([._])/_2\1/')
                            local base_filename=$(echo "$normalized" | sed "s/${srr}//")
                            mv "$fq_file" "${fastq_path}${rename}${base_filename}"
                            echo "  ✓ fasterq-dump: ${srr} → ${rename}${base_filename} (quality may be missing)"
                        done
                    else
                        echo "Error: Both sff-dump and fasterq-dump failed for $srr" >&2
                    fi
                    rm -rf "$temp_fq"
                fi
            fi

            rm -rf "${temp_dl_path:?}/$srr"
        else
            echo "Warning: Unknown accession format: $srr" >&2
        fi
    done < "${base_dir}${acc_file}"

    rm -rf "$temp_dl_path"
    echo "Process finished (SFF mode). Output directory: $fastq_path"
}

Amplicon_Common_MakeManifestFileForQiime2() {
    cd $dataset_path
    local temp_file_path=$dataset_path"temp/temp_file/"
    mkdir -p "$temp_file_path"
    local _dp="${dataset_path%/}"
    dataset_name="${_dp##*/}"
    find "$fastq_path" -type f -name "*.fastq*" > $temp_file_path$dataset_name"-file.txt"
    if [ "$sequence_type" = "single" ]; then
        python "${SCRIPTS}/py_16s.py" mk_manifest_SE --FilePath $temp_file_path$dataset_name"-file.txt"
    else
        python "${SCRIPTS}/py_16s.py" mk_manifest_PE --FilePath $temp_file_path$dataset_name"-file.txt"
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

    echo "🔹 Processing dataset: $dataset_name"

    # --- SINGLE-END ---
    if [ "${sequence_type:-paired}" = "single" ]; then
        echo "🧬 Importing single-end reads..."
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path "$paired_manifest" \
            --output-path "${qza_path}${dataset_name}.qza" \
            --input-format SingleEndFastqManifestPhred33V2
        echo "✅ Single-end import completed."
        return
    fi

    # --- PAIRED-END ---
    # Import as paired-end; DADA2 denoise-paired handles merging internally
    echo "🧬 Importing paired-end reads..."
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$paired_manifest" \
        --output-path "${qza_path}${dataset_name}.qza" \
        --input-format PairedEndFastqManifestPhred33V2
    echo "✅ Paired-end import completed."
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
        echo "The analysis ran successfully! Now only keeping final files to save space."
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

        rm -rf "${dataset_path}ori_fastq/" 2>/dev/null || true
        rm -rf "${dataset_path}working_fastq/" 2>/dev/null || true
        
        return 0
        
    # Check for vsearch output (454 pipeline)
    elif [ -f "${dataset_path%/}/${dataset_name}-table-vsearch.qza" ]; then
        echo "The analysis ran successfully! Now only keeping final files to save space."
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
    local base_dir="${1%/}/"
    local acc_file="$2"
    local fastq_path="${base_dir}ori_fastq/"
    local dataset_name="${base_dir%/}"
    dataset_name="${dataset_name##*/}"
    local raw_counts_file="${base_dir}${dataset_name}_raw_read_counts.tsv"
    : > "$raw_counts_file"

    echo ">>> Counting raw reads per sample..."
    while IFS=$'\t' read -r srr rename _; do
        [[ -z "$srr" || -z "$rename" ]] && continue
        local total_lines=0
        for fq in "${fastq_path}${rename}"*.fastq*; do
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
    done < "${base_dir}${acc_file}"

    echo "✓ Raw read counts saved to $raw_counts_file"
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
    echo "📊 Detected QZA type: $qza_type"

    if [[ "$qza_type" == *"PairedEnd"* ]]; then
        # ── PAIRED-END: use dada2 denoise-paired ──
        echo "📊 Using dada2 denoise-paired for paired-end data..."
        local need_compute=true
        [[ -n "$start_in" && -n "$end_in" ]] && need_compute=false

        local start_f="" end_f="" start_r="" end_r=""

        if $need_compute; then
            echo "Computing trim positions from QIIME 2 visualization..."
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
            echo "Forward: $fwd_result"
            IFS=',' read -r start_f end_f <<< "$fwd_result"

            # Reverse trim positions
            find "$qf_view_path" -type f -name 'reverse-seven-number-summaries.tsv' -exec cp -f {} "${qf_trim_pos_path%/}/" \;
            local rev_tsv="${qf_trim_pos_path%/}/reverse-seven-number-summaries.tsv"
            if [[ -s "$rev_tsv" ]]; then
                local rev_result
                rev_result="$(python "${SCRIPTS}/py_16s.py" trim_pos_deblur --FilePath "$rev_tsv")"
                echo "Reverse: $rev_result"
                IFS=',' read -r start_r end_r <<< "$rev_result"
            else
                echo "⚠️ No reverse summary found, using forward positions for reverse"
                start_r="$start_f"
                end_r="$end_f"
            fi

            rm -rf "$qf_view_path"
            echo "Computed: forward start=$start_f, end=$end_f; reverse start=$start_r, end=$end_r"
        fi

        local start="${start_in:-$start_f}"
        local end="${end_in:-$end_f}"
        local start_rev="${start_in:-$start_r}"
        local end_rev="${end_in:-$end_r}"
        echo "Using trim positions: forward start=$start, end=$end; reverse start=$start_rev, end=$end_rev"
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
        echo "📊 Using dada2 denoise-pyro for single-end data..."
        local need_compute=true
        [[ -n "$start_in" && -n "$end_in" ]] && need_compute=false
        echo "start_in: $start_in ; end_in: $end_in"
        echo "need_compute is: $need_compute"

        local final_start="" final_end=""

        if $need_compute; then
            echo "Computing trim positions from QIIME 2 visualization..."
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
            echo "$trim_pos_result"
            IFS=',' read -r final_start final_end <<< "$trim_pos_result"
            echo "Computed: start=$final_start, end=$final_end"
        fi

        local start="${start_in:-$final_start}"
        local end="${end_in:-$final_end}"
        echo "Using trim positions: start=$start, end=$end"
        echo "$start $end" > "${qf_trim_pos_path%/}/Trim_position.txt"

        qiime dada2 denoise-pyro \
            --i-demultiplexed-seqs "$qza_file" \
            --p-trunc-len "$end" \
            --p-trim-left "$start" \
            --o-representative-sequences "${denoising_path%/}/${dataset_name}-rep-seqs-denoising.qza" \
            --o-table "${denoising_path%/}/${dataset_name}-table-denoising.qza" \
            --o-denoising-stats "${denoising_path%/}/${dataset_name}-denoising-stats.qza"
    fi
    echo "denoising complete."
}
################################################################################
#                        LS454 PLATFORM FUNCTIONS                              #
################################################################################
# Functions for 454 pyrosequencing data processing
# 454 pipeline: QC (length filter) → Deduplication → Chimera removal → OTU clustering (97%)
# Quality scores from fasterq-dump are unreliable for 454, so q-score filtering is disabled.

Amplicon_LS454_QualityControlForQZA() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local qza_path="${dataset_path%/}/temp/step_03_qza_import/"
    local quality_filter_path="${dataset_path%/}/temp/step_04_qza_import_QualityFilter/"
    mkdir -p "$quality_filter_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"

    # 454 quality scores are unreliable (dummy values from fasterq-dump),
    # so disable q-score filtering (--p-min-quality 0) and only apply:
    #   - length filtering (keep reads >= 85% of median length)
    #   - ambiguous base removal (remove reads containing N)
    qiime quality-filter q-score \
        --i-demux "${qza_path%/}/${dataset_name}.qza" \
        --p-min-quality 0 \
        --p-min-length-fraction 0.85 \
        --p-max-ambiguous 0 \
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
    qiime vsearch dereplicate-sequences \
        --i-sequences "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --o-dereplicated-table "${dedupicate_path%/}/${dataset_name}-table.qza" \
        --o-dereplicated-sequences "${dedupicate_path%/}/${dataset_name}-repseq.qza"

    # Remove singletons (features appearing only once)
    qiime feature-table filter-features \
        --i-table "${dedupicate_path%/}/${dataset_name}-table.qza" \
        --p-min-frequency 2 \
        --o-filtered-table "${dedupicate_path%/}/${dataset_name}-table-nosingle.qza"

    # Sync sequences with filtered table
    qiime feature-table filter-seqs \
        --i-data "${dedupicate_path%/}/${dataset_name}-repseq.qza" \
        --i-table "${dedupicate_path%/}/${dataset_name}-table-nosingle.qza" \
        --o-filtered-data "${dedupicate_path%/}/${dataset_name}-repseq-nosingle.qza"
}

Amplicon_LS454_ChimerasRemoval() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    local dedupicate_path="${dataset_path%/}/temp/step_05_dedupicate/"
    local chimeras_path="${dataset_path%/}/temp/step_06_ChimerasRemoval/"
    mkdir -p "$chimeras_path"

    # De novo chimera detection
    qiime vsearch uchime-denovo \
        --i-sequences "${dedupicate_path%/}/${dataset_name}-repseq-nosingle.qza" \
        --i-table "${dedupicate_path%/}/${dataset_name}-table-nosingle.qza" \
        --o-chimeras "${chimeras_path%/}/${dataset_name}-chimeras.qza" \
        --o-nonchimeras "${chimeras_path%/}/${dataset_name}-nonchimeras.qza" \
        --o-stats "${chimeras_path%/}/${dataset_name}-uchime-stats.qza"

    # Keep only non-chimeric features in the table
    qiime feature-table filter-features \
        --i-table "${dedupicate_path%/}/${dataset_name}-table-nosingle.qza" \
        --m-metadata-file "${chimeras_path%/}/${dataset_name}-nonchimeras.qza" \
        --o-filtered-table "${chimeras_path%/}/${dataset_name}-table-nochimera.qza"

    # Filter sequences to non-chimeras
    qiime feature-table filter-seqs \
        --i-data "${dedupicate_path%/}/${dataset_name}-repseq-nosingle.qza" \
        --m-metadata-file "${chimeras_path%/}/${dataset_name}-nonchimeras.qza" \
        --o-filtered-data "${chimeras_path%/}/${dataset_name}-repseq-nochimera.qza"
}

Amplicon_LS454_ClusterDenovo() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    local chimeras_path="${dataset_path%/}/temp/step_06_ChimerasRemoval/"

    # OTU clustering at 97% identity
    qiime vsearch cluster-features-de-novo \
        --i-table "${chimeras_path%/}/${dataset_name}-table-nochimera.qza" \
        --i-sequences "${chimeras_path%/}/${dataset_name}-repseq-nochimera.qza" \
        --p-perc-identity 0.97 \
        --o-clustered-table "${dataset_path%/}/${dataset_name}-table-vsearch.qza" \
        --o-clustered-sequences "${dataset_path%/}/${dataset_name}-rep-seqs-vsearch.qza"
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

    # Validate that a primer was detected
    if [[ -z "$detected_primer" ]]; then
        echo "❌ ERROR: No primer detected for PacBio denoise-ccs."
        echo "   dada2 denoise-ccs requires --p-front to orient CCS reads."
        return 1
    fi
    echo "  Using detected primer for --p-front: $detected_primer"

    # DADA2 denoise-ccs: handles read orientation, primer removal, denoising,
    # and chimera removal in one step
    qiime dada2 denoise-ccs \
        --i-demultiplexed-seqs "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --p-front "$detected_primer" \
        --p-min-len 1000 \
        --p-max-len 1600 \
        --o-table "${denoising_path%/}/${dataset_name}-table-denoising.qza" \
        --o-representative-sequences "${denoising_path%/}/${dataset_name}-rep-seqs-denoising.qza" \
        --o-denoising-stats "${denoising_path%/}/${dataset_name}-denoising-stats.qza" \
        --verbose
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

