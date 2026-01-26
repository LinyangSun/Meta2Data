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
Download_CRR_Aspera() {
    local crr=$1
    local target_dir=$2
    local rename_prefix=$3
    local max_retries=3
    local ascp_key=".aspera01.openssh"
    
    echo "[CNCB] Processing $crr"

    # 1. Ensure Aspera key exists
    if [[ ! -f "$ascp_key" ]]; then
        wget -q --user-agent="Mozilla/5.0" "https://ngdc.cncb.ac.cn/gsa/file/downFile?fileName=download/aspera01.openssh" -O "$ascp_key"
        chmod 600 "$ascp_key"
    fi

    # 2. Map CRR to parent CRA ID
    local cra=$(wget -qO- --user-agent="Mozilla/5.0" "https://ngdc.cncb.ac.cn/gsa/search?searchTerm=${crr}" | grep -v "example" | grep -oe "CRA[0-9]\+" | uniq | head -n 1)
    
    if [[ -z "$cra" ]]; then
        echo "Error: Could not map $crr to a CRA project." >&2
        return 1
    fi

    # 3. Fetch project MD5 list and path prefix
    local md5_file="${cra}_md5.txt"
    local path_prefix=""
    if [[ ! -s "$md5_file" ]]; then
        path_prefix=$(wget -qO- --post-data="searchTerm=${cra}&totalDatas=1&downLoadCount=1" --user-agent="Mozilla/5.0" "https://ngdc.cncb.ac.cn/gsa/search/getRunInfoByCra" | grep -oe "gsa[0-9]\+/${cra}\|gsa/${cra}" | uniq | head -n 1)
        wget -q "https://download.cncb.ac.cn/${path_prefix}/md5sum.txt" -O "$md5_file"
    else
        # Get path prefix from existing md5 file or fetch it
        path_prefix=$(wget -qO- --post-data="searchTerm=${cra}&totalDatas=1&downLoadCount=1" --user-agent="Mozilla/5.0" "https://ngdc.cncb.ac.cn/gsa/search/getRunInfoByCra" | grep -oe "gsa[0-9]\+/${cra}\|gsa/${cra}" | uniq | head -n 1)
    fi

    echo "  Using path prefix: $path_prefix"

    # 4. Identify files belonging to the CRR accession
    local -a files
    while IFS= read -r filepath; do
        files+=("$(basename "$filepath")")
    done < <(grep "$crr" "$md5_file" | awk '{print $2}')

    for filename in "${files[@]}"; do
        local expected_md5=$(grep "$filename" "$md5_file" | awk '{print tolower($1)}')
        local success=false

        # Try Aspera first
        echo "  Trying Aspera for $filename..."
        for ((i=1; i<=max_retries; i++)); do
            # Aspera download (Quiet mode)
            if command -v ascp >/dev/null 2>&1 && ascp -P 33001 -i "$ascp_key" -QT -l 200M -k1 -d "aspera01@download.cncb.ac.cn:gsa/${cra}/${crr}/${filename}" . > /dev/null 2>&1; then
                # MD5 Verification
                if [[ -f "$filename" ]]; then
                    local current_md5=$(md5sum "$filename" | awk '{print $1}')
                    if [[ "$current_md5" == "$expected_md5" ]]; then
                        # Standardize filename: _r1/_R1 → _1, _r2/_R2 → _2
                        local new_filename=$(echo "$filename" | sed -E 's/_[rR]1([._])/_1\1/; s/_[rR]2([._])/_2\1/')
                        # Strip CRR ID from filename to avoid duplication
                        local base_filename=$(echo "$new_filename" | sed "s/${crr}[_-]//")
                        mv "$filename" "${target_dir}/${rename_prefix}-${base_filename}"
                        success=true
                        echo "  ✓ Downloaded via Aspera: $filename → ${rename_prefix}-${base_filename}"
                        break
                    else
                        rm -f "$filename"
                    fi
                fi
            fi
            [[ $i -lt $max_retries ]] && sleep 2
        done

        # Fallback to wget if Aspera failed
        if [[ "$success" = false ]]; then
            echo "  Aspera failed, trying wget for $filename..."
            for ((i=1; i<=max_retries; i++)); do
                local wget_url="https://download.cncb.ac.cn/${path_prefix}/${crr}/${filename}"
                echo "    Downloading from: $wget_url"
                if wget -q --user-agent="Mozilla/5.0" "$wget_url" -O "$filename"; then
                    # MD5 Verification
                    if [[ -f "$filename" ]]; then
                        local current_md5=$(md5sum "$filename" | awk '{print $1}')
                        if [[ "$current_md5" == "$expected_md5" ]]; then
                            # Standardize filename: _r1/_R1 → _1, _r2/_R2 → _2
                            local new_filename=$(echo "$filename" | sed -E 's/_[rR]1([._])/_1\1/; s/_[rR]2([._])/_2\1/')
                            # Strip CRR ID from filename to avoid duplication
                            local base_filename=$(echo "$new_filename" | sed "s/${crr}[_-]//")
                            mv "$filename" "${target_dir}/${rename_prefix}-${base_filename}"
                            success=true
                            echo "  ✓ Downloaded via wget: $filename → ${rename_prefix}-${base_filename}"
                            break
                        else
                            echo "  MD5 mismatch, retrying..."
                            rm -f "$filename"
                        fi
                    fi
                fi
                [[ $i -lt $max_retries ]] && sleep 2
            done
        fi

        if [[ "$success" = false ]]; then
            echo "Error: Failed to download $filename after trying both Aspera and wget." >&2
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

    # Check dependencies
    command -v prefetch >/dev/null 2>&1 || { echo "Error: 'prefetch' not found." >&2; return 1; }
    command -v fasterq-dump >/dev/null 2>&1 || { echo "Error: 'fasterq-dump' not found." >&2; return 1; }
    command -v ascp >/dev/null 2>&1 || { echo "Error: 'ascp' not found." >&2; return 1; }
    
    local base_dir="${dir_path%/}/"
    local fastq_path="${base_dir}ori_fastq/"
    local temp_dl_path="${base_dir}temp1/"
    mkdir -p "$fastq_path" "$temp_dl_path"

    # Phase 1: Data Acquisition
    while IFS=$'\t' read -r srr rename _; do
        [[ -z "$srr" || -z "$rename" ]] && continue

        # Route by ID type (CRR vs SRR/ERR/DRR)
        if [[ "$srr" =~ ^CRR ]]; then
            # CRR files go directly to ori_fastq/ with rename prefix
            Download_CRR_Aspera "$srr" "$fastq_path" "$rename"

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
                            # Strip SRR/ERR/DRR ID from filename to avoid duplication
                            local base_filename=$(echo "$fq_basename" | sed "s/${srr}[_-]//; s/^_//")
                            # Standardize to _1/_2 format
                            base_filename=$(echo "$base_filename" | sed -E 's/_[rR]?1([._])/_1\1/; s/_[rR]?2([._])/_2\1/')
                            mv "$fq_file" "${fastq_path}${rename}-${base_filename}"
                        done
                        rm -rf "$temp_outdir"
                    else
                        echo "Error: fasterq-dump failed for $basename_file" >&2
                    fi
                    rm -f "$file"

                elif [[ "$ext" == "fastq" || "$basename_file" =~ \.fastq\.gz$ ]]; then
                    # Already FASTQ, strip SRR/ERR/DRR ID and rename
                    local base_filename=$(echo "$basename_file" | sed "s/${srr}[_-]//; s/^_//")
                    # Standardize to _1/_2 format
                    base_filename=$(echo "$base_filename" | sed -E 's/_[rR]?1([._])/_1\1/; s/_[rR]?2([._])/_2\1/')
                    mv "$file" "${fastq_path}${rename}-${base_filename}"
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

Amplicon_Common_MakeManifestFileForQiime2() {
    cd $dataset_path
    local temp_file_path=$dataset_path"temp/temp_file/"
    local fastp_path=$dataset_path"temp/step_02_fastp/"
    mkdir -p "$temp_file_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    find $fastp_path -type f > $temp_file_path$dataset_name"-file.txt"
    if [ "$sequence_type" = "single" ]; then
        py_16s.py mk_manifest_SE --FilePath $temp_file_path$dataset_name"-file.txt"
    else
        py_16s.py mk_manifest_PE --FilePath $temp_file_path$dataset_name"-file.txt"
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
    local forward_manifest="${temp_file_path}${dataset_name}_manifest_forward.tsv"
    local merge_log="${temp_path}${dataset_name}_vsearch_merge.log"

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
    echo "🧬 Importing paired-end reads..."
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$paired_manifest" \
        --output-path "${qza_path}${dataset_name}.qza" \
        --input-format PairedEndFastqManifestPhred33V2

    echo "🔄 Merging paired-end reads..."
    # Try merging, but don't exit on failure; capture logs.
    if ! qiime vsearch merge-pairs \
        --i-demultiplexed-seqs "${qza_path}${dataset_name}.qza" \
        --o-unmerged-sequences "${qza_path}${dataset_name}_unmerged.qza" \
        --o-merged-sequences "${qza_path}${dataset_name}_join.qza" \
        &> "$merge_log"
    then
        echo "⚠️ vsearch merge failed. See log: $merge_log"
        echo "↩️ Falling back to forward reads only..."

        # Build forward-only manifest from paired manifest
        awk -F'\t' '
            BEGIN { OFS="\t" }
            NR==1 { 
                print "sample-id","absolute-filepath"
                next 
            }
            NR>1 {
                print $1, $2
            }
        ' "$paired_manifest" > "$forward_manifest"

        # Replace any existing paired artifact with forward-only import
        rm -f "${qza_path}${dataset_name}.qza"
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path "$forward_manifest" \
            --output-path "${qza_path}${dataset_name}.qza" \
            --input-format SingleEndFastqManifestPhred33V2

        echo "✅ Forward-only import completed after merge failure."
        return
    fi

    # --- CHECK MERGE QUALITY (only if merge succeeded) ---
    local unmerged_qza="${qza_path}${dataset_name}_unmerged.qza"
    local merged_qza="${qza_path}${dataset_name}_join.qza"
    local unmerged_size merged_size

    unmerged_size=$(stat -c%s "$unmerged_qza")
    merged_size=$(stat -c%s "$merged_qza")

    echo "📏 unmerged.qza = $unmerged_size bytes"
    echo "📏 merged.qza   = $merged_size bytes"

    if (( $(echo "$unmerged_size * 1.5 > $merged_size" | bc -l) )); then
        echo "⚠️ Less than ~70% merged — switching to forward reads only..."

        awk -F'\t' '
            BEGIN { OFS="\t" }
            NR==1 { 
                print "sample-id","absolute-filepath"
                next 
            }
            NR>1 {
                print $1, $2
            }
        ' "$paired_manifest" > "$forward_manifest"

        rm -f "${qza_path}${dataset_name}.qza"
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path "$forward_manifest" \
            --output-path "${qza_path}${dataset_name}.qza" \
            --input-format SingleEndFastqManifestPhred33V2

        echo "✅ Forward-only import completed."
    else
        echo "✅ Merge quality acceptable. Using merged reads."
        rm -f "${qza_path}${dataset_name}.qza"
        mv "$merged_qza" "${qza_path}${dataset_name}.qza"
    fi
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
    while getopts ":s:e:" opt; do
        case "$opt" in
            s) start_in="$OPTARG" ;;
            e) end_in="$OPTARG" ;;
            \?) echo "Invalid option: -$OPTARG" >&2; return 2 ;;
            :)  echo "Option -$OPTARG requires an argument." >&2; return 2 ;;
        esac
    done
    shift $((OPTIND-1))
    # Safety
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
    # Decide if we need to compute positions
    local need_compute=true
    [[ -n "$start_in" && -n "$end_in" ]] && need_compute=false
    echo  "start_in: $start_in ; end_in $end_in"
    local start_calc="" end_calc=""
    echo "need_compute is : $need_compute"
    if $need_compute; then
        echo "Computing trim positions from QIIME 2 visualization..."
        qiime demux summarize \
            --i-data "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
            --o-visualization "${qf_vis_path%/}/${dataset_name}_import_cutadapt_QualityFilter.qzv"
        unzip -q -o "${qf_vis_path%/}/${dataset_name}_import_cutadapt_QualityFilter.qzv" -d "$qf_view_path"
        # Fixed wrong var name from original; keep forward file
        find "$qf_view_path" -type f -name 'forward-seven-number-summaries.tsv' -exec cp -f {} "${qf_trim_pos_path%/}/" \;
        rm -rf "$qf_view_path"
        local tsv_path="${qf_trim_pos_path%/}/forward-seven-number-summaries.tsv"
        if [[ ! -s "$tsv_path" ]]; then
            echo "ERROR: Expected TSV not found: $tsv_path" >&2
            return 1
        fi
        local trim_pos_result
        trim_pos_result="$(py_16s.py trim_pos_deblur --FilePath "$tsv_path")"
        echo "$trim_pos_result"
        IFS=',' read -r final_start final_end <<< "$trim_pos_result"
        echo "Computed: start=$final_start, end=$final_end"
    fi
    # Final positions: use flags if set, otherwise computed
    local start="${start_in:-$final_start}"
    local end="${end_in:-$final_end}"
    echo "Using trim positions: start=$start, end=$end"
    echo "$start $end" > "${qf_trim_pos_path%/}/Trim_position.txt"
    # Deblur
    qiime dada2 denoise-pyro \
        --i-demultiplexed-seqs "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --p-trunc-len "$end" \
        --p-trim-left "$start" \
        --o-representative-sequences "${denoising_path%/}/${dataset_name}-rep-seqs-denoising.qza" \
        --o-table "${denoising_path%/}/${dataset_name}-table-denoising.qza" \
        --o-denoising-stats "${denoising_path%/}/${dataset_name}-denoising-stats.qza"
    echo "denoising complete."
}
################################################################################
#                        PACBIO PLATFORM FUNCTIONS                             #
################################################################################
# Functions for PacBio HiFi/CCS long-read sequencing data processing


Amplicon_Pacbio_PrimerDetectionAndQualityControl() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local fastq_path="${dataset_path%/}/ori_fastq/"
    local fastp_path="${dataset_path%/}/temp/step_02_fastp/"
    local temp_file_path="${dataset_path%/}/temp/temp_file/"
    mkdir -p "$temp_file_path" "$fastp_path"
    echo "The data is Pacbio, no need to remove primer before analysis, just copy original file to fastp folder"
    cp "${fastq_path%/}"/*.fastq "${fastp_path%/}/"

}
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
    # Paths
    local base="${dataset_path%/}"
    local quality_filter_path="${base%/}/temp/step_04_qza_import_QualityFilter/"
    local denoising_path="${base%/}/temp/step_05_denoise/"
    local qf_vis_path="${base%/}/temp/temp_file/QualityFilter_vis/"
    local qf_view_path="${base%/}/temp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${base%/}/temp/temp_file/QualityFilter_vis/qf_trim_pos/"
    mkdir -p "$denoising_path" "$qf_view_path" "$qf_trim_pos_path" "$qf_vis_path"
    # Deblur
    qiime dada2 denoise-ccs \
        --i-demultiplexed-seqs "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --p-max-len 1600 \
        --o-table "${denoising_path%/}/${dataset_name}-table-denosing.qza" \
        --o-representative-sequences "${denoising_path%/}/${dataset_name}-rep-seqs-denosing.qza" \
        --o-denoising-stats "${denoising_path%/}/${dataset_name}-denosing-stats.qza"
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

