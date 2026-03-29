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

Common_SanitizeFastq() {
    # Remove FASTQ records with sequence shorter than 50bp.
    # For PE: drops both reads if either is too short.
    # Uses $fastq_path and $sequence_type from the calling scope.
    echo ">>> Sanitizing FASTQ files (removing reads < 50bp)..."
    python3 "${SCRIPTS}/py_16s.py" sanitize_fastq \
        --input_dir "$fastq_path" \
        --min_length 50 \
        --sequence_type "$sequence_type"
}

Download_From_ENA() {
    # Download FASTQ file(s) directly from ENA (European Nucleotide Archive).
    # ENA provides pre-converted FASTQ files with original quality scores
    # (critical for DADA2 error learning).
    #
    # Args:
    #   $1 = SRR/ERR/DRR accession
    #   $2 = target directory for FASTQ files
    #   $3 = rename prefix for output files
    #
    # Output (silent on success): only prints errors/warnings to stderr.
    # Sets _ena_layout and _ena_protocol on first successful download for subsequent calls.
    # Returns: 0 on success, 1 on failure
    local srr="$1"
    local target_dir="$2"
    local rename_prefix="$3"
    local max_retries=5

    # ENA URL structure: first 6 chars of accession, then full accession
    local acc_prefix="${srr:0:6}"
    local acc_len=${#srr}

    # Determine protocol(s) to try: use cached protocol if detected, otherwise try both
    local -a protocols
    if [[ -n "$_ena_protocol" ]]; then
        protocols=("$_ena_protocol")
    else
        protocols=("https" "ftp")
    fi

    local ena_path
    if (( acc_len <= 9 )); then
        ena_path="ftp.sra.ebi.ac.uk/vol1/fastq/${acc_prefix}/${srr}"
    elif (( acc_len == 10 )); then
        ena_path="ftp.sra.ebi.ac.uk/vol1/fastq/${acc_prefix}/00${srr: -1}/${srr}"
    elif (( acc_len == 11 )); then
        ena_path="ftp.sra.ebi.ac.uk/vol1/fastq/${acc_prefix}/0${srr: -2}/${srr}"
    else
        ena_path="ftp.sra.ebi.ac.uk/vol1/fastq/${acc_prefix}/${srr: -3}/${srr}"
    fi

  for proto in "${protocols[@]}"; do
    local ena_dir="${proto}://${ena_path}"

    # Determine which file patterns to try based on detected layout.
    local -a try_files
    if [[ -n "$_ena_layout" ]]; then
        case "$_ena_layout" in
            paired)   try_files=("${srr}_1.fastq.gz" "${srr}_2.fastq.gz") ;;
            single)   try_files=("${srr}.fastq.gz") ;;
            subreads) try_files=("${srr}_subreads.fastq.gz") ;;
        esac
    else
        try_files=("${srr}_1.fastq.gz" "${srr}_2.fastq.gz" "${srr}.fastq.gz" "${srr}_subreads.fastq.gz")
    fi

    local downloaded_any=false
    local got_r1=false
    local got_r2=false

    for fname in "${try_files[@]}"; do
        if [[ "$got_r1" == true && "$got_r2" == true ]]; then
            break
        fi

        local url="${ena_dir}/${fname}"
        local out_file="${target_dir}/${fname}"
        local success=false

        for ((attempt=1; attempt<=max_retries; attempt++)); do
            local wget_err=""
            if wget_err=$(wget --timeout=60 --tries=1 "$url" -O "$out_file" 2>&1); then
                if gzip -t "$out_file" 2>/dev/null; then
                    success=true
                    break
                else
                    echo "  [ENA] Corrupt download for $fname (attempt $attempt/$max_retries)" >&2
                    rm -f "$out_file"
                fi
            else
                local http_code
                http_code=$(echo "$wget_err" | grep -oP 'HTTP request sent.*\K[0-9]{3}' | tail -1)
                if [[ -n "$http_code" && "$http_code" == "404" ]]; then
                    rm -f "$out_file"
                    break
                fi
                rm -f "$out_file"
            fi
            local wait=$(( 5 * (1 << (attempt - 1)) ))
            (( wait > 60 )) && wait=60
            [[ $attempt -lt $max_retries ]] && sleep "$wait"
        done

        if [[ "$success" == true ]]; then
            local base_filename="${fname/${srr}/}"
            base_filename="${base_filename/_subreads.fastq/.fastq}"
            mv "$out_file" "${target_dir}/${rename_prefix}${base_filename}"
            downloaded_any=true
            [[ "$fname" == "${srr}_1.fastq.gz" ]] && got_r1=true
            [[ "$fname" == "${srr}_2.fastq.gz" ]] && got_r2=true
        else
            rm -f "$out_file" 2>/dev/null
        fi
    done

    if [[ "$downloaded_any" == true ]]; then
        if [[ -z "$_ena_protocol" ]]; then
            _ena_protocol="$proto"
        fi
        if [[ -z "$_ena_layout" ]]; then
            if [[ "$got_r1" == true ]]; then
                _ena_layout="paired"
            elif [[ "$fname" == "${srr}_subreads.fastq.gz" ]]; then
                _ena_layout="subreads"
            else
                _ena_layout="single"
            fi
        fi
        return 0
    fi

  done  # end protocol loop

    echo "  [ENA] Failed: $srr (no FASTQ files available via https/ftp)" >&2
    return 1
}

Download_CRR() {
    # Download FASTQ file(s) from CNCB (China National Center for Bioinformation).
    #
    # Args:
    #   $1 = CRR accession
    #   $2 = target directory for FASTQ files
    #   $3 = rename prefix for output files
    #
    # Output: silent on success, errors to stderr.
    # Returns: 0 on success, 1 on failure
    local crr=$1
    local target_dir=$2
    local rename_prefix=$3
    local max_retries=5

    # 1. Map CRR to parent CRA ID (with retry)
    local cra=""
    for ((attempt=1; attempt<=3; attempt++)); do
        cra=$(wget -qO- --timeout=30 --user-agent="Mozilla/5.0" "https://ngdc.cncb.ac.cn/gsa/search?searchTerm=${crr}" | grep -v "example" | grep -oe "CRA[0-9]\+" | uniq | head -n 1)
        [[ -n "$cra" ]] && break
        sleep $(( 5 * attempt ))
    done

    if [[ -z "$cra" ]]; then
        echo "  [CNCB] Failed: Could not map $crr to a CRA project" >&2
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
                        local new_filename=$(echo "$filename" | sed -E 's/_[rR]1([._])/_1\1/; s/_[rR]2([._])/_2\1/')
                        local base_filename=$(echo "$new_filename" | sed "s/${crr}//")
                        mv "$filename" "${target_dir}/${rename_prefix}${base_filename}"
                        success=true
                        break
                    else
                        rm -f "$filename"
                    fi
                fi
            else
                rm -f "$filename"
            fi
            local wait=$(( 5 * (1 << (i - 1)) ))
            (( wait > 60 )) && wait=60
            [[ $i -lt $max_retries ]] && sleep "$wait"
        done

        if [[ "$success" = false ]]; then
            echo "  [CNCB] Failed: $crr/$filename after $max_retries attempts" >&2
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

    # Pre-scan accession file to count and classify accessions
    local has_ncbi_accessions=false
    local has_cncb_accessions=false
    local base_dir="${dir_path%/}"
    local total_accessions=0

    while IFS=$'\t' read -r srr _; do
        [[ -z "$srr" ]] && continue
        total_accessions=$((total_accessions + 1))
        if [[ "$srr" =~ ^CRR ]]; then
            has_cncb_accessions=true
        elif [[ "$srr" =~ ^[EDS]RR ]]; then
            has_ncbi_accessions=true
        fi
    done < "${base_dir}/${acc_file}"

    # Check dependencies
    if [[ "$has_ncbi_accessions" == true || "$has_cncb_accessions" == true ]]; then
        command -v wget >/dev/null 2>&1 || { echo "Error: 'wget' not found (required for data downloads)." >&2; return 1; }
    fi

    # Connectivity check
    if [[ "$has_ncbi_accessions" == true ]]; then
        if ! wget -q --spider --timeout=10 "https://ftp.sra.ebi.ac.uk/" 2>/dev/null \
        && ! wget -q --spider --timeout=10 "ftp://ftp.sra.ebi.ac.uk/" 2>/dev/null; then
            echo "  ENA is unreachable. Check your network connection." >&2
            return 1
        fi
    fi

    # Determine source label for summary
    local source_label="ENA"
    if [[ "$has_cncb_accessions" == true && "$has_ncbi_accessions" == true ]]; then
        source_label="ENA+CNCB"
    elif [[ "$has_cncb_accessions" == true ]]; then
        source_label="CNCB"
    fi

    local fastq_path="${base_dir}/ori_fastq"
    mkdir -p "$fastq_path"

    # Layout detection: probe file type on first SRR, reuse for all subsequent.
    local _ena_layout=""

    # Download all accessions, tracking progress
    local dl_success=0
    local dl_failed=0
    local current=0
    local -a failed_accessions=()

    while IFS=$'\t' read -r srr rename _; do
        [[ -z "$srr" || -z "$rename" ]] && continue
        current=$((current + 1))

        if [[ "$srr" =~ ^CRR ]]; then
            if Download_CRR "$srr" "$fastq_path" "$rename"; then
                dl_success=$((dl_success + 1))
            else
                dl_failed=$((dl_failed + 1))
                failed_accessions+=("$srr")
            fi

        elif [[ "$srr" =~ ^[EDS]RR ]]; then
            if Download_From_ENA "$srr" "$fastq_path" "$rename"; then
                # Silent integrity check — only report failures
                local verify_failed=false
                for fq in "${fastq_path}/${rename}"*.fastq*; do
                    [[ -f "$fq" ]] || continue
                    if ! Verify_Fastq_Integrity "$fq"; then
                        echo "  [verify] Invalid: $(basename "$fq") — removing" >&2
                        rm -f "$fq"
                        verify_failed=true
                    fi
                done
                if [[ "$verify_failed" == true ]]; then
                    local remaining_files
                    remaining_files=$(find "$fastq_path" -name "${rename}*.fastq*" -type f 2>/dev/null | wc -l)
                    if [[ "$remaining_files" -eq 0 ]]; then
                        echo "  [verify] No valid files remaining for $srr" >&2
                        dl_failed=$((dl_failed + 1))
                        failed_accessions+=("$srr")
                        continue
                    fi
                fi
                dl_success=$((dl_success + 1))
            else
                dl_failed=$((dl_failed + 1))
                failed_accessions+=("$srr")
            fi
        else
            echo "  Warning: Unknown accession format: $srr" >&2
        fi

        # Print progress at 25% milestones (for datasets with >= 8 accessions)
        if [[ "$total_accessions" -ge 8 ]]; then
            local pct=$((current * 100 / total_accessions))
            local prev_pct=$(( (current - 1) * 100 / total_accessions ))
            # Print when crossing a 25% boundary
            if [[ $((pct / 25)) -gt $((prev_pct / 25)) ]]; then
                echo "  [download] Progress: ${current}/${total_accessions} (${pct}%)"
            fi
        fi
    done < "${base_dir}/${acc_file}"

    # Orphan file cleanup (silent unless orphans found)
    local orphan_count=0
    for r1 in "${fastq_path}/"*_1.fastq*; do
        [[ -f "$r1" ]] || continue
        local prefix="${r1%_1.fastq*}"
        local has_r2=false
        for r2 in "${prefix}_2.fastq"*; do
            [[ -f "$r2" ]] && has_r2=true && break
        done
        $has_r2 || continue
        for orphan in "${prefix}.fastq" "${prefix}.fastq.gz"; do
            if [[ -f "$orphan" ]]; then
                rm -f "$orphan"
                orphan_count=$((orphan_count + 1))
            fi
        done
    done

    # Count total files and size
    local total_files
    total_files=$(find "$fastq_path" -type f -name '*.fastq*' 2>/dev/null | wc -l | tr -d ' ')
    local total_size
    total_size=$(du -sh "$fastq_path" 2>/dev/null | cut -f1 | tr -d ' ')

    # Detect layout label
    local layout_label="${_ena_layout:-unknown}"

    # Print summary
    if [[ "$dl_failed" -eq 0 ]]; then
        echo "  [${source_label}] Downloaded ${dl_success}/${total_accessions} | Layout: ${layout_label} | ${total_files} files (${total_size})"
    else
        echo "  [${source_label}] Downloaded ${dl_success}/${total_accessions} | ${dl_failed} FAILED | ${total_files} files (${total_size})"
        for acc in "${failed_accessions[@]}"; do
            echo "    Failed: $acc" >&2
        done
        return 1
    fi

    if [[ "$orphan_count" -gt 0 ]]; then
        echo "  Cleaned up $orphan_count orphan file(s)"
    fi

    return 0
}

Amplicon_Common_MakeManifestFileForQiime2() {
    cd "${dataset_path%/}" || { echo "ERROR: Cannot access dataset path: $dataset_path"; exit 1; }
    local temp_file_path="${dataset_path%/}/tmp/temp_file"
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

    local temp_path="${dataset_path%/}/tmp/"
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
    local quality_filter_path="${dataset_path%/}/tmp/step_04_qza_import_QualityFilter/"
    local denoising_path="${dataset_path%/}/tmp/step_05_denoise/"
    local qf_vis_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/"
    local qf_view_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/qf_trim_pos/"
    local qc_vis="${dataset_path%/}/tmp/temp_file/qc_vis/"
    local denoising_vis="${dataset_path%/}/tmp/temp_file/denoising_vis/"
    
    cd "$dataset_path" || { echo "Error: dataset_path not found"; return 1; }
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    quality_filter_path="${dataset_path%/}/tmp/step_04_qza_import_QualityFilter/"
    denoising_path="${dataset_path%/}/tmp/step_05_denoise/"
    qf_trim_pos_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/qf_trim_pos/"
    qc_vis="${dataset_path%/}/tmp/temp_file/qc_vis/"
    denoising_vis="${dataset_path%/}/tmp/temp_file/denoising_vis/"
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

Count_Feature_Table_Reads() {
    # Count total reads in a QIIME2 feature table.
    # Args: $1 = path to feature-table.qza
    # Prints: total read count (integer)
    local table_qza="$1"
    python3 -c "
import tempfile, subprocess, os, sys
from biom import load_table
table_path = sys.argv[1]
with tempfile.TemporaryDirectory() as tmpdir:
    subprocess.run(['qiime', 'tools', 'export', '--input-path', table_path, '--output-path', tmpdir], check=True, capture_output=True)
    table = load_table(os.path.join(tmpdir, 'feature-table.biom'))
    print(int(sum(table.sum(axis='sample'))))
" "$table_qza"
}

Count_Raw_Reads_Total() {
    # Count total raw reads from raw_read_counts.tsv
    # Args: $1 = path to raw_read_counts.tsv
    # Prints: total raw read count (integer)
    awk '{sum+=$3} END{print sum+0}' "$1"
}

################################################################################
#                       ILLUMINA PLATFORM FUNCTIONS                            #
################################################################################
# Functions for Illumina short-read sequencing data processing



Amplicon_Illumina_QualityControlForQZA(){
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local qza_path="${dataset_path%/}/tmp/step_03_qza_import/"
    local quality_filter_path="${dataset_path%/}/tmp/step_04_qza_import_QualityFilter/"
    local qf_vis_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/"
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
    local denoising_path="${base%/}/tmp/step_05_denoise/"
    local qf_vis_path="${base%/}/tmp/temp_file/QualityFilter_vis/"
    local qf_view_path="${base%/}/tmp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${base%/}/tmp/temp_file/QualityFilter_vis/qf_trim_pos/"
    mkdir -p "$denoising_path" "$qf_view_path" "$qf_trim_pos_path" "$qf_vis_path"
    local dataset_name="${base##*/}"

    # Auto-detect input: use quality-filtered data if available (454, Ion Torrent),
    # otherwise use imported data directly. Illumina skips quality-filter to let
    # DADA2 handle quality via its error model, avoiding redundant double-filtering.
    local qza_file=""
    local qf_file="${base%/}/tmp/step_04_qza_import_QualityFilter/${dataset_name}_QualityFilter.qza"
    local import_file="${base%/}/tmp/step_03_qza_import/${dataset_name}.qza"

    if [[ -f "$qf_file" ]]; then
        qza_file="$qf_file"
    else
        qza_file="$import_file"
    fi

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
        # ── SINGLE-END ──
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

        if [[ "${platform:-}" == "ION_TORRENT" ]]; then
            # Ion Torrent: homopolymer-aware error model (denoise-pyro)
            qiime dada2 denoise-pyro \
                --i-demultiplexed-seqs "$qza_file" \
                --p-trunc-len "$end" \
                --p-trim-left "$start" \
                --o-representative-sequences "${denoising_path%/}/${dataset_name}-rep-seqs-denoising.qza" \
                --o-table "${denoising_path%/}/${dataset_name}-table-denoising.qza" \
                --o-denoising-stats "${denoising_path%/}/${dataset_name}-denoising-stats.qza"
        else
            # Illumina: substitution error model (denoise-single)
            qiime dada2 denoise-single \
                --i-demultiplexed-seqs "$qza_file" \
                --p-trunc-len "$end" \
                --p-trim-left "$start" \
                --o-representative-sequences "${denoising_path%/}/${dataset_name}-rep-seqs-denoising.qza" \
                --o-table "${denoising_path%/}/${dataset_name}-table-denoising.qza" \
                --o-denoising-stats "${denoising_path%/}/${dataset_name}-denoising-stats.qza"
        fi
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
    local qza_path="${dataset_path%/}/tmp/step_03_qza_import/"
    local quality_filter_path="${dataset_path%/}/tmp/step_04_qza_import_QualityFilter/"
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
    local quality_filter_path="${dataset_path%/}/tmp/step_04_qza_import_QualityFilter/"
    local dedupicate_path="${dataset_path%/}/tmp/step_05_dedupicate/"
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
    local dedupicate_path="${dataset_path%/}/tmp/step_05_dedupicate/"
    local chimeras_path="${dataset_path%/}/tmp/step_06_ChimerasRemoval/"
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
    local chimeras_path="${dataset_path%/}/tmp/step_06_ChimerasRemoval/"
    local cluster_path="${dataset_path%/}/tmp/step_07_cluster/"
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
    local cluster_path="${dataset_path%/}/tmp/step_07_cluster/"

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
    local qza_path="${dataset_path%/}/tmp/step_03_qza_import/"
    local quality_filter_path="${dataset_path%/}/tmp/step_04_qza_import_QualityFilter/"
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
#                    DEGRADED QUALITY SCORE FUNCTIONS                          #
################################################################################
# Functions for data with binned/unreliable quality scores.
# When DADA2 cannot learn an error model (too few unique Q values),
# this VSEARCH-based pipeline is used instead.
# Reuses LS454 dereplicate/chimera/cluster functions downstream.

Amplicon_DegradedQ_QualityControlForQZA() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local qza_path="${dataset_path%/}/tmp/step_03_qza_import/"
    local quality_filter_path="${dataset_path%/}/tmp/step_04_qza_import_QualityFilter/"
    mkdir -p "$quality_filter_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"

    # Quality scores are unreliable (binned/dummy), so skip q-score filtering.
    # N bases already filtered upstream (max_n=1). Apply length filter only.
    qiime quality-filter q-score \
        --i-demux "${qza_path%/}/${dataset_name}.qza" \
        --p-min-quality 0 \
        --p-min-length-fraction 0.85 \
        --p-max-ambiguous 1 \
        --o-filtered-sequences "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --o-filter-stats "${quality_filter_path%/}/${dataset_name}_filter-stats.qza" \
        --verbose
}

# ── Degraded Quality: Direct VSEARCH derep (bypasses QIIME2) ────────────────
# Replaces: ImportFastqToQiime2 + QualityControlForQZA + LS454_Deduplication + ExportForVsearch
# Reads preprocessed FASTQ directly → vsearch --derep_fulllength → size-annotated FASTA

Amplicon_DegradedQ_DirectDerep() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local vsearch_path="${dataset_path%/}/tmp/step_06_vsearch_cli/"
    local threads="${THREADS_PER_DATASET:-4}"
    mkdir -p "$vsearch_path"

    echo ">>> Direct dereplication via vsearch (bypassing QIIME2)..."
    python3 "${SCRIPTS}/py_16s.py" derep_fastq_for_vsearch \
        --input_dir "$fastq_path" \
        --output_fasta "${vsearch_path%/}/derep_sized.fasta" \
        --threads "$threads"
}

# ── Degraded Quality: VSEARCH CLI pipeline ──────────────────────────────────
# These functions export QIIME2 artifacts → run vsearch natively → import back.
# Pipeline: ExportForVsearch → VsearchDenoise → MapReadsToZotus → ImportResults

Amplicon_DegradedQ_ExportForVsearch() {
    # Export QIIME2 dereplicated artifacts to size-annotated FASTA for vsearch.
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    local dedupicate_path="${dataset_path%/}/tmp/step_05_dedupicate/"
    local vsearch_path="${dataset_path%/}/tmp/step_06_vsearch_cli/"
    mkdir -p "$vsearch_path"

    echo ">>> Exporting dereplicated data for vsearch CLI pipeline..."
    python3 "${SCRIPTS}/py_16s.py" export_derep_for_vsearch \
        --repseq_qza "${dedupicate_path%/}/${dataset_name}-repseq.qza" \
        --table_qza  "${dedupicate_path%/}/${dataset_name}-table.qza" \
        --output_fasta "${vsearch_path%/}/derep_sized.fasta"
}

Amplicon_DegradedQ_VsearchDenoise() {
    # Run native vsearch denoising: 99% pre-cluster → UNOISE3 → uchime3.
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local vsearch_path="${dataset_path%/}/tmp/step_06_vsearch_cli/"
    local threads="${THREADS_PER_DATASET:-4}"

    echo ">>> Step B1: Filtering singletons (minsize=2)..."
    vsearch --sortbysize "${vsearch_path%/}/derep_sized.fasta" \
        --output "${vsearch_path%/}/derep_minsize2.fasta" \
        --minsize 2

    echo ">>> Step B2: 99% pre-clustering..."
    vsearch --cluster_size "${vsearch_path%/}/derep_minsize2.fasta" \
        --id 0.99 \
        --centroids "${vsearch_path%/}/preclust_99.fasta" \
        --sizein --sizeout \
        --threads "$threads"

    echo ">>> Step C: UNOISE3 denoising (minsize=2)..."
    vsearch --cluster_unoise "${vsearch_path%/}/preclust_99.fasta" \
        --centroids "${vsearch_path%/}/zotus.fasta" \
        --sizein --sizeout \
        --minsize 2

    echo ">>> Chimera removal (uchime3_denovo)..."
    vsearch --uchime3_denovo "${vsearch_path%/}/zotus.fasta" \
        --nonchimeras "${vsearch_path%/}/zotus_nochim.fasta" \
        --sizein --sizeout
}

Amplicon_DegradedQ_MapReadsToZotus() {
    # Map all preprocessed reads back to ZOTUs to build the OTU table.
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    local vsearch_path="${dataset_path%/}/tmp/step_06_vsearch_cli/"
    local manifest="${dataset_path%/}/tmp/temp_file/${dataset_name}_manifest.tsv"
    local threads="${THREADS_PER_DATASET:-4}"

    echo ">>> Relabeling reads with sample IDs..."
    python3 "${SCRIPTS}/py_16s.py" relabel_reads_for_mapping \
        --manifest_path "$manifest" \
        --output_fasta "${vsearch_path%/}/all_reads_labeled.fasta" \
        --threads "$threads"

    echo ">>> Mapping reads to ZOTUs (97% identity)..."
    vsearch --usearch_global "${vsearch_path%/}/all_reads_labeled.fasta" \
        --db "${vsearch_path%/}/zotus_nochim.fasta" \
        --id 0.97 \
        --otutabout "${vsearch_path%/}/otu_table.tsv" \
        --sizein \
        --threads "$threads"
}

Amplicon_DegradedQ_ImportResults() {
    # Import vsearch results back into QIIME2 artifacts.
    # Output paths are aligned with Amplicon_LS454_FilterLowFreqOTUs expectations.
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    local vsearch_path="${dataset_path%/}/tmp/step_06_vsearch_cli/"
    local cluster_path="${dataset_path%/}/tmp/step_07_cluster/"
    local manifest="${dataset_path%/}/tmp/temp_file/${dataset_name}_manifest.tsv"
    mkdir -p "$cluster_path"

    echo ">>> Importing vsearch results into QIIME2..."
    python3 "${SCRIPTS}/py_16s.py" import_vsearch_to_qiime2 \
        --zotu_fasta    "${vsearch_path%/}/zotus_nochim.fasta" \
        --otu_table_tsv "${vsearch_path%/}/otu_table.tsv" \
        --manifest_path "$manifest" \
        --output_table_qza  "${cluster_path%/}/${dataset_name}-table-clustered.qza" \
        --output_repseq_qza "${cluster_path%/}/${dataset_name}-repseq-clustered.qza"
}

################################################################################
#                        PACBIO PLATFORM FUNCTIONS                             #
################################################################################
# Functions for PacBio HiFi/CCS long-read sequencing data processing


Amplicon_Pacbio_QualityControlForQZA() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local qza_path="${dataset_path%/}/tmp/step_03_qza_import/"
    local quality_filter_path="${dataset_path%/}/tmp/step_04_qza_import_QualityFilter/"
    local qf_vis_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/"
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
    local quality_filter_path="${base}/tmp/step_04_qza_import_QualityFilter/"
    local denoising_path="${base}/tmp/step_05_denoise/"
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
    local denoising_path="${base}/tmp/step_05_denoise"
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

