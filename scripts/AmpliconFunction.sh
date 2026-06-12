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
#

set -e


################################################################################
#                         COMMON FUNCTIONS                                     #
################################################################################

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
    # Download FASTQ file(s) from ENA using pre-fetched URL map.
    # Requires _ena_url_map file (set by Common_SRADownloadToFastq_MultiSource).
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

    # Look up exact URLs from the filereport map
    local url_line
    url_line=$(awk -v acc="$srr" '$1 == acc {print $2}' "$_ena_url_map")

    if [[ -z "$url_line" ]]; then
        echo "  [ENA] No filereport entry for $srr" >&2
        return 1
    fi

    # url_line is semicolon-separated list of FTP paths from ENA
    local downloaded_any=false

    IFS=';' read -ra urls <<< "$url_line"
    for ftp_path in "${urls[@]}"; do
        [[ -z "$ftp_path" ]] && continue
        local fname
        fname=$(basename "$ftp_path")
        local url="ftp://${ftp_path}"
        local out_file="${target_dir}/${fname}"
        local success=false

        # Prefer paired files; skip single .fastq.gz if we already have _1 + _2
        if [[ "$fname" == "${srr}.fastq.gz" ]]; then
            if [[ -f "${target_dir}/${rename_prefix}_1.fastq.gz" && -f "${target_dir}/${rename_prefix}_2.fastq.gz" ]]; then
                continue
            fi
        fi

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
        else
            rm -f "$out_file" 2>/dev/null
        fi
    done

    if [[ "$downloaded_any" == true ]]; then
        return 0
    fi

    echo "  [ENA] Failed: $srr (download failed for all URLs)" >&2
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


Download_From_NCBI() {
    # Download FASTQ file(s) from NCBI SRA via sra-toolkit (prefetch + fasterq-dump).
    # Used as a fallback when ENA has not mirrored the run's data yet.
    # NCBI imposes per-IP rate limits; callers should not parallelize this.
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

    # Place the SRA workdir alongside the dataset's tmp/, never inside
    # ori_fastq (downstream tooling globs ori_fastq for *.fastq*).
    local tmp_root
    tmp_root="$(dirname "$target_dir")/tmp/sra_dl"
    mkdir -p "$tmp_root"
    local tmp_dir
    tmp_dir=$(mktemp -d "${tmp_root}/${srr}.XXXXXX") || {
        echo "  [NCBI] Failed to create tmp dir for $srr" >&2
        return 1
    }

    # Cushion to ease NCBI per-IP rate limits between sequential accessions.
    # NCBI without an API key allows 3 req/s; 0.3s stays under that.
    sleep 0.3

    local success=false
    local attempt
    for ((attempt=1; attempt<=max_retries; attempt++)); do
        if prefetch --max-size u --output-directory "$tmp_dir" "$srr" >/dev/null 2>&1; then
            # Resolve the .sra path (prefetch layout: <tmp>/<srr>/<srr>.sra).
            local sra_file="${tmp_dir}/${srr}/${srr}.sra"
            [[ -f "$sra_file" ]] || sra_file="${tmp_dir}/${srr}.sra"
            # Pass the file path if we found it, else let fasterq-dump resolve.
            local dump_target="$srr"
            [[ -f "$sra_file" ]] && dump_target="$sra_file"

            if fasterq-dump --split-files --skip-technical --threads 2 \
                            -O "$tmp_dir" "$dump_target" >/dev/null 2>&1; then
                success=true
                break
            fi
            echo "  [NCBI] fasterq-dump failed for $srr (attempt $attempt/$max_retries)" >&2
        else
            echo "  [NCBI] prefetch failed for $srr (attempt $attempt/$max_retries)" >&2
        fi

        local wait=$(( 5 * (1 << (attempt - 1)) ))
        (( wait > 60 )) && wait=60
        [[ $attempt -lt $max_retries ]] && sleep "$wait"
    done

    if ! $success; then
        rm -rf "$tmp_dir"
        echo "  [NCBI] Failed: $srr after $max_retries attempts" >&2
        return 1
    fi

    # Move + gzip + rename FASTQ outputs. fasterq-dump emits ${srr}.fastq (SE),
    # ${srr}_1.fastq + ${srr}_2.fastq (PE), or ${srr}_subreads.fastq (PacBio).
    local downloaded_any=false
    shopt -s nullglob
    for fq in "${tmp_dir}/${srr}"*.fastq; do
        [[ -f "$fq" ]] || continue
        local fname
        fname=$(basename "$fq")
        local stem="${fname%.fastq}"
        local suffix="${stem#${srr}}"
        # Match Download_From_ENA: PacBio _subreads is folded into the base name.
        suffix="${suffix/_subreads/}"
        local out_path="${target_dir}/${rename_prefix}${suffix}.fastq.gz"
        if gzip -c "$fq" > "$out_path"; then
            downloaded_any=true
        else
            rm -f "$out_path"
        fi
        rm -f "$fq"
    done
    shopt -u nullglob

    rm -rf "$tmp_dir"

    if ! $downloaded_any; then
        echo "  [NCBI] Failed: $srr (no FASTQ output produced)" >&2
        return 1
    fi
    return 0
}

Common_SRADownloadToFastq_MultiSource() {
    local dir_path="" acc_file="" bioproject=""
    OPTIND=1
    while getopts ":d:a:b:" opt; do
        case $opt in
            d) dir_path=$OPTARG ;;
            a) acc_file=$OPTARG ;;
            b) bioproject=$OPTARG ;;
            ?) echo "Unknown Parameter: -$OPTARG" >&2; return 1 ;;
        esac
    done

    if [[ -z "$dir_path" || -z "$acc_file" ]]; then
        echo "Usage: Common_SRADownloadToFastq_MultiSource -d <dir> -a <accession_tsv> [-b <bioproject>]" >&2
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

    # ENA reachability is informational: if ENA is down we can still serve
    # NCBI accessions via the SRA-Toolkit fallback below.
    local ena_reachable=true
    if [[ "$has_ncbi_accessions" == true ]]; then
        if ! wget -q --spider --timeout=10 "https://www.ebi.ac.uk/ena/portal/api/" 2>/dev/null; then
            ena_reachable=false
            echo "  [ENA] API unreachable; will route NCBI accessions through SRA Toolkit." >&2
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

    # Fetch ENA filereport to get exact download URLs for NCBI accessions.
    # Skipped when ENA is unreachable — gap detection below will see an empty
    # map and route the whole dataset through the NCBI SRA Toolkit.
    local _ena_url_map="${base_dir}/.ena_url_map.tsv"
    if [[ "$has_ncbi_accessions" == true && "$ena_reachable" == true && -n "$bioproject" ]]; then
        echo "  [ENA] Fetching filereport for ${bioproject}..."
        local filereport_url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${bioproject}&result=read_run&fields=run_accession,fastq_ftp&format=tsv"
        if wget -q --timeout=30 "$filereport_url" -O "${_ena_url_map}.raw" 2>/dev/null; then
            # Extract run_accession and fastq_ftp columns, skip header
            awk -F'\t' 'NR>1 && $1!="" && $2!="" {print $1"\t"$2}' "${_ena_url_map}.raw" > "$_ena_url_map"
            rm -f "${_ena_url_map}.raw"
            local map_count
            map_count=$(wc -l < "$_ena_url_map" | tr -d ' ')
            echo "  [ENA] Got URLs for ${map_count} runs"
        else
            echo "  [ENA] Filereport query failed; falling back to NCBI SRA Toolkit." >&2
            : > "$_ena_url_map"
        fi
    else
        : > "$_ena_url_map"
    fi

    # Detect ENA mirroring gap: any input [EDS]RR accession lacking a populated
    # fastq_ftp URL in ENA's filereport. ENA stores the metadata record before
    # mirroring the actual reads from NCBI SRA, so a run can exist on NCBI with
    # full data while ENA reports read_count=0 and empty FTP fields.
    # When this happens, we cannot recover by retrying ENA — switch the entire
    # dataset to NCBI SRA Toolkit and discard any partial ENA downloads to keep
    # the dataset's source uniform.
    local ncbi_fallback_mode=false
    if [[ "$has_ncbi_accessions" == true ]]; then
        local ena_gap_count=0
        while IFS=$'\t' read -r srr _; do
            [[ -z "$srr" ]] && continue
            [[ ! "$srr" =~ ^[EDS]RR ]] && continue
            if ! awk -v acc="$srr" '$1 == acc && $2 != "" {found=1; exit} END{exit !found}' "$_ena_url_map" 2>/dev/null; then
                ena_gap_count=$((ena_gap_count + 1))
            fi
        done < "${base_dir}/${acc_file}"

        if (( ena_gap_count > 0 )); then
            if ! command -v prefetch >/dev/null 2>&1 || ! command -v fasterq-dump >/dev/null 2>&1; then
                echo "  [ENA] ${ena_gap_count} accession(s) lack mirrored FASTQ in ENA, but NCBI SRA Toolkit (prefetch, fasterq-dump) is not available." >&2
                echo "  [ENA] Install sra-toolkit to enable NCBI fallback for unmirrored runs." >&2
                return 1
            fi
            # Fail fast if NCBI is unreachable rather than burning hours on
            # 5 backoff retries per accession before discovering the network is down.
            if ! wget -q --spider --timeout=10 "https://trace.ncbi.nlm.nih.gov/" 2>/dev/null; then
                echo "  [ENA] ${ena_gap_count} accession(s) need NCBI fallback, but trace.ncbi.nlm.nih.gov is unreachable." >&2
                return 1
            fi
            echo "  [ENA] ${ena_gap_count} accession(s) lack mirrored FASTQ in ENA — switching dataset to NCBI SRA Toolkit"
            # Wipe any partial ENA downloads so the dataset has a single source.
            find "$fastq_path" -mindepth 1 -delete 2>/dev/null || true
            ncbi_fallback_mode=true
            source_label="NCBI"
        fi
    fi

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
            local _dl_ok=false
            if $ncbi_fallback_mode; then
                Download_From_NCBI "$srr" "$fastq_path" "$rename" && _dl_ok=true
            else
                Download_From_ENA "$srr" "$fastq_path" "$rename" && _dl_ok=true
            fi
            if $_dl_ok; then
                # Silent integrity check — only report failures.
                # Enumerate exact suffixes to avoid prefix collisions
                # (e.g. rename="MNB_04" matching files of "MNB_04317").
                local _candidate_files=(
                    "${fastq_path}/${rename}.fastq"
                    "${fastq_path}/${rename}.fastq.gz"
                    "${fastq_path}/${rename}_1.fastq"
                    "${fastq_path}/${rename}_1.fastq.gz"
                    "${fastq_path}/${rename}_2.fastq"
                    "${fastq_path}/${rename}_2.fastq.gz"
                )
                local verify_failed=false
                for fq in "${_candidate_files[@]}"; do
                    [[ -f "$fq" ]] || continue
                    if ! Verify_Fastq_Integrity "$fq"; then
                        echo "  [verify] Invalid: $(basename "$fq") — removing" >&2
                        rm -f "$fq"
                        verify_failed=true
                    fi
                done
                if [[ "$verify_failed" == true ]]; then
                    local remaining_files=0
                    for fq in "${_candidate_files[@]}"; do
                        [[ -f "$fq" ]] && remaining_files=$((remaining_files + 1))
                    done
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

    # Retry failed accessions after a cooldown period.
    # ENA's FTP enforces a per-host connection limit so we wait 120s; NCBI's
    # SRA endpoints don't, so 30s is enough when running in NCBI fallback mode.
    if [[ ${#failed_accessions[@]} -gt 0 ]]; then
        local cooldown=120
        $ncbi_fallback_mode && cooldown=30
        echo "  [${source_label}] Retrying ${#failed_accessions[@]} failed sample(s) in ${cooldown}s..."
        sleep "$cooldown"
        local -a still_failed=()
        for failed_srr in "${failed_accessions[@]}"; do
            local failed_rename
            failed_rename=$(awk -F'\t' -v acc="$failed_srr" '$1 == acc {print $2}' "${base_dir}/${acc_file}")
            local retry_ok=false
            if [[ -n "$failed_rename" ]]; then
                if $ncbi_fallback_mode; then
                    Download_From_NCBI "$failed_srr" "$fastq_path" "$failed_rename" && retry_ok=true
                else
                    Download_From_ENA "$failed_srr" "$fastq_path" "$failed_rename" && retry_ok=true
                fi
            fi
            if $retry_ok; then
                echo "  [${source_label}] Retry OK: $failed_srr"
                dl_success=$((dl_success + 1))
                dl_failed=$((dl_failed - 1))
            else
                still_failed+=("$failed_srr")
            fi
        done
        failed_accessions=("${still_failed[@]}")
    fi

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

    # Print summary
    if [[ "$dl_failed" -eq 0 ]]; then
        echo "  [${source_label}] Downloaded ${dl_success}/${total_accessions} | ${total_files} files (${total_size})"
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

    rm -f "$_ena_url_map"
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
    local qf_trim_pos_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/qf_trim_pos/"
    local qc_vis="${dataset_path%/}/tmp/temp_file/qc_vis/"
    local denoising_vis="${dataset_path%/}/tmp/temp_file/denoising_vis/"

    cd "$dataset_path" || { echo "Error: dataset_path not found"; return 1; }
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    mkdir -p "$qc_vis" "$denoising_vis"
    
    # Check for denoising output
    if [ -d "$denoising_path" ] && [ -f "${denoising_path%/}/${dataset_name}-table-denoising.qza" ]; then
        cp "${denoising_path%/}/${dataset_name}-rep-seqs-denoising.qza" "${dataset_path%/}/${dataset_name}-${MODE}-final-rep-seqs.qza"
        cp "${denoising_path%/}/${dataset_name}-table-denoising.qza" "${dataset_path%/}/${dataset_name}-${MODE}-final-table.qza"
        
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
        [[ "${KEEP_INTERMEDIATE:-0}" == "1" ]] || find "$dataset_path" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +

        rm -f "${dataset_path%/}/"{denoising.log,fastp.html,fastp.json}

        rm -rf "${dataset_path%/}/ori_fastq" 2>/dev/null || true
        rm -rf "${dataset_path%/}/working_fastq" 2>/dev/null || true

        return 0

    # Check for vsearch output (454 pipeline)
    elif [ -f "${dataset_path%/}/${dataset_name}-table-vsearch.qza" ]; then
        mv "${dataset_path%/}/${dataset_name}-table-vsearch.qza" "${dataset_path%/}/${dataset_name}-${MODE}-final-table.qza"
        mv "${dataset_path%/}/${dataset_name}-rep-seqs-vsearch.qza" "${dataset_path%/}/${dataset_name}-${MODE}-final-rep-seqs.qza"
        [[ "${KEEP_INTERMEDIATE:-0}" == "1" ]] || find "$dataset_path" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
        
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
            --o-denoising-stats "${denoising_path%/}/${dataset_name}-denoising-stats.qza" \
            --p-n-threads "$cpu"
    else
        # ── SINGLE-END ──
        if [[ "${platform:-}" == "ION_TORRENT" ]]; then
            # Ion Torrent: variable-length reads → skip trim_pos_deblur.
            # denoise-pyro handles quality internally via --p-trunc-q and --p-max-ee.
            # trunc-len 0 = no fixed truncation (appropriate for variable-length data).
            local start="${start_in:-0}"
            echo "$start 0" > "${qf_trim_pos_path%/}/Trim_position.txt"

            qiime dada2 denoise-pyro \
                --i-demultiplexed-seqs "$qza_file" \
                --p-trunc-len 0 \
                --p-trim-left "$start" \
                --o-representative-sequences "${denoising_path%/}/${dataset_name}-rep-seqs-denoising.qza" \
                --o-table "${denoising_path%/}/${dataset_name}-table-denoising.qza" \
                --o-denoising-stats "${denoising_path%/}/${dataset_name}-denoising-stats.qza" \
                --p-n-threads "$cpu"
        else
            # Illumina: compute trim positions from QC visualization
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

            # Illumina: substitution error model (denoise-single)
            qiime dada2 denoise-single \
                --i-demultiplexed-seqs "$qza_file" \
                --p-trunc-len "$end" \
                --p-trim-left "$start" \
                --o-representative-sequences "${denoising_path%/}/${dataset_name}-rep-seqs-denoising.qza" \
                --o-table "${denoising_path%/}/${dataset_name}-table-denoising.qza" \
                --o-denoising-stats "${denoising_path%/}/${dataset_name}-denoising-stats.qza" \
                --p-n-threads "$cpu"
        fi
    fi
}
################################################################################
#                        LS454 PLATFORM FUNCTIONS                              #
################################################################################
# Functions for 454 pyrosequencing data processing
# 454 pipeline: QC (length filter) → Deduplication → Chimera removal → OTU clustering (97%) → Filter low-freq OTUs
# Quality scores from fasterq-dump are unreliable for 454, so q-score filtering is disabled.

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
    local _derep_stdout
    _derep_stdout=$(python3 "${SCRIPTS}/py_16s.py" derep_fastq_for_vsearch \
        --input_dir "$fastq_path" \
        --output_fasta "${vsearch_path%/}/derep_sized.fasta" \
        --threads "$threads")
    local _derep_rc=$?
    echo "$_derep_stdout"
    [[ $_derep_rc -ne 0 ]] && return $_derep_rc

    DEREP_INPUT_READS=$(printf '%s\n' "$_derep_stdout" | awk -F= '/^DEREP_INPUT_READS=/{print $2; exit}')
    DEREP_UNIQUE_SEQS=$(printf '%s\n' "$_derep_stdout" | awk -F= '/^DEREP_UNIQUE_SEQS=/{print $2; exit}')
    export DEREP_INPUT_READS DEREP_UNIQUE_SEQS
}

# ── Degraded Quality: Abundance sanity check ────────────────────────────────
# Detects datasets that lack the abundance signal de novo ASV/OTU calling
# requires. Trigger: single-sample dataset where ≥95% of reads are singletons
# after dereplication. In that regime, UNOISE3/VSEARCH minsize=2 wipes the
# data to zero (every read looks unique because of sequencing error), and there
# are no other samples to provide cross-sample replication. Any downstream
# table would be statistically meaningless, so abort early.
#
# Returns 0 if data is OK to process, 99 if dataset should be skipped.
Amplicon_DegradedQ_AbundanceSanityCheck() {
    local n_samples="${1:-0}"
    local input_reads="${DEREP_INPUT_READS:-0}"
    local unique_seqs="${DEREP_UNIQUE_SEQS:-0}"

    [[ "$input_reads" -le 0 ]] && return 0

    # parts-per-1000 to avoid floating-point in bash
    local ratio_ppt
    ratio_ppt=$(python3 -c "print(int(${unique_seqs} * 1000 / ${input_reads}))")

    if [[ "$n_samples" -eq 1 ]] && [[ "$ratio_ppt" -gt 950 ]]; then
        local ratio_pct
        ratio_pct=$(python3 -c "print(f'{${unique_seqs}/${input_reads}*100:.2f}')")
        echo ""
        echo "============================================================"
        echo "UNTRUSTWORTHY DATA DETECTED — aborting analysis"
        echo "============================================================"
        echo "  Samples in dataset:   1"
        echo "  Input reads:          ${input_reads}"
        echo "  Unique sequences:     ${unique_seqs}"
        echo "  Unique ratio:         ${ratio_pct}%"
        echo ""
        echo "  Reason: single-sample dataset where >=95% of reads are singletons."
        echo "  De novo ASV/OTU calling needs duplicate reads to separate signal"
        echo "  from sequencing error; with no abundance signal, UNOISE3 has"
        echo "  nothing to denoise, and there is no cross-sample replication to"
        echo "  fall back on. Any downstream taxonomy table would be"
        echo "  statistically meaningless."
        echo ""
        echo "  Skipping this dataset and moving to the next BioProject."
        echo "============================================================"
        return 99
    fi
    return 0
}

# ── Degraded Quality: VSEARCH CLI pipeline ──────────────────────────────────
# These functions export QIIME2 artifacts → run vsearch natively → import back.
# Pipeline: ExportForVsearch → VsearchDenoise → MapReadsToZotus → ImportResults

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

    # --strand follows $OTU_STRAND (plus for short reads; both for PacBio).
    # On --strand both the 99% precluster merges forward/reverse-complement
    # copies (summing --sizein) BEFORE UNOISE3's minsize=2 survival decision,
    # so denoising sees the correct total abundance instead of a split count.
    echo ">>> Step B2: 99% pre-clustering (strand=${OTU_STRAND:-plus})..."
    vsearch --cluster_size "${vsearch_path%/}/derep_minsize2.fasta" \
        --id 0.99 \
        --strand "${OTU_STRAND:-plus}" \
        --centroids "${vsearch_path%/}/preclust_99.fasta" \
        --sizein --sizeout \
        --threads "$threads"

    echo ">>> Step C: UNOISE3 denoising (minsize=2, strand=${OTU_STRAND:-plus})..."
    vsearch --cluster_unoise "${vsearch_path%/}/preclust_99.fasta" \
        --centroids "${vsearch_path%/}/zotus.fasta" \
        --strand "${OTU_STRAND:-plus}" \
        --sizein --sizeout \
        --minsize 2

    echo ">>> Chimera removal (uchime3_denovo)..."
    vsearch --uchime3_denovo "${vsearch_path%/}/zotus.fasta" \
        --nonchimeras "${vsearch_path%/}/zotus_nochim.fasta" \
        --sizein --sizeout
}

################################################################################
#                     UNIFIED OTU BACK-END (vsearch)                           #
################################################################################
# Shared genus-level OTU back-end used by --otu across platforms.
# Chain (short-read / pooled):
#   <platform preprocess> -> DirectDerep -> AbundanceSanityCheck
#     -> DegradedQ_VsearchDenoise (99% precluster -> UNOISE3 -> uchime3)
#     -> Amplicon_OTU_ClusterFast97  (B1: collapse ZOTUs to 97% OTU centroids)
#     -> Amplicon_OTU_MapBack        (B3: strand-aware read mapping -> OTU table)
#     -> Amplicon_OTU_ImportResults  -> LS454_FilterLowFreqOTUs -> FinalFilesCleaning
# Strand is controlled by $OTU_STRAND (plus for short reads, both for long reads).

Amplicon_OTU_ClusterFast97() {
    # B1: collapse UNOISE3 ZOTUs (zotus_nochim.fasta) to genus-level 97% OTU
    # centroids. Without this, rep-seqs stay at ASV granularity (a sequence
    # cloud) and only collapse at the abundance-table level — which reintroduces
    # the ASV cloud when merging across datasets. cluster_fast sorts by length.
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local vsearch_path="${dataset_path%/}/tmp/step_06_vsearch_cli/"
    local threads="${THREADS_PER_DATASET:-4}"

    echo ">>> Clustering ZOTUs to 97% OTUs (cluster_fast, strand=${OTU_STRAND:-plus})..."
    vsearch --cluster_fast "${vsearch_path%/}/zotus_nochim.fasta" \
        --id 0.97 \
        --centroids "${vsearch_path%/}/otus_97.fasta" \
        --sizein --sizeout \
        --strand "${OTU_STRAND:-plus}" \
        --threads "$threads"
}

Amplicon_OTU_MapBack() {
    # Map all preprocessed reads back to the 97% OTUs to build the OTU table.
    # Strand-aware (B3): short reads = plus; PacBio/ONT = both.
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

    echo ">>> Mapping reads to 97% OTUs (id=0.97, strand=${OTU_STRAND:-plus})..."
    vsearch --usearch_global "${vsearch_path%/}/all_reads_labeled.fasta" \
        --db "${vsearch_path%/}/otus_97.fasta" \
        --id 0.97 \
        --strand "${OTU_STRAND:-plus}" \
        --otutabout "${vsearch_path%/}/otu_table.tsv" \
        --sizein \
        --threads "$threads"
}

Amplicon_OTU_ImportResults() {
    # Import the 97% OTU rep-seqs + OTU table back into QIIME2 artifacts.
    # Output paths align with Amplicon_LS454_FilterLowFreqOTUs expectations.
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    local vsearch_path="${dataset_path%/}/tmp/step_06_vsearch_cli/"
    local cluster_path="${dataset_path%/}/tmp/step_07_cluster/"
    local manifest="${dataset_path%/}/tmp/temp_file/${dataset_name}_manifest.tsv"
    mkdir -p "$cluster_path"

    echo ">>> Importing 97% OTU results into QIIME2..."
    python3 "${SCRIPTS}/py_16s.py" import_vsearch_to_qiime2 \
        --zotu_fasta    "${vsearch_path%/}/otus_97.fasta" \
        --otu_table_tsv "${vsearch_path%/}/otu_table.tsv" \
        --manifest_path "$manifest" \
        --output_table_qza  "${cluster_path%/}/${dataset_name}-table-clustered.qza" \
        --output_repseq_qza "${cluster_path%/}/${dataset_name}-repseq-clustered.qza"
}

Amplicon_Illumina_OTU_Preprocess() {
    # Produce clean, full-amplicon single-end reads for the pooled OTU back-end.
    #   normal quality : per-sample PE merge (vsearch --fastq_mergepairs);
    #                    if merge rate < OTU_MERGE_MIN, fall back to forward-only
    #                    R1; then maxee filter (vsearch --fastq_filter).
    #   binned quality : reuse degraded_quality_preprocess (5' trim + truncate +
    #                    N filter, forward-only) — maxee is unreliable on binned Q.
    # On return: $fastq_path -> clean dir; sequence_type=single.
    # Inputs from caller scope: fastp_path, quality_status, sequence_type.
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local clean_path="${dataset_path%/}/tmp/step_02d_otu_preprocess"
    local threads="${THREADS_PER_DATASET:-4}"
    local maxee="${OTU_MAXEE:-1.0}"
    local merge_min="${OTU_MERGE_MIN:-0.5}"

    rm -rf "$clean_path"
    mkdir -p "$clean_path"

    if [[ "$quality_status" == "degraded_binned" ]]; then
        echo ">>> OTU preprocess (binned/degraded quality): forward-only truncation path..."
        python3 "${SCRIPTS}/py_16s.py" degraded_quality_preprocess \
            --input_dir "$fastp_path" \
            --output_dir "$clean_path" \
            --trim_front 15 --truncate_length 0 \
            --max_n 1 --min_length 50 \
            --sequence_type "$sequence_type" \
            --threads "$threads"
    else
        echo ">>> OTU preprocess (normal quality): merge pairs + maxee ${maxee}..."
        local r1 r2 sample merged n_in n_merged
        shopt -s nullglob
        if [[ "$sequence_type" == "paired" ]]; then
            for r1 in "${fastp_path}/"*_1.fastq*; do
                [[ -f "$r1" ]] || continue
                r2="${r1/_1.fastq/_2.fastq}"
                sample=$(basename "$r1"); sample="${sample%%_1.fastq*}"
                if [[ -f "$r2" ]]; then
                    merged="${clean_path}/${sample}_merged.fastq"
                    if vsearch --fastq_mergepairs "$r1" --reverse "$r2" \
                               --fastq_qmax 93 --fastq_qmaxout 93 \
                               --fastqout "$merged" --threads "$threads" --quiet 2>/dev/null \
                       && [[ -s "$merged" ]]; then
                        n_in=$(( $(zcat -f "$r1" | wc -l) / 4 ))
                        n_merged=$(( $(wc -l < "$merged") / 4 ))
                        if python3 -c "import sys; sys.exit(0 if ($n_in>0 and $n_merged/$n_in >= $merge_min) else 1)"; then
                            vsearch --fastq_filter "$merged" --fastq_maxee "$maxee" --fastq_qmax 93 \
                                    --fastqout "${clean_path}/${sample}.fastq" --threads "$threads" --quiet
                            rm -f "$merged"
                            continue
                        fi
                        echo "    ${sample}: low merge rate (${n_merged}/${n_in}) → forward-only R1"
                        rm -f "$merged"
                    else
                        echo "    ${sample}: merge failed → forward-only R1"
                        rm -f "$merged" 2>/dev/null || true
                    fi
                fi
                # Fallback (low merge rate / merge failed / no R2): forward-only R1
                vsearch --fastq_filter "$r1" --fastq_maxee "$maxee" --fastq_qmax 93 \
                        --fastqout "${clean_path}/${sample}.fastq" --threads "$threads" --quiet
            done
        else
            for r1 in "${fastp_path}/"*.fastq*; do
                [[ -f "$r1" ]] || continue
                sample=$(basename "$r1"); sample="${sample%%.fastq*}"
                vsearch --fastq_filter "$r1" --fastq_maxee "$maxee" --fastq_qmax 93 \
                        --fastqout "${clean_path}/${sample}.fastq" --threads "$threads" --quiet
            done
        fi
        shopt -u nullglob
    fi

    fastq_path="$clean_path"
    sequence_type="single"
    export fastq_path sequence_type
}

Amplicon_OTU_RunPooledChain() {
    # Shared pooled OTU back-end (Illumina / 454 / Ion / PacBio).
    #   derep → UNOISE3(99% precluster) → uchime3 → cluster_fast 97%
    #   → strand-aware map-back → import → low-freq filter → finalize.
    # ONT does NOT use this — it keeps its own racon-polished, error-tolerant
    # back-half (map at id=0.90; see Amplicon_ONT_*).
    # Caller must have set (in scope): fastq_path (clean reads dir),
    #   sequence_type=single, OTU_STRAND (plus|both), n_srr (sample count),
    #   and have $dataset_ID / $low_quality_log available (run.sh scope).
    Amplicon_Common_MakeManifestFileForQiime2
    Amplicon_DegradedQ_DirectDerep

    # Untrustworthy single-sample, all-singleton data → skip (rc 99).
    if ! Amplicon_DegradedQ_AbundanceSanityCheck "${n_srr:-0}"; then
        local untrust_marker="${dataset_path%/}/${dataset_ID}-UNTRUSTABLE.txt"
        {
            echo "Dataset:           ${dataset_ID}"
            echo "Reason:            single sample with >=95% singleton reads"
            echo "Samples:           ${n_srr}"
            echo "Input reads:       ${DEREP_INPUT_READS}"
            echo "Unique sequences:  ${DEREP_UNIQUE_SEQS}"
            echo "Detected at:      $(date '+%Y-%m-%d %H:%M:%S')"
        } > "$untrust_marker"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $dataset_ID - UNTRUSTABLE - 1 sample, ${DEREP_UNIQUE_SEQS}/${DEREP_INPUT_READS} unique" >> "$low_quality_log"
        exit 99
    fi

    Amplicon_DegradedQ_VsearchDenoise   # 99% precluster → UNOISE3 → uchime3
    Amplicon_OTU_ClusterFast97          # B1: collapse ZOTUs → 97% OTUs
    Amplicon_OTU_MapBack                # B3: strand-aware read mapping
    Amplicon_OTU_ImportResults
    Amplicon_LS454_FilterLowFreqOTUs
    Amplicon_Common_FinalFilesCleaning
}

Amplicon_IonTorrent_OTU_Preprocess() {
    # Ion Torrent OTU preprocess: strip 5' 10 bp (signal instability — the
    # vsearch equivalent of DADA2 denoise-pyro --p-trim-left 10) + maxee filter.
    # No fixed truncation (B6). SE reads. Reads input from $fastp_path.
    # On return: fastq_path -> clean dir; sequence_type=single.
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local in_dir="${fastp_path:-${dataset_path%/}/tmp/step_02_fastp}"
    local clean_path="${dataset_path%/}/tmp/step_02d_otu_preprocess"
    local threads="${THREADS_PER_DATASET:-4}"
    # Ion Torrent reads (~200-400 bp) carry homopolymer indel errors, so a flat
    # maxee 1.0 (≤1 expected error / read) discards the bulk of real reads
    # (~6% retention observed). Use a more permissive Ion-specific default (2.0);
    # tune via ION_OTU_MAXEE. (Illumina keeps OTU_MAXEE=1.0.)
    local maxee="${ION_OTU_MAXEE:-2.0}"
    local stripleft="${ION_OTU_STRIPLEFT:-10}"

    rm -rf "$clean_path"; mkdir -p "$clean_path"
    echo ">>> Ion OTU preprocess: stripleft ${stripleft} + maxee ${maxee}..."
    local fq sample
    shopt -s nullglob
    for fq in "${in_dir}/"*.fastq*; do
        [[ -f "$fq" ]] || continue
        sample=$(basename "$fq"); sample="${sample%%.fastq*}"
        vsearch --fastq_filter "$fq" \
            --fastq_stripleft "$stripleft" \
            --fastq_maxee "$maxee" \
            --fastq_qmax 93 \
            --fastqout "${clean_path}/${sample}.fastq" \
            --threads "$threads" --quiet
    done
    shopt -u nullglob

    fastq_path="$clean_path"; sequence_type="single"; export fastq_path sequence_type
}

Amplicon_Pacbio_OTU_Preprocess() {
    # PacBio CCS OTU preprocess: length window + per-base-rate maxee (B5).
    # CCS reads are high-accuracy full-length 16S (~1500 bp). Keep a length
    # window (no truncation, B4); filter by error RATE (flat maxee too strict at
    # 1500 bp); raise fastq_qmax (CCS Q can exceed vsearch's default 41). Reads
    # occur in both orientations → caller sets OTU_STRAND=both. Input from
    # $adapter_removed_path (PacBio has no separate primer-trim step).
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local in_dir="${adapter_removed_path:-${dataset_path%/}/tmp/step_01_adapter_removed}"
    local clean_path="${dataset_path%/}/tmp/step_02d_otu_preprocess"
    local threads="${THREADS_PER_DATASET:-4}"
    local maxee_rate="${OTU_MAXEE_RATE:-0.01}"
    local minlen="${PACBIO_OTU_MINLEN:-1000}"
    local maxlen="${PACBIO_OTU_MAXLEN:-2000}"

    rm -rf "$clean_path"; mkdir -p "$clean_path"
    echo ">>> PacBio OTU preprocess: len[${minlen},${maxlen}] + maxee_rate ${maxee_rate}..."
    local fq sample
    shopt -s nullglob
    for fq in "${in_dir}/"*.fastq*; do
        [[ -f "$fq" ]] || continue
        case "$(basename "$fq")" in fastp.*) continue ;; esac   # skip fastp.json/html
        sample=$(basename "$fq"); sample="${sample%%.fastq*}"
        vsearch --fastq_filter "$fq" \
            --fastq_qmax 93 \
            --fastq_minlen "$minlen" --fastq_maxlen "$maxlen" \
            --fastq_maxee_rate "$maxee_rate" \
            --fastqout "${clean_path}/${sample}.fastq" \
            --threads "$threads" --quiet
    done
    shopt -u nullglob

    fastq_path="$clean_path"; sequence_type="single"; export fastq_path sequence_type
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
        --p-n-threads "$cpu"
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
# Oxford Nanopore long-read amplicon processing.
#
# Faithful port of the ONT-AmpSeq algorithm (Eskildsen et al.) into the
# Meta2Data framework, with two deliberate adaptations:
#
#   1. Length window is AUTO-DETECTED (read-length peak finding via
#      py_16s.py detect_length_window) instead of the user-set 1200-1600bp
#      ONT-AmpSeq requires. Primers are auto-trimmed upstream by
#      entropy_primer_detect.py. This keeps the pipeline fully automatic and
#      amplicon-agnostic (16S / ITS / 18S all auto-fit their own peak).
#
#   2. The abundance table is built by MAPPING ALL READS BACK to the OTU
#      centroids (vsearch --usearch_global, as in the DegradedQ pipeline) so
#      it carries true per-sample read abundance. ONT-AmpSeq's
#      `cluster_fast --otutabout` instead reports polished-centroid counts,
#      which destroy the abundance signal and are not comparable with the
#      read-abundance feature tables every other Meta2Data platform produces.
#
# Algorithm:
#   chopper length/quality filter (auto window)
#     -> per-sample fasta + vsearch --cluster_unoise --minsize 1 (light denoise)
#     -> concatenate per-sample centroids
#     -> per-sample minimap2 map-ont + racon consensus polishing (ONT error fix)
#     -> per-sample vsearch --sortbysize --sample <id> + merge
#     -> vsearch --cluster_fast --id <ONT_OTU_IDENTITY> --relabel OTU_  (OTU seqs)
#     -> map all reads back to OTUs at <ONT_MAP_IDENTITY> (ONT-error-tolerant,
#        default 0.90) -> true read-abundance table
#     -> import to QIIME2 (reuses import_vsearch_to_qiime2)
#     -> Amplicon_LS454_FilterLowFreqOTUs + Amplicon_Common_FinalFilesCleaning
#
# DADA2 is intentionally not used: ONT error profiles are unsupported by
# DADA2's error model. Quality is handled by chopper + abundance/consensus
# methods (UNOISE3 + racon), mirroring the DegradedQ rationale.
#
# Tunable env vars (all optional, sensible defaults):
#   ONT_QUALITY            chopper mean-quality cutoff           (default 20)
#   ONT_LENGTH_TOLERANCE   length window half-width fraction     (default 0.15)
#   ONT_LENGTH_FLOOR       hard lower length bound, bp           (default 200)
#   ONT_OTU_IDENTITY       final cluster_fast (OTU) identity     (default 0.97)
#   ONT_MAP_IDENTITY       read->OTU mapping identity            (default 0.90)
#                          (lower than OTU identity to tolerate raw ONT error)

# ── Per-sample length-window + quality filter via chopper ───────────────────
# In : $fastq_path  (primer-trimmed reads, one FASTQ per sample, stem = sample id)
# Out: tmp/step_03_chopper/<stem>.fastq  (+ .chopper_done marker; sets ONT_FASTQ_DIR)
Amplicon_ONT_ChopperFilter() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path" || { echo "❌ [ONT] cannot cd $dataset_path"; return 1; }
    trimmed_path="${dataset_path%/}"; dataset_name="${trimmed_path##*/}"
    local in_dir="${fastq_path%/}"
    local out_dir="${dataset_path%/}/tmp/step_03_chopper"
    local threads="${THREADS_PER_DATASET:-${cpu:-4}}"
    # chopper's own docs default -q to 0 and use -q 10 as the practical example;
    # ONT-AmpSeq's config uses Q20, but that is tuned for high-quality (R10/Zymo)
    # data and discards ~99% of typical public ONT reads (mean Q often <20, e.g.
    # R9). Default to Q10 (chopper's documented practical value) for general
    # public-data reanalysis; residual error is handled by racon polish + 97%
    # clustering. Override per-run via ONT_QUALITY (e.g. =20 to match ONT-AmpSeq).
    local qual="${ONT_QUALITY:-10}"
    local tol="${ONT_LENGTH_TOLERANCE:-0.15}"
    local floor="${ONT_LENGTH_FLOOR:-200}"
    mkdir -p "$out_dir"

    echo ">>> [ONT] Auto-detecting amplicon length window..."
    local win lo hi peak
    win=$(python3 "${SCRIPTS}/py_16s.py" detect_length_window \
        --input_dir "$in_dir" --tolerance "$tol" --floor "$floor" \
        --max_sample_reads 10000) || { echo "❌ [ONT] length window detection failed"; return 1; }
    lo=$(printf '%s\n' "$win" | awk -F= '/^LENGTH_LO=/{print $2; exit}')
    hi=$(printf '%s\n' "$win" | awk -F= '/^LENGTH_HI=/{print $2; exit}')
    peak=$(printf '%s\n' "$win" | awk -F= '/^PEAK=/{print $2; exit}')
    [[ -n "$lo" && -n "$hi" ]] || { echo "❌ [ONT] could not parse length window"; return 1; }
    echo "  [ONT] Length window: ${lo}-${hi} bp (peak ~${peak}), quality >= Q${qual}"

    local any=false fq stem
    for fq in "${in_dir}/"*.fastq*; do
        [[ -f "$fq" ]] || continue
        stem=$(basename "$fq"); stem="${stem%.gz}"; stem="${stem%.fastq}"
        # pipefail (in a subshell so it stays local) makes a failing zcat on a
        # corrupt .gz abort the sample instead of being masked by chopper's exit 0.
        if ! ( set -o pipefail
               { if [[ "$fq" == *.gz ]]; then zcat "$fq"; else cat "$fq"; fi; } \
                 | chopper -q "$qual" --minlength "$lo" --maxlength "$hi" -t "$threads" \
                 > "${out_dir}/${stem}.fastq" 2>>"${dataset_path%/}/tmp/ont_chopper.log" ); then
            echo "❌ [ONT] chopper failed on ${stem}"; return 1
        fi
        if [[ -s "${out_dir}/${stem}.fastq" ]]; then
            any=true
        else
            echo "  [ONT] ${stem}: 0 reads after filtering, dropping" >&2
            rm -f "${out_dir}/${stem}.fastq"
        fi
    done
    [[ "$any" == true ]] || { echo "❌ [ONT] no reads passed chopper filtering"; return 1; }
    touch "${out_dir}/.chopper_done"
    export ONT_FASTQ_DIR="$out_dir"
}

# ── Per-sample UNOISE3 denoise; concatenate centroids across samples ────────
# Reads are relabeled uniquely per sample (>stem.N) during fasta conversion to
# avoid cross-sample read-name collisions in the merged file (racon aborts on
# duplicate sequence names).
# In : $ONT_FASTQ_DIR  Out: tmp/step_06_ont/persample/<stem>_cluster.fasta
#                           tmp/step_06_ont/combined_centroids.fasta
Amplicon_ONT_ClusterPerSample() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path" || return 1
    trimmed_path="${dataset_path%/}"; dataset_name="${trimmed_path##*/}"
    local in_dir="${ONT_FASTQ_DIR:-${dataset_path%/}/tmp/step_03_chopper}"
    local work="${dataset_path%/}/tmp/step_06_ont"
    local ps="${work}/persample"
    local threads="${THREADS_PER_DATASET:-${cpu:-4}}"
    mkdir -p "$ps"
    : > "${work}/combined_centroids.fasta"

    local fq stem fasta cent
    for fq in "${in_dir}/"*.fastq; do
        [[ -f "$fq" ]] || continue
        stem=$(basename "$fq" .fastq)
        fasta="${ps}/${stem}.fasta"
        cent="${ps}/${stem}_cluster.fasta"
        awk -v s="$stem" 'NR%4==1{c++; print ">" s "." c} NR%4==2{print}' "$fq" > "$fasta"
        [[ -s "$fasta" ]] || { echo "  [ONT] ${stem}: empty after fasta conversion, skipping" >&2; continue; }
        # --strand both: ONT amplicon reads occur in both orientations; without
        # it forward and reverse-complement copies of one sequence would form
        # separate centroids.
        vsearch --cluster_unoise "$fasta" \
            --minsize 1 --sizeout \
            --strand both \
            --threads "$threads" \
            --centroids "$cent" --quiet 2>>"${work}/ont_vsearch.log" \
            || { echo "❌ [ONT] cluster_unoise failed on ${stem}"; return 1; }
        cat "$cent" >> "${work}/combined_centroids.fasta"
    done
    [[ -s "${work}/combined_centroids.fasta" ]] || { echo "❌ [ONT] no centroids produced"; return 1; }
}

# ── Per-sample minimap2 (map-ont) + racon consensus polishing ───────────────
# ref = that sample's centroids; reads = all samples' centroids (combined),
# exactly as ONT-AmpSeq rules 04-05. Targets without read support are dropped
# by racon (removes spurious variants); a sample whose racon output is empty
# falls back to its raw centroids.
# In : tmp/step_06_ont/persample/*_cluster.fasta + combined_centroids.fasta
# Out: tmp/step_06_ont/polish/<stem>_polished.fasta
Amplicon_ONT_PolishRacon() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path" || return 1
    local work="${dataset_path%/}/tmp/step_06_ont"
    local ps="${work}/persample"
    local mapd="${work}/map" pol="${work}/polish"
    local combined="${work}/combined_centroids.fasta"
    local threads="${THREADS_PER_DATASET:-${cpu:-4}}"
    mkdir -p "$mapd" "$pol"

    local cent stem sam
    for cent in "${ps}/"*_cluster.fasta; do
        [[ -f "$cent" ]] || continue
        stem=$(basename "$cent" _cluster.fasta)
        sam="${mapd}/${stem}.sam"
        # -K / -f match the ONT-AmpSeq reference (config: K=500M, f=0.0002):
        #   -K minibatch size; -f filter out the most frequent minimizers.
        minimap2 -ax map-ont \
            -K "${ONT_MINIMAP2_K:-500M}" -f "${ONT_MINIMAP2_F:-0.0002}" \
            --secondary=no -t "$threads" \
            "$cent" "$combined" > "$sam" 2>>"${work}/ont_minimap2.log" \
            || { echo "❌ [ONT] minimap2 failed on ${stem}"; return 1; }
        if racon -t "$threads" "$combined" "$sam" "$cent" \
                > "${pol}/${stem}_polished.fasta" 2>>"${work}/ont_racon.log" \
                && [[ -s "${pol}/${stem}_polished.fasta" ]]; then
            :
        else
            echo "  [ONT] ${stem}: racon empty/failed, using raw centroids" >&2
            cp "$cent" "${pol}/${stem}_polished.fasta"
        fi
    done
    [[ -n "$(ls -A "$pol" 2>/dev/null)" ]] || { echo "❌ [ONT] no polished output"; return 1; }
}

# ── Per-sample relabel (tag sample id, sort by size) + merge ────────────────
# In : tmp/step_06_ont/polish/*_polished.fasta
# Out: tmp/step_06_ont/merged_polished_relabeled.fasta
Amplicon_ONT_RelabelMerge() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path" || return 1
    local work="${dataset_path%/}/tmp/step_06_ont"
    local pol="${work}/polish" rel="${work}/relabel"
    local threads="${THREADS_PER_DATASET:-${cpu:-4}}"
    mkdir -p "$rel"
    : > "${work}/merged_polished_relabeled.fasta"

    local pf stem
    for pf in "${pol}/"*_polished.fasta; do
        [[ -f "$pf" ]] || continue
        stem=$(basename "$pf" _polished.fasta)
        vsearch --sortbysize "$pf" \
            --sample "$stem" \
            --threads "$threads" \
            --output "${rel}/${stem}_relabeled.fasta" --quiet 2>>"${work}/ont_vsearch.log" \
            || { echo "❌ [ONT] relabel/sortbysize failed on ${stem}"; return 1; }
        cat "${rel}/${stem}_relabeled.fasta" >> "${work}/merged_polished_relabeled.fasta"
    done
    [[ -s "${work}/merged_polished_relabeled.fasta" ]] \
        || { echo "❌ [ONT] merged relabeled file empty"; return 1; }
}

# ── Cluster merged polished seqs at fixed identity -> OTU representative seqs ─
# In : tmp/step_06_ont/merged_polished_relabeled.fasta
# Out: tmp/step_06_ont/otus.fasta  (relabeled OTU_1..N, with ;size=)
Amplicon_ONT_ClusterID() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path" || return 1
    local work="${dataset_path%/}/tmp/step_06_ont"
    local threads="${THREADS_PER_DATASET:-${cpu:-4}}"
    local id="${ONT_OTU_IDENTITY:-0.97}"
    # --strand both: collapse forward and reverse-complement representatives of
    # the same amplicon into one OTU (ONT reads are not strand-normalized).
    vsearch --cluster_fast "${work}/merged_polished_relabeled.fasta" \
        --id "$id" \
        --strand both \
        --threads "$threads" \
        --relabel OTU_ --sizein --sizeout \
        --centroids "${work}/otus.fasta" --quiet 2>>"${work}/ont_vsearch.log" \
        || { echo "❌ [ONT] cluster_fast failed"; return 1; }
    [[ -s "${work}/otus.fasta" ]] || { echo "❌ [ONT] no OTU centroids produced"; return 1; }
    echo "  [ONT] OTUs defined: $(grep -c '^>' "${work}/otus.fasta")"
}

# ── Build TRUE read-abundance table: map all reads back to OTUs (97%) + import ─
# Reuses relabel_reads_for_mapping (sample-labeled read fasta from manifest)
# and import_vsearch_to_qiime2, identical to the DegradedQ pipeline.
# In : tmp/step_06_ont/otus.fasta + manifest (built by caller)
# Out: tmp/step_07_cluster/<ds>-table-clustered.qza, -repseq-clustered.qza
Amplicon_ONT_MapReadsToOTUs() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path" || return 1
    trimmed_path="${dataset_path%/}"; dataset_name="${trimmed_path##*/}"
    local work="${dataset_path%/}/tmp/step_06_ont"
    local cluster_path="${dataset_path%/}/tmp/step_07_cluster"
    local manifest="${dataset_path%/}/tmp/temp_file/${dataset_name}_manifest.tsv"
    local threads="${THREADS_PER_DATASET:-${cpu:-4}}"
    mkdir -p "$cluster_path"

    [[ -f "$manifest" ]] || { echo "❌ [ONT] manifest not found: $manifest"; return 1; }

    echo ">>> [ONT] Relabeling reads for mapping..."
    python3 "${SCRIPTS}/py_16s.py" relabel_reads_for_mapping \
        --manifest_path "$manifest" \
        --output_fasta "${work}/all_reads_labeled.fasta" \
        --threads "$threads" || { echo "❌ [ONT] relabel_reads_for_mapping failed"; return 1; }

    # OTUs are clustered at 97% (accurate, racon-polished), but raw ONT reads
    # carry ~5-10% native error, so mapping them back at 97% would discard most
    # of them. Map at a lower, ONT-error-tolerant identity (default 0.90, which
    # recovers ~97% of reads in testing). Tunable via ONT_MAP_IDENTITY.
    #   --strand both        : ONT reads occur in both orientations.
    #   --maxaccepts/--maxrejects : bound the per-read search for the best OTU.
    #     Fully exhaustive (0/0) assigns every read to its GLOBALLY best OTU but
    #     is O(reads x OTUs) — at realistic ONT depth (e.g. ~19k reads x ~1.4k
    #     OTUs) it runs for tens of minutes per dataset. Defaults below give a
    #     near-best assignment orders of magnitude faster; set both to 0 via
    #     ONT_MAP_MAXACCEPTS/ONT_MAP_MAXREJECTS for the exhaustive behaviour.
    local map_id="${ONT_MAP_IDENTITY:-0.90}"
    local map_acc="${ONT_MAP_MAXACCEPTS:-8}"
    local map_rej="${ONT_MAP_MAXREJECTS:-64}"
    echo ">>> [ONT] Mapping all reads to OTUs (id=${map_id}, maxaccepts=${map_acc}/maxrejects=${map_rej}) for abundance table..."
    vsearch --usearch_global "${work}/all_reads_labeled.fasta" \
        --db "${work}/otus.fasta" \
        --id "$map_id" \
        --strand both \
        --maxaccepts "$map_acc" --maxrejects "$map_rej" \
        --otutabout "${work}/otu_table.tsv" \
        --threads "$threads" --quiet 2>>"${work}/ont_vsearch.log" \
        || { echo "❌ [ONT] usearch_global mapping failed"; return 1; }

    echo ">>> [ONT] Importing OTU table + rep-seqs into QIIME2..."
    python3 "${SCRIPTS}/py_16s.py" import_vsearch_to_qiime2 \
        --zotu_fasta "${work}/otus.fasta" \
        --otu_table_tsv "${work}/otu_table.tsv" \
        --manifest_path "$manifest" \
        --output_table_qza "${cluster_path}/${dataset_name}-table-clustered.qza" \
        --output_repseq_qza "${cluster_path}/${dataset_name}-repseq-clustered.qza" \
        || { echo "❌ [ONT] import_vsearch_to_qiime2 failed"; return 1; }
}

################################################################################
#                              END OF FILE                                     #
################################################################################

