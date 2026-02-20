#!/usr/bin/env bash
Amplicon_Illumina_DenosingDeblur() {
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
    local quality_filter_path="${base%/}/tmp/step_04_qza_import_QualityFilter/"
    local deblur_path="${base%/}/tmp/step_05_denoise/deblur/"
    local qf_vis_path="${base%/}/tmp/temp_file/QualityFilter_vis/"
    local qf_view_path="${base%/}/tmp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${base%/}/tmp/temp_file/QualityFilter_vis/qf_trim_pos/"
    mkdir -p "$deblur_path" "$qf_view_path" "$qf_trim_pos_path" "$qf_vis_path"
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
    qiime deblur denoise-16S \
        --i-demultiplexed-seqs "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --p-trim-length "$end" \
        --p-left-trim-len "$start" \
        --o-representative-sequences "${deblur_path%/}/${dataset_name}-rep-seqs-deblur.qza" \
        --o-table "${deblur_path%/}/${dataset_name}-table-deblur.qza" \
        --p-sample-stats \
        --o-stats "${deblur_path%/}/${dataset_name}-deblur-stats.qza"
    echo "Deblur complete."
}