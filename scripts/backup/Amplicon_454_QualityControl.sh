#!/usr/bin/env bash
Amplicon_454_QualityControl() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local fastp_path="${dataset_path%/}/ori_fastq/"
    local temp_file_path="${dataset_path%/}/temp/temp_file/"
    local qza_path="${dataset_path%/}/temp/step_03_qza_import/"
    local qc_path="${dataset_path%/}/temp/step_04_quality_control/"
    local qf_view_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/qf_trim_pos/"
    local trimmed_path dataset_name trim_pos_result trim_length f_len r_len shorter overlap
    mkdir -p "$temp_file_path" "$qc_path" "$qf_trim_pos_path" "$qf_view_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    # Quality visualization
    qiime demux summarize \
        --i-data "${qza_path%/}/${dataset_name}.qza" \
        --o-visualization "${qf_view_path%/}/${dataset_name}_import_cutadapt_QualityFilter.qzv" \
        || return 1
    # Extract trim positions
    unzip -q -o "${qf_view_path%/}/${dataset_name}_import_cutadapt_QualityFilter.qzv" -d "${qf_view_path%/}/unzipped/"
    find "${qf_view_path%/}/unzipped/" -type f -name 'forward-seven-number-summaries.tsv' -exec cp {} "$qf_trim_pos_path" \;
    trim_pos_result=$(py_16s.py trim_pos_454 --FilePath "${qf_trim_pos_path%/}/forward-seven-number-summaries.tsv")
    trim_length=$(echo "$trim_pos_result" | cut -d',' -f1)
    echo "Trim length: $trim_length" | tee "${qf_trim_pos_path%/}/Trim_position.txt"
    f_len=$(echo -n "$F_primer" | wc -c)
    r_len=$(echo -n "$R_primer" | wc -c)
    shorter=$(( f_len < r_len ? f_len : r_len ))
    overlap=$(echo "$shorter * 7 / 10" | bc)
    [ "$overlap" -lt 8 ] && overlap=8

    qiime cutadapt trim-single \
            --i-demultiplexed-sequences "${qza_path%/}/${dataset_name}.qza" \
            --p-front "$F_primer" \
            --p-adapter "$REVRC" \
            --p-error-rate 0.10 \
            --p-minimum-length "$trim_length" \
            --o-trimmed-sequences "${qc_path%/}/${dataset_name}_trim_single.qza"
}