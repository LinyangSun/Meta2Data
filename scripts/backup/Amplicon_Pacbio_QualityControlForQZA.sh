#!/usr/bin/env bash
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