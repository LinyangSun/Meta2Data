#!/usr/bin/env bash
Amplicon_Pacbio_DenosingDada2() {
    
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
    # Deblur
    qiime dada2 denoise-ccs \
        --i-demultiplexed-seqs "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --p-max-len 1600 \
        --p-front ${F_primer} \
        --p-adapter ${R_primer} \
        --o-table "${deblur_path%/}/${dataset_name}-table-deblur.qza" \
        --o-representative-sequences "${deblur_path%/}/${dataset_name}-rep-seqs-deblur.qza" \
        --o-denoising-stats "${deblur_path%/}/${dataset_name}-deblur-stats.qza"
}