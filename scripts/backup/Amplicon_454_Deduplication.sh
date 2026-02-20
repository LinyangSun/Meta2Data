#!/usr/bin/env bash
Amplicon_454_Deduplication() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    # Extract dataset name
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    local fastp_path="${dataset_path%/}/ori_fastq/"
    local temp_file_path="${dataset_path%/}/tmp/temp_file/"
    local qza_path="${dataset_path%/}/tmp/step_03_qza_import/"
    local qc_path="${dataset_path%/}/tmp/step_04_quality_control/"
    local qf_view_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/qf_trim_pos/"
    local dedupicate_path="${dataset_path%/}/tmp/step_05_dedupicate/"
    mkdir -p "$dedupicate_path"
    
    qiime vsearch dereplicate-sequences \
        --i-sequences "${qc_path%/}/${dataset_name}_trim_single.qza" \
        --o-dereplicated-table "${dedupicate_path%/}/${dataset_name}_trim_single-table.qza" \
        --o-dereplicated-sequences "${dedupicate_path%/}/${dataset_name}_trim_single-repseq.qza"
    
    qiime feature-table filter-features \
        --i-table "${dedupicate_path%/}/${dataset_name}_trim_single-table.qza" \
        --p-min-frequency 2 \
        --o-filtered-table "${dedupicate_path%/}/${dataset_name}_trim_single-table-nosingletable.qza"
    
    qiime feature-table filter-seqs \
        --i-data "${dedupicate_path%/}/${dataset_name}_trim_single-repseq.qza" \
        --i-table "${dedupicate_path%/}/${dataset_name}_trim_single-table-nosingletable.qza" \
        --o-filtered-data "${dedupicate_path%/}/${dataset_name}_trim_single-repseq-nosingle.qza"
    
}