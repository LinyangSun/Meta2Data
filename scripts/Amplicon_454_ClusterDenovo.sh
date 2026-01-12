#!/usr/bin/env bash
Amplicon_454_ClusterDenovo() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local fastp_path="${dataset_path%/}/ori_fastq/"
    local temp_file_path="${dataset_path%/}/temp/temp_file/"
    local qza_path="${dataset_path%/}/temp/step_03_qza_import/"
    local qc_path="${dataset_path%/}/temp/step_04_quality_control/"
    local qf_view_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/qf_trim_pos/"
    local dedupicate_path="${dataset_path%/}/temp/step_05_dedupicate/"
    local chemeras_path="${dataset_path%/}/temp/step_06_ChimerasRemovel/"
    mkdir -p "$temp_file_path" "$fastp_path"
    
    qiime vsearch cluster-features-de-novo \
          --i-table "${chemeras_path%/}/${dataset_name}_trim_single-table-nochimeratable.qza" \
          --i-sequences "${chemeras_path%/}/${dataset_name}_trim_single-nochimerarepseq.qza" \
          --p-perc-identity 0.97 \
          --o-clustered-table "${dataset_path%/}/${dataset_name}-table-vsearch.qza" \
          --o-clustered-sequences "${dataset_path%/}/${dataset_name}-rep-seqs-vsearch.qza"
    
}