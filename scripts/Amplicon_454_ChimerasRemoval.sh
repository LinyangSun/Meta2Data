#!/usr/bin/env bash
Amplicon_454_ChimerasRemoval() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    # Extract dataset name
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    local fastp_path="${dataset_path%/}/ori_fastq/"
    local temp_file_path="${dataset_path%/}/temp/temp_file/"
    local qza_path="${dataset_path%/}/temp/step_03_qza_import/"
    local qc_path="${dataset_path%/}/temp/step_04_quality_control/"
    local qf_view_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/qf_trim_pos/"
    local dedupicate_path="${dataset_path%/}/temp/step_05_dedupicate/"
    local chemeras_path="${dataset_path%/}/temp/step_06_ChimerasRemovel/"
    mkdir -p "$chemeras_path"
    
    qiime vsearch uchime-denovo \
          --i-sequences "${dedupicate_path%/}/${dataset_name}_trim_single-repseq-nosingle.qza" \
          --i-table "${dedupicate_path%/}/${dataset_name}_trim_single-table-nosingletable.qza" \
          --o-chimeras "${chemeras_path%/}/${dataset_name}_trim_single-chimeras.qza" \
          --o-nonchimeras "${chemeras_path%/}/${dataset_name}_trim_single-nochimeras.qza" \
          --o-stats "${chemeras_path%/}/${dataset_name}_trim_single-uchimstats.qza"
# Keep only non-chimeric features in the table
    qiime feature-table filter-features \
          --i-table "${dedupicate_path%/}/${dataset_name}_trim_single-table-nosingletable.qza" \
          --m-metadata-file "${chemeras_path%/}/${dataset_name}_trim_single-nochimeras.qza" \
          --o-filtered-table "${chemeras_path%/}/${dataset_name}_trim_single-table-nochimeratable.qza"
# Filter sequences to non-chimeras
    qiime feature-table filter-seqs \
         --i-data "${dedupicate_path%/}/${dataset_name}_trim_single-repseq-nosingle.qza" \
         --m-metadata-file "${chemeras_path%/}/${dataset_name}_trim_single-nochimeras.qza" \
         --o-filtered-data "${chemeras_path%/}/${dataset_name}_trim_single-nochimerarepseq.qza"
}