#!/usr/bin/env bash

Amplicon_454_ImportToQiime2() {
    
    cd $dataset_path
    local fastq_path=$dataset_path"ori_fastq/"
    local fastp_path=$dataset_path"temp/step_02_fastp/"
    local temp_file_path=$dataset_path"temp/temp_file/"
    mkdir -p "$temp_file_path" "$fastp_path"

    echo "It's 454 data, directly import to qiime2"
    
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    find $fastq_path -type f > $temp_file_path$dataset_name"-file.txt"
    py_16s.py mk_manifest_SE --FilePath $temp_file_path$dataset_name"-file.txt"
    
    local temp_path=$dataset_path"temp/"
    local temp_file_path=$dataset_path"temp/temp_file/"
    local qza_path=$dataset_path"temp/step_03_qza_import/"
    mkdir -p "$temp_file_path" "$qza_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    
    
    
    qiime tools import \
        --type 'SampleData[SequencesWithQuality]' \
        --input-path $temp_file_path$dataset_name"_manifest.tsv" \
        --output-path $qza_path$dataset_name".qza" \
        --input-format SingleEndFastqManifestPhred33V2


}