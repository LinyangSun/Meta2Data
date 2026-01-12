#!/usr/bin/env bash

Amplicon_Common_MakeManifestFileForQiime2() {
    cd $dataset_path
    local temp_file_path=$dataset_path"temp/temp_file/"
    local fastp_path=$dataset_path"temp/step_02_fastp/"
    mkdir -p "$temp_file_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    find $fastp_path -type f > $temp_file_path$dataset_name"-file.txt"
    if [ "$sequence_type" = "single" ]; then
        py_16s.py mk_manifest_SE --FilePath $temp_file_path$dataset_name"-file.txt"
    else
        py_16s.py mk_manifest_PE --FilePath $temp_file_path$dataset_name"-file.txt"
    fi
}