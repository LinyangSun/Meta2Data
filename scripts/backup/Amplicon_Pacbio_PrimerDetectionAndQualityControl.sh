#!/usr/bin/env bash

Amplicon_Pacbio_PrimerDetectionAndQualityControl() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local fastq_path="${dataset_path%/}/ori_fastq/"
    local fastp_path="${dataset_path%/}/temp/step_02_fastp/"
    local temp_file_path="${dataset_path%/}/temp/temp_file/"
    mkdir -p "$temp_file_path" "$fastp_path"
    echo "The data is Pacbio, no need to remove primer before analysis, just copy original file to fastp folder"
    cp "${fastq_path%/}"/*.fastq "${fastp_path%/}/"

}