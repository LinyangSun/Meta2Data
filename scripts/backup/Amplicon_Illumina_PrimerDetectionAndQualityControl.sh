#!/usr/bin/env bash

Amplicon_Illumina_PrimerDetectionAndQualityControl() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local fastq_path="${dataset_path%/}/ori_fastq/"
    local fastp_path="${dataset_path%/}/tmp/step_02_fastp/"
    local temp_file_path="${dataset_path%/}/tmp/temp_file/"
    mkdir -p "$temp_file_path" "$fastp_path"

    echo "start fastq2fastp"

    if [ "$sequence_type" = "single" ]; then
        first_file=$(find "$fastq_path" -type f -name "*.fastq" | head -n 1)
        head -10000 "$first_file" > "${temp_file_path%/}/f10000_R1.fastq"
        seqkit grep -m 4 -p ${F_primer} "${temp_file_path%/}/f10000_R1.fastq" > "${temp_file_path%/}/primer_detection.txt"
        file_size=$(stat -c %s "${temp_file_path%/}/primer_detection.txt")
        if [ $file_size -gt 0 ]; then
            echo "###################"
            echo "Primer exist... removeing"
            echo "###################"
            f_len=$(echo -n ${F_primer} | wc -c)
            find "$fastq_path" -type f -name "*.fastq"| while read -r file; do
                filename=$(basename "$file")
                read=$file
                out="${fastp_path%/}/$filename"
                fastp -Q -G -i $read -o $out --trim_front1 $f_len
            done
        else
            echo "###################"
            echo "Primer has been removed..."
            echo "###################"
            find "$fastq_path" -type f -name "*.fastq"| while read -r file; do
                filename=$(basename "$file")
                read=$file
                out="${fastp_path%/}/$filename"
                fastp -Q -G -i $read -o $out
            done
        fi
    else
        first_file=$(find "$fastq_path" -type f -name "*_1.fastq" | head -n 1)
        head -10000 "$first_file" > "${temp_file_path%/}/f10000_R1.fastq"
        seqkit grep -m 4 -p ${F_primer} "${temp_file_path%/}/f10000_R1.fastq" > "${temp_file_path%/}/primer_detection.txt"
        file_size=$(stat -c %s "${temp_file_path%/}/primer_detection.txt")
        if [ $file_size -gt 0 ]; then
            echo "###################"
            echo "Primer exist... removeing"
            echo "###################"
            f_len=$(echo -n ${F_primer} | wc -c)
            r_len=$(echo -n ${R_primer} | wc -c)
            find "$fastq_path" -type f -name "*_1.fastq"| while read -r file; do
                filename=$(basename "$file")
                base_name=$(echo $filename | sed 's/_1\.fastq//')
                read1="${fastq_path%/}/${base_name}_1.fastq"
                read2="${fastq_path%/}/${base_name}_2.fastq"
                out1="${fastp_path%/}/${base_name}_1.fastq"
                out2="${fastp_path%/}/${base_name}_2.fastq"
                fastp -Q -G -i $read1 -o $out1 -I $read2 -O $out2 --trim_front1 $f_len --trim_front2 $r_len
            done
        else
            echo "###################"
            echo "Primer has been removed..."
            echo "###################"
            find "$fastq_path" -type f -name "*_1.fastq"| while read -r file; do
                filename=$(basename "$file")
                base_name=$(echo $filename | sed 's/_1\.fastq//')
                read1="${fastq_path%/}/${base_name}_1.fastq"
                read2="${fastq_path%/}/${base_name}_2.fastq"
                out1="${fastp_path%/}/${base_name}_1.fastq"
                out2="${fastp_path%/}/${base_name}_2.fastq"
                fastp -Q -G -i $read1 -o $out1 -I $read2 -O $out2
            done
        fi
    fi
}