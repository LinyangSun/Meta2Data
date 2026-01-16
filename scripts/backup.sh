

Amplicon_Illumina_DenosingDeblur() {
    # Parse flags: -s <start>, -e <end>
    local start_in="" end_in="" opt
    while getopts ":s:e:" opt; do
        case "$opt" in
            s) start_in="$OPTARG" ;;
            e) end_in="$OPTARG" ;;
            \?) echo "Invalid option: -$OPTARG" >&2; return 2 ;;
            :)  echo "Option -$OPTARG requires an argument." >&2; return 2 ;;
        esac
    done
    shift $((OPTIND-1))
    # Safety
    set -euo pipefail
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    # Paths
    local base="${dataset_path%/}"
    local quality_filter_path="${base%/}/temp/step_04_qza_import_QualityFilter/"
    local deblur_path="${base%/}/temp/step_05_denoise/deblur/"
    local qf_vis_path="${base%/}/temp/temp_file/QualityFilter_vis/"
    local qf_view_path="${base%/}/temp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${base%/}/temp/temp_file/QualityFilter_vis/qf_trim_pos/"
    mkdir -p "$deblur_path" "$qf_view_path" "$qf_trim_pos_path" "$qf_vis_path"
    local dataset_name="${base##*/}"
    # Decide if we need to compute positions
    local need_compute=true
    [[ -n "$start_in" && -n "$end_in" ]] && need_compute=false
    echo  "start_in: $start_in ; end_in $end_in"
    local start_calc="" end_calc=""
    echo "need_compute is : $need_compute"
    if $need_compute; then
        echo "Computing trim positions from QIIME 2 visualization..."
        qiime demux summarize \
            --i-data "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
            --o-visualization "${qf_vis_path%/}/${dataset_name}_import_cutadapt_QualityFilter.qzv"
        unzip -q -o "${qf_vis_path%/}/${dataset_name}_import_cutadapt_QualityFilter.qzv" -d "$qf_view_path"
        # Fixed wrong var name from original; keep forward file
        find "$qf_view_path" -type f -name 'forward-seven-number-summaries.tsv' -exec cp -f {} "${qf_trim_pos_path%/}/" \;
        rm -rf "$qf_view_path"
        local tsv_path="${qf_trim_pos_path%/}/forward-seven-number-summaries.tsv"
        if [[ ! -s "$tsv_path" ]]; then
            echo "ERROR: Expected TSV not found: $tsv_path" >&2
            return 1
        fi
        local trim_pos_result
        trim_pos_result="$(py_16s.py trim_pos_deblur --FilePath "$tsv_path")"
        echo "$trim_pos_result"
        IFS=',' read -r final_start final_end <<< "$trim_pos_result"
        echo "Computed: start=$final_start, end=$final_end"
    fi
    # Final positions: use flags if set, otherwise computed
    local start="${start_in:-$final_start}"
    local end="${end_in:-$final_end}"
    echo "Using trim positions: start=$start, end=$end"
    echo "$start $end" > "${qf_trim_pos_path%/}/Trim_position.txt"
    # Deblur
    qiime deblur denoise-16S \
        --i-demultiplexed-seqs "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --p-trim-length "$end" \
        --p-left-trim-len "$start" \
        --o-representative-sequences "${deblur_path%/}/${dataset_name}-rep-seqs-deblur.qza" \
        --o-table "${deblur_path%/}/${dataset_name}-table-deblur.qza" \
        --p-sample-stats \
        --o-stats "${deblur_path%/}/${dataset_name}-deblur-stats.qza"
    echo "Deblur complete."
}

Amplicon_Illumina_PrimerDetectionAndQualityControl() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local fastq_path="${dataset_path%/}/ori_fastq/"
    local fastp_path="${dataset_path%/}/temp/step_02_fastp/"
    local temp_file_path="${dataset_path%/}/temp/temp_file/"
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



###############################################################################
#                         454 PLATFORM FUNCTIONS                               #
################################################################################
# Functions for Roche 454 pyrosequencing data processing


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
Amplicon_454_Deduplication() {
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