#!/bin/bash
#
# AmpliconFunction.sh - Unified Amplicon Processing Functions
#
# This file consolidates all amplicon sequencing processing functions
# for different platforms (454, Illumina, PacBio, ONT) and common utilities.
#
# Functions are organized by platform and processing stage:
#   - Common Functions: Shared utilities across all platforms
#   - 454 Functions: Roche 454 pyrosequencing pipeline
#   - Illumina Functions: Illumina short-read pipeline
#   - PacBio Functions: PacBio HiFi/CCS long-read pipeline
#   - ONT Functions: Oxford Nanopore long-read pipeline (TODO)
#
# Usage: Source this file in your pipeline script
#   source scripts/AmpliconFunction.sh
#
# Author: Meta2Data Project
# Date: $(date +%Y-%m-%d)
#

set -e


################################################################################
#                         COMMON FUNCTIONS                                     #
################################################################################
# Functions shared across all sequencing platforms


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
Amplicon_Common_ImportFastqToQiime2() {
    set -u
    cd "$dataset_path" || { echo "❌ Cannot access dataset path: $dataset_path"; exit 1; }

    local temp_path="${dataset_path%/}/temp/"
    local temp_file_path="${temp_path}temp_file/"
    local qza_path="${temp_path}step_03_qza_import/"
    mkdir -p "$temp_file_path" "$qza_path"

    local dataset_name="${dataset_path%/}"
    dataset_name="${dataset_name##*/}"

    local paired_manifest="${temp_file_path}${dataset_name}_manifest.tsv"
    local forward_manifest="${temp_file_path}${dataset_name}_manifest_forward.tsv"
    local merge_log="${temp_path}${dataset_name}_vsearch_merge.log"

    echo "🔹 Processing dataset: $dataset_name"

    # --- SINGLE-END ---
    if [ "${sequence_type:-paired}" = "single" ]; then
        echo "🧬 Importing single-end reads..."
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path "$paired_manifest" \
            --output-path "${qza_path}${dataset_name}.qza" \
            --input-format SingleEndFastqManifestPhred33V2
        echo "✅ Single-end import completed."
        return
    fi

    # --- PAIRED-END ---
    echo "🧬 Importing paired-end reads..."
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$paired_manifest" \
        --output-path "${qza_path}${dataset_name}.qza" \
        --input-format PairedEndFastqManifestPhred33V2

    echo "🔄 Merging paired-end reads..."
    # Try merging, but don't exit on failure; capture logs.
    if ! qiime vsearch merge-pairs \
        --i-demultiplexed-seqs "${qza_path}${dataset_name}.qza" \
        --o-unmerged-sequences "${qza_path}${dataset_name}_unmerged.qza" \
        --o-merged-sequences "${qza_path}${dataset_name}_join.qza" \
        &> "$merge_log"
    then
        echo "⚠️ vsearch merge failed. See log: $merge_log"
        echo "↩️ Falling back to forward reads only..."

        # Build forward-only manifest from paired manifest
        awk -F'\t' '
            BEGIN { OFS="\t" }
            NR==1 { 
                print "sample-id","absolute-filepath"
                next 
            }
            NR>1 {
                print $1, $2
            }
        ' "$paired_manifest" > "$forward_manifest"

        # Replace any existing paired artifact with forward-only import
        rm -f "${qza_path}${dataset_name}.qza"
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path "$forward_manifest" \
            --output-path "${qza_path}${dataset_name}.qza" \
            --input-format SingleEndFastqManifestPhred33V2

        echo "✅ Forward-only import completed after merge failure."
        return
    fi

    # --- CHECK MERGE QUALITY (only if merge succeeded) ---
    local unmerged_qza="${qza_path}${dataset_name}_unmerged.qza"
    local merged_qza="${qza_path}${dataset_name}_join.qza"
    local unmerged_size merged_size

    unmerged_size=$(stat -c%s "$unmerged_qza")
    merged_size=$(stat -c%s "$merged_qza")

    echo "📏 unmerged.qza = $unmerged_size bytes"
    echo "📏 merged.qza   = $merged_size bytes"

    if (( $(echo "$unmerged_size * 1.5 > $merged_size" | bc -l) )); then
        echo "⚠️ Less than ~70% merged — switching to forward reads only..."

        awk -F'\t' '
            BEGIN { OFS="\t" }
            NR==1 { 
                print "sample-id","absolute-filepath"
                next 
            }
            NR>1 {
                print $1, $2
            }
        ' "$paired_manifest" > "$forward_manifest"

        rm -f "${qza_path}${dataset_name}.qza"
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path "$forward_manifest" \
            --output-path "${qza_path}${dataset_name}.qza" \
            --input-format SingleEndFastqManifestPhred33V2

        echo "✅ Forward-only import completed."
    else
        echo "✅ Merge quality acceptable. Using merged reads."
        rm -f "${qza_path}${dataset_name}.qza"
        mv "$merged_qza" "${qza_path}${dataset_name}.qza"
    fi
}
Amplicon_Common_FinalFilesCleaning() {
    dataset_path="${dataset_path%/}/"
    local quality_filter_path="${dataset_path%/}/temp/step_04_qza_import_QualityFilter/"
    local deblur_path="${dataset_path%/}/temp/step_05_denoise/deblur/"
    local qf_vis_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/"
    local qf_view_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/qf_trim_pos/"
    local qc_vis="${dataset_path%/}/temp/temp_file/qc_vis/"
    local denosing_vis="${dataset_path%/}/temp/temp_file/denosing_vis/"
    
    cd "$dataset_path" || { echo "Error: dataset_path not found"; return 1; }
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    quality_filter_path="${dataset_path%/}/temp/step_04_qza_import_QualityFilter/"
    deblur_path="${dataset_path%/}/temp/step_05_denoise/deblur/"
    qf_trim_pos_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/qf_trim_pos/"
    qc_vis="${dataset_path%/}/temp/temp_file/qc_vis/"
    denosing_vis="${dataset_path%/}/temp/temp_file/denosing_vis/"
    mkdir -p "$qc_vis" "$denosing_vis"
    
    # Check for Deblur output
    if [ -d "$deblur_path" ] && [ -f "${deblur_path%/}/${dataset_name}-table-deblur.qza" ]; then
        echo "The analysis ran successfully! Now only keeping final files to save space."
        cp "${deblur_path%/}/${dataset_name}-rep-seqs-deblur.qza" "${dataset_path%/}/${dataset_name}-final-rep-seqs.qza"
        cp "${deblur_path%/}/${dataset_name}-table-deblur.qza" "${dataset_path%/}/${dataset_name}-final-table.qza"
        
        # Copy trim position files if they exist
        if [ -f "${qf_trim_pos_path%/}/Trim_position.txt" ]; then
            cp "${qf_trim_pos_path%/}/Trim_position.txt" "${dataset_path%/}/${dataset_name}-TrimPosition.txt"
        fi
        if [ -f "${qf_trim_pos_path%/}/forward-seven-number-summaries.tsv" ]; then
            cp "${qf_trim_pos_path%/}/forward-seven-number-summaries.tsv" "${dataset_path%/}/${dataset_name}-TrimPositionOriginalFile.tsv"
        fi
        
        # Extract stats
        if [ -f "${quality_filter_path%/}/${dataset_name}_filter-stats.qza" ]; then
            unzip -q "${quality_filter_path%/}/${dataset_name}_filter-stats.qza" -d "$qc_vis"
            find "$qc_vis" -type f -name 'stats.csv' -exec cp {} "${dataset_path%/}/${dataset_name}-QCStats.csv" \;
        fi
        
        if [ -f "${deblur_path%/}/${dataset_name}-deblur-stats.qza" ]; then
            unzip -q "${deblur_path%/}/${dataset_name}-deblur-stats.qza" -d "$denosing_vis"
            find "$denosing_vis" -type f -name 'stats.csv' -exec cp {} "${dataset_path%/}/${dataset_name}-DenosingStats.csv" \;
        fi
        
        # Remove temporary directories but keep dataset_path itself
        find "$dataset_path" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
        
        # Remove unwanted logs
        rm -f "${dataset_path%/}/"{deblur.log,fastp.html,fastp.json,wild-v4-1_sra.txt}
        
        return 0
        
    # Check for vsearch output (454 pipeline)
    elif [ -f "${dataset_path%/}/${dataset_name}-table-vsearch.qza" ]; then
        echo "The analysis ran successfully! Now only keeping final files to save space."
        mv "${dataset_path%/}/${dataset_name}-table-vsearch.qza" "${dataset_path%/}/${dataset_name}-final-table.qza"
        mv "${dataset_path%/}/${dataset_name}-rep-seqs-vsearch.qza" "${dataset_path%/}/${dataset_name}-final-rep-seqs.qza"
        find "$dataset_path" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
        
        return 0
        
    else
        echo "❌ ERROR: The analysis failed! The final denoising output does not exist."
        echo "   Expected files:"
        echo "     - ${deblur_path%/}/${dataset_name}-table-deblur.qza (Illumina)"
        echo "     - ${dataset_path%/}/${dataset_name}-table-vsearch.qza (454)"
        echo ""
        echo "   Check previous steps for errors."
        
        # Don't clean up temp files so we can debug
        return 1
    fi
}


################################################################################
#                       ILLUMINA PLATFORM FUNCTIONS                            #
################################################################################
# Functions for Illumina short-read sequencing data processing



Amplicon_Illumina_QualityControlForQZA(){
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local qza_path="${dataset_path%/}/temp/step_03_qza_import/"
    local quality_filter_path="${dataset_path%/}/temp/step_04_qza_import_QualityFilter/"
    local qf_vis_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/"
    mkdir -p "$quality_filter_path" "$qf_vis_path"
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    qiime quality-filter q-score \
        --i-demux "${qza_path%/}/${dataset_name}.qza" \
        --p-min-quality 20 \
        --o-filtered-sequences "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --o-filter-stats "${quality_filter_path%/}/${dataset_name}_filter-stats.qza" \
        --verbose
}

Amplicon_Illumina_DenosingDada2() {
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
    qiime dada2 denoise-pyro \
        --i-demultiplexed-seqs "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --p-trunc-len "$end" \
        --p-trim-left "$start" \
        --o-representative-sequences "${deblur_path%/}/${dataset_name}-rep-seqs-denoising.qza" \
        --o-table "${deblur_path%/}/${dataset_name}-table-denoising.qza" \
        --o-denoising-stats "${deblur_path%/}/${dataset_name}-denoising-stats.qza"
    echo "denoising complete."
}
################################################################################
#                        PACBIO PLATFORM FUNCTIONS                             #
################################################################################
# Functions for PacBio HiFi/CCS long-read sequencing data processing


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
Amplicon_Pacbio_QualityControlForQZA() {
    dataset_path="${dataset_path%/}/"
    cd "$dataset_path"
    local qza_path="${dataset_path%/}/temp/step_03_qza_import/"
    local quality_filter_path="${dataset_path%/}/temp/step_04_qza_import_QualityFilter/"
    local qf_vis_path="${dataset_path%/}/temp/temp_file/QualityFilter_vis/"
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
Amplicon_Pacbio_DenosingDada2() {
    
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
    # Deblur
    qiime dada2 denoise-ccs \
        --i-demultiplexed-seqs "${quality_filter_path%/}/${dataset_name}_QualityFilter.qza" \
        --p-max-len 1600 \
        --o-table "${deblur_path%/}/${dataset_name}-table-deblur.qza" \
        --o-representative-sequences "${deblur_path%/}/${dataset_name}-rep-seqs-deblur.qza" \
        --o-denoising-stats "${deblur_path%/}/${dataset_name}-deblur-stats.qza"
}

################################################################################
#                         ONT PLATFORM FUNCTIONS                               #
################################################################################
# Functions for Oxford Nanopore long-read sequencing data processing
# TODO: These functions are placeholders and need to be implemented

# Amplicon_ONT_PrimerDetectionAndQualityControl() {
#     echo "TODO: Implement ONT primer detection and quality control"
#     echo "  - Adapter trimming for ONT reads"
#     echo "  - Quality filtering for ONT specific error profiles"
#     echo "  - Primer detection and trimming"
#     return 1
# }

# Amplicon_ONT_QualityControlForQZA() {
#     echo "TODO: Implement ONT quality control for QZA"
#     echo "  - Import ONT reads to QIIME2"
#     echo "  - Quality filtering based on ONT characteristics"
#     return 1
# }

# Amplicon_ONT_DenosingMethod() {
#     echo "TODO: Implement ONT denoising method"
#     echo "  - Evaluate denoising approaches for ONT (DADA2? Custom?)"
#     echo "  - Handle ONT-specific error patterns"
#     echo "  - Generate ASV/OTU table"
#     return 1
# }


################################################################################
#                              END OF FILE                                     #
################################################################################

