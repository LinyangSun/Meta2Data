#!/usr/bin/env bash
Amplicon_Common_FinalFilesCleaning() {
    dataset_path="${dataset_path%/}/"
    local quality_filter_path="${dataset_path%/}/tmp/step_04_qza_import_QualityFilter/"
    local deblur_path="${dataset_path%/}/tmp/step_05_denoise/deblur/"
    local qf_vis_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/"
    local qf_view_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/qf_view/"
    local qf_trim_pos_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/qf_trim_pos/"
    local qc_vis="${dataset_path%/}/tmp/temp_file/qc_vis/"
    local denosing_vis="${dataset_path%/}/tmp/temp_file/denosing_vis/"
    
    cd "$dataset_path" || { echo "Error: dataset_path not found"; return 1; }
    trimmed_path="${dataset_path%/}"
    dataset_name="${trimmed_path##*/}"
    quality_filter_path="${dataset_path%/}/tmp/step_04_qza_import_QualityFilter/"
    deblur_path="${dataset_path%/}/tmp/step_05_denoise/deblur/"
    qf_trim_pos_path="${dataset_path%/}/tmp/temp_file/QualityFilter_vis/qf_trim_pos/"
    qc_vis="${dataset_path%/}/tmp/temp_file/qc_vis/"
    denosing_vis="${dataset_path%/}/tmp/temp_file/denosing_vis/"
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