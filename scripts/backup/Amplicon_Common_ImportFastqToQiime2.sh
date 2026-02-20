#!/usr/bin/env bash
Amplicon_Common_ImportFastqToQiime2() {
    set -u
    cd "$dataset_path" || { echo "❌ Cannot access dataset path: $dataset_path"; exit 1; }

    local temp_path="${dataset_path%/}/tmp/"
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