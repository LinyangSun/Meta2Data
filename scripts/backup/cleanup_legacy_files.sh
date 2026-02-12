#!/bin/bash
#
# cleanup_legacy_files.sh - Remove legacy dash-separated files
#
# Removes old files with dash separator, keeping only underscore versions

set -e

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <directory>"
    echo ""
    echo "Example: $0 /scratch/project_2009135/linyang/test/PRJCA040882/ori_fastq/"
    exit 1
fi

target_dir="${1%/}/"

if [[ ! -d "$target_dir" ]]; then
    echo "Error: Directory not found: $target_dir"
    exit 1
fi

echo "Scanning directory: $target_dir"
echo "Looking for legacy dash-separated files..."
echo "----------------------------------------"

removed_count=0

# Find files with dash separator pattern: PROJECT-RUNID-X.fastq
find "$target_dir" -type f \( -name "*-[CDE]RR[0-9]*-[12].fastq*" -o -name "*-[A-Z]RR[0-9]*-[12].fastq*" \) | while read -r filepath; do
    filename=$(basename "$filepath")

    # Check if underscore version exists
    underscore_version=$(echo "$filename" | sed 's/-\([CDE]RR[0-9]\+\)-/_\1_/')

    if [[ -f "${target_dir}${underscore_version}" ]]; then
        echo "  Removing legacy: $filename"
        echo "  (Underscore version exists: $underscore_version)"
        rm -f "$filepath"
        ((removed_count++))
    else
        echo "  Skipping: $filename (no underscore version found)"
    fi
done

echo "----------------------------------------"
echo "Removed: $removed_count legacy files"
echo "Done!"
