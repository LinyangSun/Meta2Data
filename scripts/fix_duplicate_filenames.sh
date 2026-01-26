#!/bin/bash
#
# fix_duplicate_filenames.sh - Remove duplicate Run IDs from filenames
#
# Usage: fix_duplicate_filenames.sh <directory>
#
# Fixes filenames like:
#   PRJCA040882_CRR1878501_CRR1878501_1.fastq.gz
#   → PRJCA040882_CRR1878501_1.fastq.gz
#
#   PRJCA040882-CRR1878501-CRR1878501_1.fastq.gz (legacy)
#   → PRJCA040882_CRR1878501_1.fastq.gz
#
#   PRJNA123_SRR456_SRR456_1.fastq
#   → PRJNA123_SRR456_1.fastq

set -e

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <directory>"
    echo ""
    echo "Example: $0 /path/to/dataset/ori_fastq/"
    exit 1
fi

target_dir="${1%/}/"

if [[ ! -d "$target_dir" ]]; then
    echo "Error: Directory not found: $target_dir"
    exit 1
fi

echo "Scanning directory: $target_dir"
echo "----------------------------------------"

renamed_count=0
skipped_count=0

# Find all fastq files
find "$target_dir" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | while read -r filepath; do
    filename=$(basename "$filepath")
    dirname=$(dirname "$filepath")

    # Pattern: PROJECT_RUNID_RUNID_1.fastq(.gz) or PROJECT-RUNID-RUNID_1.fastq(.gz)
    # Extract the run ID that appears twice
    if [[ "$filename" =~ ^([^_-]+)[_-]([CDE]RR[0-9]+)[_-]\2(_[12]\.(fastq|fastq\.gz))$ ]]; then
        project="${BASH_REMATCH[1]}"
        runid="${BASH_REMATCH[2]}"
        suffix="${BASH_REMATCH[3]}"

        new_filename="${project}_${runid}${suffix}"

        echo "  Renaming: $filename"
        echo "        → $new_filename"

        mv "$filepath" "${dirname}/${new_filename}"
        ((renamed_count++))

    elif [[ "$filename" =~ ^([^_-]+)[_-]([A-Z]RR[0-9]+)[_-]\2(\.(fastq|fastq\.gz))$ ]]; then
        # Single-end pattern without _1/_2
        project="${BASH_REMATCH[1]}"
        runid="${BASH_REMATCH[2]}"
        suffix="${BASH_REMATCH[3]}"

        new_filename="${project}_${runid}${suffix}"

        echo "  Renaming: $filename"
        echo "        → $new_filename"

        mv "$filepath" "${dirname}/${new_filename}"
        ((renamed_count++))
    else
        ((skipped_count++))
    fi
done

echo "----------------------------------------"
echo "Summary:"
echo "  Renamed: $renamed_count files"
echo "  Skipped: $skipped_count files (no duplication detected)"
echo "Done!"
