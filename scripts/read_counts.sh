# Read-count hooks are inactive when helpers are used outside the PIP runner.
Audit_Counts() {
    [[ -n "${READ_COUNTS_STATE:-}" ]] || return 0
    python3 "${SCRIPTS}/read_counts.py" "$@" --state "$READ_COUNTS_STATE"
}
Audit_Fastq() {
    Audit_Counts fastq --stage "$1" --input "$2" --parent "${3:-}" "${@:4}"
}
Audit_Table() {
    Audit_Counts table --stage "$1" --input "$2" --parent "${3:-}"
}
Audit_Fasta() {
    Audit_Counts fasta --stage "$1" --input "$2" --parent "${3:-}"
}
Audit_Dada2() {
    Audit_Counts stats --stage dada2 --input "${dataset_path%/}/tmp/step_05_denoise/${dataset_name}-denoising-stats.qza"
}
Audit_QC() {
    Audit_Counts stats --stage qc --input "${dataset_path%/}/tmp/step_04_qza_import_QualityFilter/${dataset_name}_filter-stats.qza"
}
Audit_Exit() {
    local rc="$1" status=failed
    case "$rc" in 0) status=success ;; 98|99) status=skipped ;; esac
    Audit_Counts status --input "$status"
}
