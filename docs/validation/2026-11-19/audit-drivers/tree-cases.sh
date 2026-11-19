#!/usr/bin/env bash
set -euo pipefail
M2D_RUN=/workspaces/meta2data-validation
trap 'm2d_status=$?; printf "%s\n" "$m2d_status" > "$M2D_RUN/tree-audit-01.exit"' EXIT
export TMPDIR="$M2D_RUN/tmp" PATH="/workspaces/Meta2Data/bin:$PATH" PYTHONNOUSERSITE=1 MPLBACKEND=agg UNIFRAC_USE_GPU=N
for tree_case in singleV sepp; do
  case_dir="$M2D_RUN/tree-audit-01/$tree_case"
  mkdir -p "$case_dir/results"
  ln -s "$M2D_RUN/results/db" "$case_dir/results/db"
  tree_args=()
  if [[ "$tree_case" == singleV ]]; then tree_args+=(--singleV); fi
  cd "$case_dir"
  /opt/conda/bin/conda run --no-capture-output --prefix "$M2D_RUN/env-release-frozen" Meta2Data AmpliconTAXA -i "$M2D_RUN/workflow-audit-01/pip-dada2/results/pip/pe_explicit" --dada2 --classifier silva-144 -t 2 "${tree_args[@]}" > "$M2D_RUN/logs/tree-audit-01-$tree_case.log" 2>&1
  printf '%s PASS\n' "$tree_case"
done
