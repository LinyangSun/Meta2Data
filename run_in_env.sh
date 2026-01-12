#!/bin/bash
# Wrapper script to run commands in Meta2Data conda environment
# Usage: ./run_in_env.sh <command>

set -e

# Initialize conda for bash
if [ -f ~/anaconda3/etc/profile.d/conda.sh ]; then
    source ~/anaconda3/etc/profile.d/conda.sh
elif [ -f ~/miniconda3/etc/profile.d/conda.sh ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
else
    echo "Error: Cannot find conda initialization script"
    exit 1
fi

# Activate Meta2Data environment
conda activate Meta2Data

# Load bash aliases (needed for iseq alias)
shopt -s expand_aliases
source ~/.bashrc 2>/dev/null || true
source ~/.zshrc 2>/dev/null || true

# Also manually define iseq alias as fallback
alias iseq='docker run --platform linux/amd64 --rm -v $(pwd):/workspace -w /workspace iseq iseq'

# Run the command passed as arguments
exec bash -c "$*"
