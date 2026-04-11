# Meta2Data - Bioinformatics Pipeline

Bioinformatics tool for automated amplicon sequencing data processing from NCBI databases. Downloads, processes, and analyzes metagenomic data via CLI.

## ⚠️ IMPORTANT: Setup Environment First

**Run these commands immediately after reading this file:**

```bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate Meta2Data

```

This ensures you have the correct environment configured for all Meta2Data operations.

## Project Structure

```
bin/                    # CLI entry points (Meta2Data, MetaDL, MetaProcessingRun.sh)
scripts/               # Platform-specific pipelines (Amplicon_454/Illumina/Pacbio_*.sh)
test/                  # Test files and data
```

## Main Subcommands

- **MetaDL** - NCBI metadata search and download
- **AmpliconPIP** - Download and process amplicon sequencing (main pipeline)
- **Evaluate** - Summarize processing results
- **ShortreadsPIP** - (TBD)

## Quick Start

```bash
conda activate Meta2Data
./bin/Meta2Data AmpliconPIP -m metadata.csv -o output_dir -t 8
```

## AmpliconPIP Parameters

**Required:**
- `-m, --metadata` - Input CSV file

**Optional:**
- `-o, --output` - Output directory (default: metadata file directory)
- `-t, --threads` - CPU threads (default: 4)
- `--skip-download` - Skip SRA download step
- `--keep-tmp` - Keep temporary files
- `--col-*` - Custom CSV column names (e.g., `--col-bioproject`, `--col-platform`, `--col-sra`, `--col-biosample`, `--col-forward`, `--col-reverse`, `--col-region`)

**Example:**
```bash
Meta2Data AmpliconPIP -m data.csv -t 8 --col-bioproject "ProjectID" --col-platform "Platform"
```

## Processing Pipeline

Workflow: SRA download → Quality control → Platform-specific processing → Qiime2 import → Cleanup

**Platform-specific steps:**
- **454**: Chimera removal, cluster denovo, post-clustering low-freq OTU filter
- **Illumina**: Denoising with Deblur
- **PacBio**: Denoising with DADA2

## Environment Setup

**Platform:** Linux only (x86_64 or aarch64).

Meta2Data is a self-contained tool that lives entirely inside the cloned
repository. There is **no `pip install` step** — Meta2Data is never copied
into any conda env. Users expose it to their shell by adding `<repo>/bin`
to `PATH`.

QIIME2 is a **runtime dependency of `AmpliconPIP` / `ggCOMBO` only**, not
of Meta2Data itself. `MetaDL` runs on any Python 3 interpreter.

**Install flow:**

```bash
# 1. Clone and put <repo>/bin on PATH.
git clone https://github.com/LinyangSun/Meta2Data.git
echo 'export PATH="'"$PWD"'/Meta2Data/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# 2. (AmpliconPIP / ggCOMBO only) Install and activate QIIME2 2024.10
#    amplicon following QIIME2's official instructions.
conda activate qiime2-amplicon-2024.10
```

That's it. `vsearch` + `fastp` + Python packages all auto-install on the
first subcommand run (see below).

**Native binaries (vsearch / fastp) are auto-installed at runtime.**
`Meta2Data-AmpliconPIP` and `Meta2Data-ggCOMBO` source
`scripts/check_dependencies.sh` and call `meta2data_ensure_vendor_binaries`
before doing any work. If `<repo>/vendor/bin/vsearch` or
`<repo>/vendor/bin/fastp` is missing, it invokes
`scripts/install_binaries.sh`, which pulls prebuilt Linux static binaries
from upstream releases (no compilation, no conda). The installer is
idempotent. `MetaDL` does **not** call this helper because it does not
need these binaries. Users who prefer to install eagerly can still run
`bash scripts/install_binaries.sh` manually — see its `--help` for
`--prefix` and `--force` overrides.

**Python packages are auto-installed at runtime.** Every Meta2Data entry
script (`Meta2Data-MetaDL`, `Meta2Data-AmpliconPIP`, `Meta2Data-ggCOMBO`)
sources `scripts/check_dependencies.sh` and calls
`meta2data_ensure_python_deps` before doing any work. That helper probes
the active `python3` for `biopython`, `pandas`, `numpy`, and `requests`,
and runs `python3 -m pip install` for anything missing — so packages
always land in whichever Python env the user currently has on PATH.

Both auto-install helpers honor `META2DATA_SKIP_DEP_CHECK=1` for users
who manage their environment by hand (HPC login nodes with read-only
prefixes, offline machines, etc.).

**Software used at runtime:**
- QIIME2 2024.10 amplicon distribution — **only needed for `AmpliconPIP` /
  `ggCOMBO`**; provides `qiime` plugins, pandas, numpy, requests,
  biom-format, cutadapt, and its own vsearch for plugin use
- `vsearch` 2.30.0 — installed by `install_binaries.sh`; used directly by the
  DegradedQ pipeline (UNOISE3, uchime3, usearch_global)
- `fastp` 0.24.0 — installed by `install_binaries.sh`; adapter removal and
  quality filtering for all Illumina / 454 / Ion Torrent / PacBio inputs
- `biopython`, `pandas`, `numpy`, `requests` — auto-installed into the
  active `python3` on first run. `biopython` is used by both
  `metadata_downloader.py` (NCBI Entrez for `MetaDL`) and `py_16s.py`
  (`SeqIO` parsing + NCBI platform detection inside `AmpliconPIP`)
- `wget`, `unzip`, `tar`, `gzip`, `md5sum` — expected to be present on any
  standard Linux distribution



## Development

**Shell Script Standards:**
- Shebang: `#!/bin/bash`
- Error handling: `set -e`
- Parameter validation upfront
- Help messages use heredocs

**Common Patterns:**
- Argument parsing: while loops with case statements
- Directory resolution: `$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)`
- Script calls: Absolute paths via `${SCRIPTS}/script_name.sh`

**Add New Subcommand:**
1. Create `bin/Meta2Data-SubcommandName`
2. Add case entry in `bin/Meta2Data` dispatcher
3. Create processing script in `scripts/` if needed

**Key Files:**
- `bin/Meta2Data` - CLI dispatcher
- `bin/MetaDL` - Metadata download handler
- `bin/MetaProcessingRun.sh` - Processing orchestration
- `scripts/Function_Import.sh` - Shared functions

## Notes

- Preserve scientific accuracy in bioinformatics operations
- Maintain shell script best practices (quoting, error handling)
- Follow existing naming conventions
- Scripts are modular and pipeline-oriented
- Output directories created automatically
- **IMPORTANT**: Debug and test only in `Meta2Data` conda for consistency
