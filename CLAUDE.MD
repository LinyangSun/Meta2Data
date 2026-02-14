# Meta2Data - Bioinformatics Pipeline

Bioinformatics tool for automated amplicon sequencing data processing from NCBI databases. Downloads, processes, and analyzes metagenomic data via CLI.

## ⚠️ IMPORTANT: Setup Environment First

**Run these commands immediately after reading this file:**

```bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate Meta2Data

```

This ensures you have the correct environment and iSeq alias configured for all Meta2Data operations.

never do:
cd /test && source ~/anaconda3/etc/profile.d/conda.sh && conda activate Meta2Data 2>&1 && iseq -i 

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

### Conda (Recommended)

**Environment Details:**
- **Name**: `Meta2Data`
- **Python**: 3.10
- **Base**: QIIME2 2024.10 (Amplicon distribution)

**Included Software:**
- QIIME2 2024.10 with all plugins
- vsearch (deduplication, clustering, chimera removal)
- fastp (quality control, adapter trimming)
- seqkit (FASTQ/FASTA processing)
- sra-tools (SRA download and conversion)
- q2-greengenes2 (taxonomy database)
- pandas, numpy, biopython, biom-format (data processing)
- iseq (SRA/ENA/DDBJ/GSA download - see below for setup)

**Activate:**
```bash

conda activate Meta2Data
```



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
