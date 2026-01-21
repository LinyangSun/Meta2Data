# Meta2Data

**Automated Bioinformatics Pipeline for Amplicon Sequencing Data Processing from NCBI Databases**

Meta2Data is a comprehensive command-line tool for downloading, processing, and analyzing metagenomic amplicon sequencing data from public databases (NCBI, CNCB/GSA). It integrates metadata retrieval, SRA data download, quality control, and QIIME2-based analysis into a single, automated workflow.

## Features

- **Metadata Download**: Search and download metadata from NCBI and CNCB/GSA databases with parallel processing and checkpoint/resume capability
- **Multi-Platform Support**: Automatic detection and processing of Illumina, PacBio, Ion Torrent, and 454 sequencing platforms
- **Smart Primer Detection**: Automatic primer detection and trimming for amplicon data
- **QIIME2 Integration**: Seamless integration with QIIME2 2024.10 for downstream analysis
- **Parallel Processing**: Concurrent processing with configurable thread counts
- **Robust Error Handling**: Checkpoint system, retry mechanism, and comprehensive logging
- **Taxonomy Assignment**: Optional GreenGenes2 taxonomy classification and phylogenetic tree generation

## Installation

### Option 1: Conda Environment (Recommended)

The conda environment includes all dependencies (QIIME2, vsearch, fastp, sra-tools, seqkit, etc.):

```bash
# Create environment from env.yml
conda env create -f env.yml

# Activate environment
conda activate Meta2Data

# Verify installation
Meta2Data --help
```

### Option 2: Pip Installation

Install Meta2Data package into an existing conda environment with QIIME2:

```bash
# Install from GitHub
pip install git+https://github.com/LinyangSun/Meta2Data.git@main

# Or install from local repository
pip install -e .

# Verify installation
Meta2Data --help
```

**Note**: Pip installation only installs Meta2Data scripts and Python dependencies. You must have QIIME2 and other bioinformatics tools (vsearch, fastp, sra-tools, seqkit) already installed in your environment.

## Requirements

### System Requirements
- **OS**: Linux (tested on Ubuntu/CentOS)
- **Memory**: Minimum 8GB RAM, 16GB+ recommended for large datasets
- **Storage**: Varies by dataset size (SRA files can be large)

### Software Dependencies
- Python ≥ 3.10
- QIIME2 2024.10 (Amplicon distribution)
- vsearch 2.22+
- fastp
- sra-tools (including Aspera for fast downloads)
- seqkit
- q2-greengenes2 (for taxonomy assignment)

### Python Packages
- biopython
- pandas
- numpy
- openpyxl
- requests

All dependencies are included in the conda environment (`env.yml`).

## Usage

Meta2Data provides several subcommands for different stages of the workflow:

### Overview

```bash
Meta2Data <command> [options]

Available commands:
    MetaDL         Enhanced metadata download with parallel processing (NCBI + CNCB)
    AmpliconPIP    Download and process amplicon sequencing data
    Evaluate       Summarize processing results
    ShortreadsPIP  (In development)
```

### 1. MetaDL: Enhanced Metadata Download

Download metadata from NCBI and CNCB databases with parallel processing and checkpoint/resume capability.

**Two modes available:**

#### Mode 1: BioProject ID Input

```bash
# Basic usage
Meta2Data MetaDL \
    -i bioproject_ids/ \
    -o metadata_output/

# With NCBI API key (faster, 8 parallel workers)
Meta2Data MetaDL \
    -i bioproject_ids/ \
    -o metadata_output/ \
    -k YOUR_NCBI_API_KEY
```

**Input format**: Directory containing `.txt` files with one BioProject ID per line
- Supports: PRJNA*, PRJEB*, PRJDB*, PRJCA*, CRA*, etc.

#### Mode 2: Keyword Search

```bash
# Search by keywords
Meta2Data MetaDL \
    -o metadata_output/ \
    --keywords \
    --field "16S rRNA" "amplicon" \
    --organism "gut microbiome" "bacteria"

# With optional terms
Meta2Data MetaDL \
    -o metadata_output/ \
    --keywords \
    --field "metagenome" \
    --organism "soil" \
    --opt "Illumina"
```

**Output files:**
- `all_metadata_merged.csv` - Final merged dataset (NCBI + CNCB)
- `ncbi_merged_*.csv` - NCBI data
- `cncb_combined.csv` - CNCB data
- `logs/metadl_v2_*.log` - Detailed execution log
- `checkpoints/download_state.json` - Resume state for interrupted runs

### 2. AmpliconPIP: Amplicon Data Processing

Download SRA data and process amplicon sequencing data through platform-specific pipelines.

#### Basic Usage

```bash
# Basic processing (with standard column names)
Meta2Data AmpliconPIP \
    -m metadata.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
    -t 8

# Specify output directory
Meta2Data AmpliconPIP \
    -m metadata.csv \
    -o /output/path/ \
    --col-bioproject Bioproject \
    --col-sra Run \
    -t 8

# Enable GreenGenes2 taxonomy assignment (requires backbone and taxonomy files)
Meta2Data AmpliconPIP \
    -m metadata.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
    -t 8 \
    --gg2 \
    --i-backbone /path/to/backbone.qza \
    --i-reference-taxonomy /path/to/taxonomy.qza

# Test mode with sample data (no column names needed)
Meta2Data AmpliconPIP --test -t 8
```

#### Custom Column Names

If your CSV uses different column names:

```bash
Meta2Data AmpliconPIP \
    -m metadata.csv \
    --col-bioproject "ProjectID" \
    --col-sra "SRA_Accession" \
    -t 8
```

#### Metadata CSV Format

**Required columns** (customizable with `--col-*` options):
- `Bioproject` - BioProject accession (e.g., PRJNA123456)
- `Run` - SRA run accession (e.g., SRR123456)

**Example CSV:**
```csv
Bioproject,Run
PRJNA12345,SRR123456
PRJNA12345,SRR123457
PRJNA67890,SRR234567
```

#### Processing Pipeline

The pipeline automatically:
1. **Downloads** SRA data (via Aspera/FTP)
2. **Detects** sequencing platform (Illumina, PacBio, Ion Torrent, 454)
3. **Identifies** and trims primers
4. **Processes** data using platform-specific methods:
   - **Illumina/Ion Torrent**: DADA2 denoising
   - **PacBio**: DADA2 with PacBio parameters
   - **454**: Vsearch clustering and chimera removal
5. **Generates** QIIME2 artifacts (`.qza` files)
6. **Assigns** taxonomy (optional, with `--gg2`)

#### Output Structure

```
<output_dir>/
├── datasets_ID.txt                    # Generated dataset list
├── <dataset_ID>/                      # One directory per dataset
│   ├── <dataset_ID>_sra.txt          # SRA accession list
│   ├── ori_fastq/                     # Downloaded FASTQ files
│   ├── <dataset_ID>-final-rep-seqs.qza       # Final sequences
│   └── <dataset_ID>-final-table.qza          # Final feature table
├── final/                             # (if --gg2 enabled)
│   └── merged/                        # Merged results with taxonomy
├── failed_datasets.log
├── success_datasets.log
└── skipped_datasets.log
```

### 3. Evaluate: Summarize Results

(Feature in development - check `Meta2Data Evaluate --help` for details)

## Command-Line Options

### MetaDL Options

```
Required:
    -o, --output DIR              Output directory

Mode 1: BioProject Input
    -i, --input DIR               Directory with BioProject ID txt files

Mode 2: Keyword Search
    --keywords                    Enable keyword search mode
    --field "term1" "term2"       Search field terms
    --organism "term1" ...        Organism terms
    --opt "term1" ...             Optional additional terms

Optional:
    -e, --email EMAIL             Email (auto-generated if not provided)
    -k, --api-key KEY             NCBI API key (enables 8 workers)
    -w, --max-workers NUM         Max parallel workers (default: 3 or 8)
    -h, --help                    Show help
```

### AmpliconPIP Options

```
Required (unless --test is used):
    -m, --metadata FILE           Input metadata CSV file
    --col-bioproject NAME         Column name for BioProject in CSV
    --col-sra NAME                Column name for SRA accession in CSV

Optional:
    -o, --output DIR              Output directory (default: metadata file dir)
    -t, --threads INT             CPU threads (default: 4)
    --test                        Use test data (test/ampliconpiptest.csv)
    -h, --help                    Show help

GreenGenes2 Taxonomy Options:
    --gg2                         Enable GreenGenes2 taxonomy assignment
    --i-backbone FILE             Backbone tree file (required if --gg2 enabled)
    --i-reference-taxonomy FILE   Reference taxonomy file (required if --gg2 enabled)
```

## Examples

### Complete Workflow

```bash
# 1. Activate environment
conda activate Meta2Data

# 2. Download metadata
Meta2Data MetaDL \
    -i bioproject_list/ \
    -o metadata/ \
    -k YOUR_API_KEY

# 3. Process amplicon data with taxonomy
Meta2Data AmpliconPIP \
    -m metadata/all_metadata_merged.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
    -o results/ \
    -t 16 \
    --gg2 \
    --i-backbone /path/to/backbone.qza \
    --i-reference-taxonomy /path/to/taxonomy.qza

# 4. Check results
ls results/final/merged/
```

### Process Specific Platform Data

```bash
# Process only Illumina data (filter metadata first)
grep "ILLUMINA" metadata.csv > illumina_only.csv
Meta2Data AmpliconPIP \
    -m illumina_only.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
    -t 8
```

### Resume Interrupted Download

If MetaDL is interrupted, it automatically resumes from the last checkpoint:

```bash
# Simply rerun the same command
Meta2Data MetaDL -i bioproject_ids/ -o metadata_output/
# Output: "Resuming from checkpoint..."
```

## Troubleshooting

### Command Not Found After Pip Install

If `Meta2Data` command is not found after `pip install`:

```bash
# Reinstall ensuring scripts are installed
pip uninstall Meta2Data
pip install --force-reinstall git+https://github.com/LinyangSun/Meta2Data.git@main

# Verify scripts are in PATH
which Meta2Data
```

### SRA Download Issues

If SRA downloads fail:
1. Check internet connectivity
2. Verify sra-tools installation: `which prefetch`
3. Configure Aspera for faster downloads: `vdb-config --interactive`
4. Check available disk space

### Memory Issues

For large datasets, increase memory allocation:
- Reduce parallel workers: `-t 4` instead of `-t 16`
- Process datasets individually
- Use high-memory compute nodes on HPC systems

## Development

### Project Structure

```
Meta2Data/
├── bin/                    # Command-line entry points
│   ├── Meta2Data           # Main dispatcher
│   ├── Meta2Data-MetaDL    # Metadata download v2
│   └── Meta2Data-AmpliconPIP # Amplicon pipeline
├── scripts/                # Processing scripts
│   ├── AmpliconFunction.sh # Amplicon processing functions
│   ├── metadata_downloader.py # Python metadata downloader
│   └── py_16s.py          # 16S data processing
├── test/                   # Test data
├── docs/                   # Documentation
├── env.yml                 # Conda environment
├── pyproject.toml          # Python package config
└── setup.py               # Script installation config
```

### Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Follow existing code style (bash best practices)
4. Test with sample data
5. Submit a pull request

## Citation

If you use Meta2Data in your research, please cite:

```
[Citation information to be added]
```

## License

[License information to be added]

## Contact

- **Issues**: https://github.com/LinyangSun/Meta2Data/issues
- **Repository**: https://github.com/LinyangSun/Meta2Data

## Acknowledgments

This project integrates several excellent bioinformatics tools:
- QIIME2 for microbiome analysis
- DADA2 for sequence denoising
- Vsearch for sequence clustering
- GreenGenes2 for taxonomy classification
- SRA Toolkit for data access
