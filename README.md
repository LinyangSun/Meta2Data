# Meta2Data

**Automated Bioinformatics Pipeline for Amplicon Sequencing Data Processing from NCBI Databases**

Meta2Data is a comprehensive command-line tool for downloading, processing, and analyzing metagenomic amplicon sequencing data from public databases (NCBI, CNCB/GSA). It integrates metadata retrieval, SRA data download, quality control, and QIIME2-based analysis into a single, automated workflow.

## Features

- **Metadata Download and pre-clean**: Search, download and pre-clean for metadata from NCBI and CNCB databases by keywords or provied bioproject id
- **Multi-Platform Support**: Automatic detection and processing of Illumina, PacBio, Ion Torrent, and 454 sequencing platforms(ONT not supported)
- **Smart Primer Detection**: Automatic primer detection and trimming for amplicon data
- **QIIME2 Integration**: Seamless integration with QIIME2 2024.10 for downstream analysis
- **Taxonomy Assignment**: GreenGenes2 taxonomy classification and phylogenetic tree generation

## Notice

If you are running the AmpliconPIP or ggCOMBO functions on a server, please ensure you allocate sufficient time for your task. The data download and taxonomy annotation steps can be very time-consuming.

## Installation

### Option 1: Conda Environment (Recommended)

copy the **env.yml** to your device.

env.yml : qiime2 included (QIIME2, vsearch, fastp, sra-tools, seqkit, q2-greengenes, Meta2Data, etc.)

```bash
# Create environment from env1.yml
conda env create -f env.yml

# Verify installation
Meta2Data --help
```

### Option 2: Pip Installation

Install Meta2Data package only in your current conda env:

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
- **Storage**: Varies by dataset size (SRA files could be large)

### Software Dependencies
- Python ≥ 3.10
- QIIME2 2024.10 (Amplicon distribution)
- vsearch 2.22+
- fastp
- sra-tools
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
    AmpliconPIP    Download and process 16s amplicon sequencing data (ITS AND 18S not included)
    ggCOMBO        Merge AmpliconPIP results and assign taxonomy via GreenGenes2
    ShortreadsPIP  (In development)
```

### 1. MetaDL: Enhanced Metadata Download

Download metadata from NCBI and CNCB databases with parallel processing and checkpoint/resume capability. Also included the basic column combination and format standardization. 

**Two modes available:**

#### Mode 1: BioProject ID Input

In this mode, the user might already have a list of bioproject id. the bioproject id could be separate as multiple files, but store under same dir. In this dir you should only include the bioproject id files. For each bioproject id files, it should only have one column and tab separated. please see the example under example/MetaDL/bioprojectID.txt

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
- Supports: PRJNA*, PRJEB*, PRJDB*, PRJCA*, etc.

#### Mode 2: Keyword Search

```bash
# Search by keywords
Meta2Data MetaDL \
    -o metadata_output/ \
    --keywords \
    --field "16S rRNA" \
    --organism "gut microbiome"

# With optional terms
Meta2Data MetaDL \
    -o metadata_output/ \
    --keywords \
    --field "metagenome" \
    --organism "soil" \
    --opt "Illumina"

# or use a set of keywords
field=("Bacteria" "Microbiome" "Microbes" "Metagenomics" "Metabarcoding")
organism=("bee" "apis" "bombus")
opt=("Amplicon" "16s" "skin" "gut")

Meta2Data MetaDL \
    -o metadata_output/ \
    --keywords \
    --field "${field[@]}" \
    --organism "${organism[@]}"

# With optional terms
Meta2Data MetaDL \
    -o metadata_output/ \
    --keywords \
    --field "${field[@]}" \
    --organism "${organism[@]}" \
    --opt "${opt[@]}"
```

**Output files:**
- `all_metadata_merged.csv` - Final merged dataset (NCBI + CNCB)
- `ncbi_merged_*.csv` - NCBI, EBI & DDBJ data
- `cncb_combined.csv` - CNCB data
- `logs/metadl_v2_*.log` - Detailed execution log
- `checkpoints/download_state.json` - Resume state for interrupted runs

### 2. AmpliconPIP: Amplicon Data Processing

Download SRA data and process amplicon sequencing data with provided metadata.

#### Basic Usage

```bash
# Basic processing (with standard column names)
Meta2Data AmpliconPIP \
    -m path/to/metadata.csv \
    --col-bioproject "Bioproject" \
    --col-sra "Run" \
    -t 8

# Or specify output directory
Meta2Data AmpliconPIP \
    -m path/to/metadata.csv \
    -o /output/path/ \
    --col-bioproject "Bioproject" \
    --col-sra "Run" \
    -t 8

# Test mode with built-in sample data (output to current directory)
Meta2Data AmpliconPIP --test -t 8

# Test mode with built-in data and custom output directory
Meta2Data AmpliconPIP --test -o /path/to/output/ -t 8

# Test mode with user metadata (subsets to 2 SRA per BioProject)
Meta2Data AmpliconPIP --test \
    -m metadata.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
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
               sequencing type (single / pair-end)
3. **Identifies** and trims primers
4. **Processes** data using platform-specific methods:
   - **Illumina/Ion Torrent**: DADA2 denoising
   - **PacBio**: DADA2 with PacBio parameters
5. **Generates** QIIME2 artifacts (`.qza` files)

#### Output Structure

```
<output_dir>/
├── datasets_ID.txt                    # Generated dataset list
├── <dataset_ID>/                      # One directory per dataset
│   ├── <dataset_ID>_sra.txt          # SRA accession list
│   ├── ori_fastq/                     # Downloaded FASTQ files
│   ├── <dataset_ID>-final-rep-seqs.qza       # Final sequences
│   └── <dataset_ID>-final-table.qza          # Final feature table
├── failed_datasets.log
├── success_datasets.log
├── skipped_datasets.log
└── summary.csv                       # sequencing depth for raw reads / final data
```

### 3. ggCOMBO: Merge & Taxonomy Assignment

Merge multiple AmpliconPIP dataset outputs and assign taxonomy using the GreenGenes2 database. Can be run independently from AmpliconPIP.

#### Basic Usage

```bash
# Input and output in the same directory
Meta2Data ggCOMBO \
    --db /path/to/gg2db \
    -i /path/to/amplicon_output \
    -t 8

# Separate output directory
Meta2Data ggCOMBO \
    --db /path/to/gg2db \
    -i /path/to/amplicon_output \
    -o /path/to/results \
    -t 8

# if not download database, use this
Meta2Data ggCOMBO \
    --db /path/to/gg2db \
    --dl \
    -i /path/to/amplicon_output \
    -t 8
```

#### Input

The `--input` directory must contain one or more `PRJ*` subdirectories (output of AmpliconPIP), each with:
- `<dataset_ID>-final-table.qza` — feature table
- `<dataset_ID>-final-rep-seqs.qza` — representative sequences

The `--db` directory must contain these files:
- `2024.09.backbone.full-length.nb.sklearn-1.4.2.qza` (GG2 pre-trained Naive Bayes classifier)
- `sepp-refs-gg-13-8.qza` (SEPP reference tree)

Use `--dl` to download these automatically.

#### Processing Pipeline

1. **Merge** feature tables and representative sequences across all datasets
2. **Assign** taxonomy via pre-trained Naive Bayes classifier (sklearn)
3. **Build** phylogenetic tree via SEPP fragment insertion
4. **Filter** features by tree placement

#### Output Structure

```
<output_dir>/
└── final/
    └── merged/
        ├── merged-table.qza               # Merged feature table (complete)
        ├── merged-rep-seqs.qza            # Merged representative sequences
        ├── merged-table-summary.qzv       # Table summary
        ├── merged-taxonomy.qza            # Taxonomy (sklearn NB classifier)
        ├── insertion-tree.qza             # SEPP phylogenetic tree
        ├── insertion-placements.qza       # SEPP placement details
        ├── merged-table-tree.qza          # Table filtered to tree-placed features
        └── merged-table-no-tree.qza       # Features not placed in tree
```


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
    -o, --output DIR              Output directory
                                  Default: current directory in --test mode
                                           metadata file directory in normal mode
    -t, --threads INT             CPU threads (default: 4)
    --test                        Run in test mode
                                  Without -m: use built-in test data
                                  With -m: subset metadata (2 SRA per BioProject)
                                  Output defaults to current directory
    -h, --help                    Show help
```

### ggCOMBO Options

```
Required:
    --db DIR                      Path to GreenGenes2 database directory
    -i, --input DIR               Input directory containing PRJ* dataset folders
                                  (output of a previous AmpliconPIP run)

Optional:
    -o, --output DIR              Output directory for merged results
                                  (default: same as --input)
    -t, --threads INT             CPU threads (default: 4)
    --dl                          Download GreenGenes2 database files to --db directory
    -h, --help                    Show help
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

# 3. Process amplicon data
Meta2Data AmpliconPIP \
    -m metadata/all_metadata_merged.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
    -o results/ \
    -t 16

# 4. Merge and assign taxonomy with GreenGenes2
Meta2Data ggCOMBO \
    --db /path/to/gg2db \
    --dl \
    -i results/ \
    -t 16

# 5. Check results
ls results/final/merged/
```



## Development

### Project Structure

```
Meta2Data/
├── bin/                    # Command-line entry points
│   ├── Meta2Data           # Main dispatcher
│   ├── Meta2Data-MetaDL    # Metadata download
│   ├── Meta2Data-AmpliconPIP # Amplicon pipeline
│   └── Meta2Data-ggCOMBO   # Merge & taxonomy assignment
├── scripts/                # Processing scripts
│   ├── AmpliconFunction.sh # Amplicon processing functions
│   ├── run.sh             # AmpliconPIP orchestrator
│   ├── taxonomy.sh        # ggCOMBO taxonomy pipeline
│   ├── metadata_downloader.py # Python metadata downloader
│   └── py_16s.py          # 16S data processing
├── test/                   # Test data
├── docs/                   # Documentation & reference sequences
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

