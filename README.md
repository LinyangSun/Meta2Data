# Meta2Data

**Automated Bioinformatics Pipeline for Amplicon Sequencing Data Processing from NCBI Databases**

Meta2Data is a comprehensive command-line tool for downloading, processing, and analyzing metagenomic amplicon sequencing data from public databases (NCBI, CNCB/GSA). It integrates metadata retrieval, SRA data download, quality control, and QIIME2-based analysis into a single, automated workflow.

## Features

- **Metadata Download and Pre-clean**: Search, download, and pre-clean metadata from NCBI and CNCB databases by keywords or provided BioProject ID
- **Multi-Platform Support**: Automatic detection and processing of Illumina, PacBio, Ion Torrent, and 454 sequencing platforms (ONT not supported)
- **Smart Primer Detection**: Automatic primer detection and trimming for amplicon data
- **QIIME2 Integration**: Seamless integration with QIIME2 2024.10 for downstream analysis
- **Taxonomy Assignment**: GreenGenes2 taxonomy classification and phylogenetic tree generation

## Notice

If you are running the AmpliconPIP or ggCOMBO functions on a server, please ensure you allocate sufficient time for your task. The data download and taxonomy annotation steps can be very time-consuming.

## Installation

### Step 1: Create Conda Environment

Copy the **env.yml** to your device.
The env.yml includes all dependencies (QIIME2, vsearch, fastp, git, pip, sra-tools, seqkit, etc.), but not Meta2Data itself.

```bash
conda env create -f env.yml
conda activate env
```

### Step 2: Install Meta2Data

#### Option A: Install to a custom directory (Recommended for HPC / read-only conda)

If your conda environment is read-only (e.g., TYKKY containers on HPC), install to a writable directory:

```bash
pip install --no-deps --target /path/to/Meta2Data \
    git+https://github.com/LinyangSun/Meta2Data.git@main

export PATH="/path/to/Meta2Data/bin:$PATH"
```

`--no-deps` skips dependency installation since all required packages are already provided by the conda environment.

#### Option B: Install into conda environment

```bash
pip install --no-deps git+https://github.com/LinyangSun/Meta2Data.git@main
```

#### Option C: Development install (local clone)

```bash
git clone https://github.com/LinyangSun/Meta2Data.git
cd Meta2Data
pip install --no-deps -e .
```

#### Verify installation

```bash
Meta2Data --help
```

## Requirements

### System Requirements
- **OS**: Linux (tested on Ubuntu/CentOS)
- **Memory**: Minimum 8GB RAM, 16GB+ recommended for large datasets
- **Storage**: Varies by dataset size (SRA files could be large)


## Usage

Meta2Data provides several subcommands for different stages of the workflow:

```bash
Meta2Data <command> [options]

Available commands:
    MetaDL         Search keywords combination in NCBI and CNCB. Download and preclean metadata.
    AmpliconPIP    Download and process amplicon sequencing data based on user provided metadata.
    ggCOMBO        Merge amplicon datasets and assign taxonomy using GreenGenes2, SILVA, or GSR-DB.
    ShortreadsPIP  (In development)
```

---

### MetaDL: Metadata Download

Download metadata from NCBI and CNCB databases with parallel processing and checkpoint/resume capability.

**Two modes:**

| Mode | Required Options | Description |
|------|-----------------|-------------|
| BioProject ID Input | `-i`, `-o` | Provide a directory of BioProject ID files |
| Keyword Search | `-o`, `--keywords`, `--field`, `--organism` | Search NCBI by keywords |

```
Required:
    -o, --output DIR              Output directory

BioProject Input Mode:
    -i, --input DIR               Directory with BioProject ID txt files

Keyword Search Mode:
    --keywords                    Enable keyword search mode
    --field "term1" "term2"       Search field terms
    --organism "term1" ...        Organism terms
    --opt "term1" ...             Optional additional terms

Optional:
    -k, --api-key KEY             NCBI API key (enables 8 parallel workers)
    -w, --max-workers NUM         Max parallel workers (default: 8 with key, 3 without)
    -h, --help                    Show help
```

**Output files:**
- `all_metadata_merged.csv` - Final merged dataset (NCBI + CNCB)
- `ncbi_merged_*.csv` - NCBI, EBI & DDBJ data
- `cncb_combined.csv` - CNCB data
- `logs/metadl_v2_*.log` - Detailed execution log
- `checkpoints/download_state.json` - Resume state for interrupted runs

---

### AmpliconPIP: Amplicon Data Processing

Download SRA data and process amplicon sequencing data with provided metadata.

```
Required (unless --test is used):
    -m, --metadata FILE           Input metadata CSV file
    --col-bioproject NAME         Column name for BioProject in CSV
    --col-sra NAME                Column name for SRA accession in CSV

Optional:
    -o, --output DIR              Output directory
                                  Default: current directory (--test) / metadata dir (normal)
    -t, --threads INT             Total CPU threads (default: 4)
                                  Auto-split: per-dataset threads = threads / max-parallel
    --max-parallel INT            Datasets to process in parallel (default: 2)
    --test                        Run in test mode
                                  Without -m: use built-in test data
                                  With -m: subset metadata (2 SRA per BioProject)
    -h, --help                    Show help
```

**Metadata CSV format** (column names customizable via `--col-*`):
```csv
Bioproject,Run
PRJNA12345,SRR123456
PRJNA12345,SRR123457
PRJNA67890,SRR234567
```

**Processing pipeline:**
1. Download SRA data (via Aspera/FTP)
2. Detect sequencing platform (Illumina, PacBio, Ion Torrent, 454) and layout (single/paired-end)
3. Identify and trim primers
4. Platform-specific processing (DADA2 denoising for Illumina/Ion Torrent/PacBio)
5. Generate QIIME2 artifacts (`.qza` files)

**Output structure:**
```
<output_dir>/
├── datasets_ID.txt                    # Generated dataset list
├── <dataset_ID>/                      # One directory per dataset
│   ├── <dataset_ID>_sra.txt          # SRA accession list
│   ├── ori_fastq/                     # Downloaded FASTQ files
│   ├── <dataset_ID>-final-rep-seqs.qza
│   └── <dataset_ID>-final-table.qza
├── failed_datasets.log
├── success_datasets.log
├── skipped_datasets.log
└── summary.csv                        # Sequencing depth for raw reads / final data
```

---

### ggCOMBO: Merge & Taxonomy Assignment

Merge multiple AmpliconPIP dataset outputs, assign taxonomy, and build a phylogenetic tree. Supports three taxonomy databases.

```
Required:
    --db DIR                      Path to database directory
    -i, --input DIR               Input directory containing PRJ* dataset folders
                                  (output of a previous AmpliconPIP run)

Optional:
    --db-type TYPE                Taxonomy database (default: greengenes)
                                    greengenes - GreenGenes2 2024.09
                                    silva      - SILVA 138.99
                                    gsr        - GSR-DB (Gut-Specific Reference DB)
    -o, --output DIR              Output directory (default: same as --input)
    -t, --threads INT             CPU threads (default: 4)
    --confidence FLOAT            Classifier confidence threshold (default: 0.7)
    --dl                          Download database files to --db directory
    -h, --help                    Show help
```

**Processing pipeline:**
1. Merge feature tables and representative sequences across all datasets
2. Orient sequences against database reference sequences
3. Filter feature table to oriented sequences only
4. Assign taxonomy via pre-trained Naive Bayes classifier
5. Build phylogenetic tree via SEPP fragment insertion (always GG2 reference)
6. Filter features by tree placement

**Output structure** (`DB_LABEL` = gg2, silva, or gsr):
```
<output_dir>/
└── final/
    └── merged/
        ├── merged-table.qza                 # Merged feature table (shared)
        ├── merged-rep-seqs.qza              # Merged representative sequences (shared)
        ├── merged-table-summary.qzv         # Table summary (shared)
        └── <DB_LABEL>/                      # Database-specific subfolder
            ├── oriented-rep-seqs.qza        # Oriented representative sequences
            ├── merged-table-oriented.qza    # Table filtered to oriented features
            ├── merged-taxonomy.qza          # Taxonomy classification
            ├── insertion-tree.qza           # SEPP phylogenetic tree
            ├── merged-table-tree.qza        # Table filtered to tree-placed features
            └── merged-table-no-tree.qza     # Features not placed in tree
```

---

## Examples

### Case 1: Complete Workflow — From Keyword Search to Taxonomy

Search NCBI for gut microbiome 16S studies, download & process data, then assign taxonomy.

```bash
conda activate Meta2Data

# Step 1: Search and download metadata by keywords
Meta2Data MetaDL \
    -o metadata/ \
    --keywords \
    --field "16S rRNA" "amplicon" \
    --organism "gut microbiome" \
    --opt "Illumina"

# Step 2: Process amplicon data
Meta2Data AmpliconPIP \
    -m metadata/all_metadata_merged.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
    -o results/ \
    -t 16

# Step 3: Merge and assign taxonomy (GreenGenes2, auto-download database)
Meta2Data ggCOMBO \
    --db ~/databases/gg2 \
    --dl \
    -i results/ \
    -t 16
```

### Case 2: Complete Workflow — From BioProject IDs to Taxonomy

You already have a list of BioProject IDs and want to process them end-to-end.

```bash
conda activate Meta2Data

# Step 1: Download metadata from a folder of BioProject ID files
Meta2Data MetaDL \
    -i bioproject_ids/ \
    -o metadata/ \
    -k YOUR_NCBI_API_KEY

# Step 2: Process amplicon data with 4 parallel datasets
Meta2Data AmpliconPIP \
    -m metadata/all_metadata_merged.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
    -o results/ \
    -t 32 --max-parallel 4

# Step 3: Assign taxonomy with SILVA database
Meta2Data ggCOMBO \
    --db ~/databases/silva \
    --db-type silva \
    --dl \
    -i results/ \
    -t 16
```

### Case 3: Process Amplicon Data Only (Metadata Already Prepared)

Skip the MetaDL step when you already have a metadata CSV file ready.

```bash
conda activate Meta2Data

# Custom column names matching your CSV headers
Meta2Data AmpliconPIP \
    -m my_samples.csv \
    --col-bioproject "ProjectID" \
    --col-sra "SRA_Accession" \
    -o amplicon_output/ \
    -t 8
```

### Case 4: Test Mode — Quick Validation

Verify the pipeline works before running on your full dataset.

```bash
conda activate Meta2Data

# Use built-in test data (fastest way to verify installation)
Meta2Data AmpliconPIP --test -t 8

# Or test with your own metadata (subsets to 2 SRA per BioProject)
Meta2Data AmpliconPIP --test \
    -m my_large_metadata.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
    -o test_output/ \
    -t 8
```

### Case 5: Taxonomy Assignment Only (AmpliconPIP Already Complete)

Run ggCOMBO independently on existing AmpliconPIP results, e.g., to compare databases.

```bash
conda activate Meta2Data

# GreenGenes2 (default)
Meta2Data ggCOMBO \
    --db ~/databases/gg2 \
    --db-type greengenes \
    --dl \
    -i amplicon_output/ \
    -o results_gg2/ \
    -t 16

# SILVA 138.99
Meta2Data ggCOMBO \
    --db ~/databases/silva \
    --db-type silva \
    --dl \
    -i amplicon_output/ \
    -o results_silva/ \
    -t 16

# GSR-DB (gut-specific, with lower confidence threshold)
Meta2Data ggCOMBO \
    --db ~/databases/gsr \
    --db-type gsr \
    --dl \
    --confidence 0.5 \
    -i amplicon_output/ \
    -o results_gsr/ \
    -t 16
```

### Case 6: Keyword Search with Multiple Terms

Search with broad keyword sets covering multiple organisms and fields.

```bash
conda activate Meta2Data

# Define keyword arrays
field=("Bacteria" "Microbiome" "Metagenomics" "Metabarcoding")
organism=("bee" "apis" "bombus")
opt=("Amplicon" "16s" "gut")

Meta2Data MetaDL \
    -o bee_metadata/ \
    --keywords \
    --field "${field[@]}" \
    --organism "${organism[@]}" \
    --opt "${opt[@]}" \
    -k YOUR_NCBI_API_KEY
```

### Case 7: HPC / Server Usage Tips

For long-running jobs on HPC systems, ensure sufficient time allocation.

```bash
conda activate Meta2Data

# Large-scale processing: increase threads and parallelism
Meta2Data AmpliconPIP \
    -m metadata/all_metadata_merged.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
    -o /scratch/results/ \
    -t 64 --max-parallel 8

# Run multiple taxonomy databases sequentially
Meta2Data ggCOMBO \
    --db ~/databases/gg2 \
    --db-type greengenes \
    --dl \
    -i /scratch/results/ \
    -t 16

Meta2Data ggCOMBO \
    --db ~/databases/silva \
    --db-type silva \
    --dl \
    -i /scratch/results/ \
    -t 16
# Results in /scratch/results/final/merged/gg2/ and .../silva/
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

