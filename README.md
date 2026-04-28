# Meta2Data

**Automated Bioinformatics Pipeline for microbiome Sequencing Data Processing from public Databases**

Meta2Data is a command-line tool for downloading, processing, and analyzing metabarcoding data (maybe also include metagenome in future) from public databases (NCBI, CNCB/GSA). It integrates metadata retrieval, SRA data download, quality control, and QIIME2-based analysis into a single, automated workflow.

## 3 module
1. MetaDL >> For metadata preprocessing.
2. AmpliconPIP >> For sequencing data processing.
3. ggCOMBO >> For taxonomy annotation.

## Features

- **Metadata Download and Pre-clean**: (MetaDL module) Search, download, and pre-clean metadata from NCBI and CNCB databases by keywords, BioProject ID, or BioSample ID. Auto-fetches BioProject descriptions and standardizes column names.
- **Multi-Platform Support**: (AmpliconPIP module) Automatic detection and processing of Illumina, PacBio, Ion Torrent, and 454 sequencing platforms (ONT not supported).
- **Smart Primer Detection**: (AmpliconPIP module) Automatic primer detection and trimming for amplicon data (no need provide primer detail).
- **QIIME2 Integration**: (AmpliconPIP module) Integration with QIIME2 2024.10 for downstream analysis.
- **Taxonomy Assignment**: (ggCOMBO module) Taxonomy classification (GreenGenes2, Silva and gsr supported) and phylogenetic tree generation.
- **OS**: Only for linux.
- **Others**: Parallel task supported for AmpliconPIP.

## Notice

Please ensure you allocate sufficient time for your task (1-2 days). The data download and phylogenetic tree generation steps can be very time-consuming.
The ggCOMBO do not support parallel task, you need to run AmpliconPIP and ggCOMBO at seperate task, if you prefer Parallelize the AmpliconPIP task

## Installation

Meta2Data can be installed in a local folder to avoid contaminating your QIIME2 environment and to make updates easier. Expose it to your shell by adding <repo>/bin to your PATH. QIIME2 is only required if you plan to run AmpliconPIP or ggCOMBO; MetaDL runs on any Python 3 interpreter.

### Step 1: Clone the repo and add it to PATH

```bash
git clone https://github.com/LinyangSun/Meta2Data.git
echo 'export PATH="'"$PWD"'/Meta2Data/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```
> If you also use conda, the order of lines in your rc file matters.
> Place the export PATH=$HOME/Meta2Data/bin:$PATH line before any conda activate command.


### Step 2: Install QIIME2 (only for AmpliconPIP / ggCOMBO)

See https://docs.qiime2.org/2024.10/install/ for the current recommended procedure. 

### Step 4: Verify installation

```bash
Meta2Data --help
```

## Requirements

### System Requirements
- **OS**: Linux (tested on Ubuntu/CentOS)
- **Computation Resources**: For AmpliconPIP, 8GB RAM, 4 CPU are recommended per parallel task. For ggCOMBO, 70-100GB RAM and 20CPU are recommonded



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

Download metadata from NCBI and CNCB databases with parallel processing and checkpoint/resume capability. Automatically fetches BioProject descriptions and standardizes column names (CamelCase normalization, synonym merging via dictionary).

**Three modes:**

| Mode | Required Options | Description |
|------|-----------------|-------------|
| ID Input | `-i`, `-o` | Provide a directory of txt files containing BioProject IDs (PRJ*) and/or BioSample IDs (SAM*) — they can be mixed |
| Keyword Search | `-o`, `--keywords`, `--field`, `--organism` | Search NCBI + CNCB by keywords, then download metadata for matched BioProjects |

```
Required:
    -o, --output DIR              Output directory

ID Input Mode:
    -i, --input DIR               Directory with ID txt files (BioProject and/or BioSample)

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

**Output columns:** `Run, Bioproject, Description, DesignDescription, Biosample, Experiment, ...` (core columns first, then rare columns grouped alphabetically)

**Output files:**

Main outputs:
- `all_metadata_merged.csv` — Final merged dataset with auto-fetched BioProject descriptions and SRA experiment design descriptions. Column names are standardized (CamelCase normalization + synonym merging via dictionary).
- `status.tsv` — Processing status for each input ID (`has_data` / `no_data` / `no_run_info` / `download_error`)
- `column_description.tsv` — Per-column statistics: fill rate, number of datasets covered, top 5 values, and column type (`core` / `cncb` / `rare`)
- `RecordWithoutRUNinfo.csv` — Records without SRA Run info (only generated when such records exist)

Keyword search mode only:
- `searched_keywords/combined_results.csv` — BioProject IDs matched by keyword search
- `searched_keywords/search_summary.txt` — Query summary

Internal / resume:
- `tmp/checkpoints/download_state.json` — Resume state for interrupted runs

> **Tip 1 — AI-assisted metadata screening**
>
> MetaDL automatically fetches a `Description` column for every BioProject, which summarizes each study's purpose, target organism, and experimental design. This makes the merged CSV well-suited for AI-based screening. Upload the CSV to any AI tool with a Team/collaborative workspace (Claude Team, ChatGPT Team, Gemini, etc.) and ask it to filter based on your criteria. For example:
>
> *"Here is my metadata CSV. Based on the Description and other columns, keep only gut microbiome samples from healthy human adults sequenced on Illumina with 16S amplicons. Filter the data and give me the result. You need to use 3 workers to screen the full data independently, and a leader to give a final decision. Then output with xxx xxx xxx files. Ask me anything unclear to you before starting."*

> **Tip 2 — Broader keyword search with genus-level terms**
>
> Keyword search results depend heavily on how authors annotate their BioProjects. To maximize coverage for a taxonomic group of interest, don't rely solely on high-level terms (e.g., "bee"). Instead, collect genus names from a published phylogeny or species tree for your clade (for example from tree of life), and include them as `--organism` terms. For example:
>
> ```bash
> # Instead of just "bee", also search by genus names from the Apoidea phylogeny
> organism=("bee" "Apis" "Bombus" "Megachile" "Osmia" "Andrena" "Halictus")
> Meta2Data MetaDL \
>     -o metadata/ \
>     --keywords \
>     --field "16S rRNA" "amplicon" \
>     --organism "${organism[@]}"
> ```
>
> This catches studies that only mention a genus in their BioProject metadata and would otherwise be missed.

> **Tip 3 — Quick column labeling with `column_description.tsv`**
>
> Download `column_description.tsv` and open it in Excel. Add a new column (e.g., `Label`) and tag each row as `keep`, `drop`, or `body part` ... based on the fill rate and top values. This gives you a quick overview of all available columns and makes subsequent data cleaning much faster — you can filter by your labels to decide which columns to retain before any downstream analysis.

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
3. Quality control
4. Identify and trim primers
5. Platform-specific processing (DADA2 denoising for Illumina/Ion Torrent/PacBio, Vsearch for 454)
6. Generate QIIME2 artifacts (`.qza` files)

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

### Case 1: Different way to run metaDL

```bash
conda activate qiime2-amplicon-2024.10

# Step 1: Search and download metadata by keywords

field=("16S rRNA" "amplicon")
organism=("bee" "bees")

Meta2Data MetaDL \
    -o metadata/ \
    --keywords \
    --field "${field[@]}" \
    --organism "${organism[@]}" \
    --opt "Illumina"


# Step 1: Download metadata from a folder of BioProject ID files
Meta2Data MetaDL \
    -i bioproject_ids/ \
    -o metadata/ 
```

### Case 2: Process Amplicon Data Only (Metadata Already Prepared)

Skip the MetaDL step when you already have a metadata CSV file ready.

```bash
conda activate qiime2-amplicon-2024.10

# Custom column names matching your CSV headers
Meta2Data AmpliconPIP \
    -m my_samples.csv \ # your metadata.csv
    --col-bioproject "ProjectID" \ # column name for bioproject
    --col-sra "SRA_Accession" \ # column name for sra run (the sra normally start with SRR, ERR, DRR OR CRR)
    -o amplicon_output/ \
    -t 8
```

### Case 3: Taxonomy Assignment Only (AmpliconPIP Already Complete)

Run ggCOMBO independently on existing AmpliconPIP results, e.g., to compare databases.

```bash
conda activate qiime2-amplicon-2024.10

# GreenGenes2
Meta2Data ggCOMBO \
    --db path/to/your/metafile/databases/ \ # Prepare an empty path, the pip will download database in it
    --db-type greengenes \
    --dl \
    --confidence 0.7 \
    -i path/to/your/metafile/amplicon_output/ \
    -t 16

# SILVA 138.99
Meta2Data ggCOMBO \
    --db path/to/your/metafile/databases/ \
    --db-type silva \
    --dl \
    --confidence 0.7 \
    -i path/to/your/metafile/amplicon_output/ \
    -t 16

# GSR-DB (gut-specific, with lower confidence threshold)
Meta2Data ggCOMBO \
    --db path/to/your/metafile/databases/ \
    --db-type gsr \
    --dl \
    --confidence 0.7 \
    -i amplicon_output/ \
    -t 16
```



### Case 4: Test Mode — Quick Validation

Verify the pipeline works before running on your full dataset. Every
`AmpliconPIP --test` invocation starts with a **pre-flight dependency
check** that inspects vendor binaries (`vsearch`, `fastp`), system
utilities, QIIME2 plugins, and Python packages, and aborts with a per-layer
report and exact install commands if anything is missing. Fix the missing
items and re-run — the check itself is idempotent and costs ~1 second.

```bash
conda activate qiime2-amplicon-2024.10

# Use built-in test data (fastest way to verify installation).
# The dependency check runs first; if everything is satisfied you'll see
# `[check] All dependencies satisfied.` before the test pipeline starts.
Meta2Data AmpliconPIP --test -t 8

# Or test with your own metadata (subsets to 2 SRA per BioProject)
Meta2Data AmpliconPIP --test \
    -m path/to/your/metafile/metadata.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
    -o test_output/ \
    -t 8

# Run the dependency check alone, without launching the pipeline.
bash scripts/check_dependencies.sh
```



> **Note on `--db-type gsr`** — GSR-DB's upstream ships a pre-trained
> classifier pickled against an older scikit-learn and cannot be trusted with
> QIIME2 2024.10's sklearn 1.4.2. ggCOMBO therefore downloads only the raw
> reference tarball (seqs + taxa) and trains a fresh Naive Bayes classifier
> for information of GSR, see: https://journals.asm.org/doi/10.1128/msystems.00950-23



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

