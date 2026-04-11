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

Meta2Data is a self-contained command-line tool. It lives **entirely inside
the cloned repository folder** — there is no `pip install` step and no
conda environment to create for Meta2Data itself. You expose it to your
shell by adding `<repo>/bin` to your `PATH`. QIIME2 is only required if
you plan to run `AmpliconPIP` or `ggCOMBO`; `MetaDL` runs on any Python 3
interpreter.

### Step 1: Clone the repo and add it to PATH

```bash
git clone https://github.com/LinyangSun/Meta2Data.git
echo 'export PATH="'"$PWD"'/Meta2Data/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

> **Order matters if you also use conda.** Put the `export PATH` line
> *before* any `conda activate` in your shell session or rc file. Conda
> prepends its env's `bin` to `PATH` on activation, so putting Meta2Data
> first in rc means conda later shadows it with its own `python`/`pip` —
> which is exactly what you want, and Meta2Data's own commands remain
> resolvable via the persistent `$HOME/Meta2Data/bin` entry.

### Step 2: Download vsearch and fastp (automatic on first run)

`AmpliconPIP` and `ggCOMBO` need two native binaries that QIIME2 does not
ship: `vsearch 2.30.0` and `fastp 0.24.0`. **You do not need to download
them manually** — on the first `AmpliconPIP` or `ggCOMBO` invocation,
Meta2Data detects they are missing and auto-runs
`scripts/install_binaries.sh`, which fetches prebuilt Linux static binaries
from upstream release pages (no compilation, no extra conda packages) into
`<repo>/vendor/bin/`. Subsequent runs find them in place and skip the
download.

If you prefer to download them eagerly (e.g., to verify the first run
works without an extra network round-trip, or because the target machine
is offline later), you can run the installer yourself:

```bash
cd Meta2Data
bash scripts/install_binaries.sh     # optional — otherwise runs automatically
```

The script is idempotent — re-run it anytime without side effects. If you
only intend to use `MetaDL`, you can skip this step entirely (MetaDL does
not need vsearch or fastp). Use `bash scripts/install_binaries.sh --help`
for `--prefix` and `--force`. Set `META2DATA_SKIP_DEP_CHECK=1` to disable
the automatic check on subcommand startup.

### Step 3: Install QIIME2 (only for AmpliconPIP / ggCOMBO)

`AmpliconPIP` and `ggCOMBO` invoke `qiime` plugins at runtime. Follow
QIIME2's official instructions to install the **2024.10 amplicon
distribution** and activate it before running those subcommands:

```bash
conda activate qiime2-amplicon-2024.10
```

See https://docs.qiime2.org/2024.10/install/ for the current recommended
procedure. **`MetaDL` does not need this step** — it can run against any
`python3` on your `PATH` (system, venv, or a minimal conda env with just
`python3` installed).

### Python packages — handled automatically

Meta2Data uses `biopython`, `pandas`, `numpy`, and `requests`. On the
first run of any subcommand, the entry script detects missing Python
packages against the currently active `python3` and auto-installs them
with `python3 -m pip install`. You do **not** need to install them
manually.

The install lands in whichever `python3` is first on your `PATH`:

- If you activated the QIIME2 env, packages land in the QIIME2 env.
- If you did not activate anything, they land in the system / user
  Python — same as any plain `pip install` would.
- If you made a dedicated venv for `MetaDL`, they land there.

To opt out (e.g., on a read-only HPC login node where pip will fail
anyway), set `META2DATA_SKIP_DEP_CHECK=1` and manage the packages
yourself.

### Step 4: Verify installation

```bash
# MetaDL-only users
Meta2Data MetaDL --help

# Full pipeline users (QIIME2 env active)
Meta2Data --help
bash scripts/check_dependencies.sh
```

`check_dependencies.sh` inspects four layers — pipeline binaries, system
utilities, QIIME2 plugins, and Python packages — and prints a per-layer
report of anything missing with the exact install command to fix it. The
same check runs automatically at the top of every `AmpliconPIP --test`
invocation.

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
conda activate qiime2-amplicon-2024.10

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
conda activate qiime2-amplicon-2024.10

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
conda activate qiime2-amplicon-2024.10

# Custom column names matching your CSV headers
Meta2Data AmpliconPIP \
    -m my_samples.csv \
    --col-bioproject "ProjectID" \
    --col-sra "SRA_Accession" \
    -o amplicon_output/ \
    -t 8
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
    -m my_large_metadata.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
    -o test_output/ \
    -t 8

# Run the dependency check alone, without launching the pipeline.
bash scripts/check_dependencies.sh
```

### Case 5: Taxonomy Assignment Only (AmpliconPIP Already Complete)

Run ggCOMBO independently on existing AmpliconPIP results, e.g., to compare databases.

```bash
conda activate qiime2-amplicon-2024.10

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

> **Note on `--db-type gsr`** — GSR-DB's upstream ships a pre-trained
> classifier pickled against an older scikit-learn and cannot be trusted with
> QIIME2 2024.10's sklearn 1.4.2. ggCOMBO therefore downloads only the raw
> reference tarball (seqs + taxa) and trains a fresh Naive Bayes classifier
> **locally** on first `--dl` run via
> `qiime feature-classifier fit-classifier-naive-bayes`.
>
> - **First run**: add ~15–40 minutes and up to ~20 GB peak RAM for the one-time
>   training step. Don't interrupt it — partial artifacts are cleaned up on
>   error but you'll restart from scratch.
> - **Subsequent `--dl` runs**: the existing `classifier_GSR-DB_full-16S.qza`
>   is validated and reused. Zero extra cost.
> - **QIIME2 / sklearn upgrades**: if you later upgrade to a QIIME2 release
>   that bumps scikit-learn to a new major version (e.g. 1.5.x+), delete the
>   cached classifier and re-run `--dl` to retrain against the new sklearn:
>   ```bash
>   rm ~/databases/gsr/classifier_GSR-DB_full-16S.qza
>   Meta2Data ggCOMBO --db ~/databases/gsr --db-type gsr --dl -i ... -t 16
>   ```
>
> `--db-type greengenes` and `--db-type silva` are unaffected: their upstream
> classifiers are already version-matched to QIIME2 2024.10.

### Case 6: Keyword Search with Multiple Terms

Search with broad keyword sets covering multiple organisms and fields.

```bash
conda activate qiime2-amplicon-2024.10

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
conda activate qiime2-amplicon-2024.10

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
├── bin/                        # Command-line entry points
│   ├── Meta2Data               # Main dispatcher
│   ├── Meta2Data-MetaDL        # Metadata download
│   ├── Meta2Data-AmpliconPIP   # Amplicon pipeline
│   └── Meta2Data-ggCOMBO       # Merge & taxonomy assignment
├── scripts/                    # Processing scripts
│   ├── AmpliconFunction.sh     # Amplicon processing functions
│   ├── run.sh                  # AmpliconPIP orchestrator
│   ├── taxonomy.sh             # ggCOMBO taxonomy pipeline
│   ├── metadata_downloader.py  # Python metadata downloader
│   ├── py_16s.py               # 16S data processing
│   ├── install_binaries.sh     # Downloads vsearch + fastp native binaries
│   └── check_dependencies.sh   # Pre-flight dependency check
├── vendor/                     # (gitignored) vsearch/fastp from install_binaries.sh
├── test/                       # Test data
└── docs/                       # Documentation & reference sequences
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

