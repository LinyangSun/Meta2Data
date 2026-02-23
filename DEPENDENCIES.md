# Meta2Data Dependencies

This document lists all dependencies **actually called** in the current Meta2Data workflow. Only tools and libraries that are invoked during execution are included — packages that exist in `env.yml` but are never called by the active codebase are omitted.

## Quick Start

```bash
conda env create -f env.yml
conda activate Meta2Data
pip install -e .
```

## Python (>= 3.10)

### Core Python Packages (`pyproject.toml`)

| Package | Version | Purpose |
|---|---|---|
| biopython | >= 1.84 | NCBI Entrez API access (`Bio.Entrez`), sequence parsing |
| pandas | >= 2.2.2 | Metadata manipulation and CSV/Excel handling |
| openpyxl | >= 0.0 | Excel file read/write support |
| requests | >= 2.32.3 | HTTP requests (CNCB/GSA API) |
| numpy | == 1.26.4 | Numerical computation |
| biom-format | >= 2.1.15 | BIOM feature table format support |

## Bioinformatics Tools

### Direct CLI Calls

These tools are invoked directly as shell commands in the active scripts:

| Tool | Called In | Purpose |
|---|---|---|
| **fastp** | `run.sh` | Adapter removal and quality control for FASTQ files |
| **wget** | `AmpliconFunction.sh` | Downloading FASTQ from ENA and CNCB |
| **gzip** | `AmpliconFunction.sh` | Verifying download integrity (`gzip -t`) |

### QIIME2 2024.10 (Amplicon Distribution)

Core workflow orchestration framework. The following plugins are **actually invoked** in the active codebase:

| Plugin | Called In | Invocation |
|---|---|---|
| q2cli | `AmpliconFunction.sh`, `taxonomy.sh` | `qiime tools import`, `qiime tools export` |
| q2-demux | `AmpliconFunction.sh` | `qiime demux summarize` |
| q2-dada2 | `AmpliconFunction.sh` | `qiime dada2 denoise-paired`, `denoise-single`, `denoise-pyro`, `denoise-ccs` |
| q2-quality-filter | `AmpliconFunction.sh` | `qiime quality-filter q-score` |
| q2-feature-table | `AmpliconFunction.sh`, `taxonomy.sh` | `qiime feature-table merge`, `merge-seqs`, `filter-features`, `filter-seqs`, `summarize` |
| q2-feature-classifier | `AmpliconFunction.sh`, `taxonomy.sh` | `qiime feature-classifier classify-sklearn`, `extract-reads` |
| q2-fragment-insertion | `taxonomy.sh` | `qiime fragment-insertion sepp`, `filter-features` |
| q2-vsearch | `AmpliconFunction.sh` | `qiime vsearch dereplicate-sequences`, `uchime-denovo`, `cluster-features-de-novo` |

**Note:** Tools like vsearch, cutadapt, MAFFT, SEPP, DADA2 (R), and sklearn are called **indirectly** through the QIIME2 plugins above — they are runtime dependencies of QIIME2, not invoked directly by Meta2Data scripts.

### NOT Called (present in env.yml but unused by active code)

The following tools are included in `env.yml` (as QIIME2 distribution dependencies) but are **not invoked** by any active Meta2Data script — they only appear in comments, backup files, or are never referenced:

- ~~sra-tools~~ (prefetch, fasterq-dump) — data download now uses wget + ENA/CNCB
- ~~entrez-direct~~ (esearch/efetch CLI) — NCBI access uses Bio.Entrez Python API
- ~~SeqKit~~ — only in backup scripts
- ~~BLAST~~ (blastn/blastp direct calls) — only in backup scripts
- ~~MUSCLE~~, ~~ClustalW~~ — never referenced
- ~~FastTree~~, ~~RAxML~~, ~~IQ-TREE~~ — never directly called
- ~~Bowtie2~~, ~~SAMtools~~, ~~HMMER~~, ~~SortMeRNA~~ — never referenced
- ~~jq~~, ~~pigz~~ — never referenced

## R / Bioconductor (via Conda)

R and Bioconductor packages are **runtime dependencies of QIIME2 plugins** (especially q2-dada2). They are not called directly by Meta2Data scripts but are required for DADA2 denoising to function:

| Package | Version | Required By |
|---|---|---|
| r-base | >= 4.3.3 | q2-dada2 |
| bioconductor-dada2 | >= 1.30.0 | q2-dada2 denoising engine |

## System Libraries

These are transitive dependencies required by Python packages and QIIME2:

| Library | Purpose |
|---|---|
| zlib | Compression (gzip operations) |
| OpenSSL | HTTPS connections (wget, requests) |

Other system libraries (libffi, libxml2, HDF5, etc.) are pulled in automatically by Conda as transitive dependencies.

## External Databases & APIs

### Databases (downloaded at runtime)

| Database | File | Required By |
|---|---|---|
| GreenGenes2 (GG2) | `2024.09.backbone.full-length.nb.sklearn-1.4.2.qza` | `taxonomy.sh` — Naive Bayes taxonomy classifier |
| GG2 SEPP reference | `sepp-refs-gg-13-8.qza` | `taxonomy.sh` — SEPP fragment insertion tree |

Use `Meta2Data ggCOMBO --dl` to download these automatically.

### External APIs (accessed via Python)

| API | Accessed Via | Purpose |
|---|---|---|
| NCBI Entrez | Bio.Entrez (biopython) | BioProject/BioSample metadata retrieval |
| CNCB/GSA | requests + wget | Chinese sequence data retrieval |
| ENA (European Nucleotide Archive) | wget | FASTQ file download |

## Dependencies by Subcommand

| Subcommand | Key Dependencies |
|---|---|
| **MetaDL** | biopython (Bio.Entrez), requests, pandas, openpyxl |
| **AmpliconPIP** | fastp, wget, QIIME2 (q2-dada2, q2-vsearch, q2-quality-filter, q2-demux, q2-feature-table, q2-feature-classifier) |
| **ggCOMBO** | QIIME2 (q2-feature-table, q2-feature-classifier, q2-fragment-insertion), GreenGenes2 database, wget |

## System Requirements

| Resource | Minimum | Recommended |
|---|---|---|
| OS | Linux (Ubuntu/CentOS) | Linux |
| Python | 3.10 | 3.10+ |
| RAM | 8 GB | 16 GB+ |
| Storage | 50 GB free | Varies by dataset |
| CPU | Any | Multi-core (parallelism via `-t` flag) |
