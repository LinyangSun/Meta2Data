# Meta2Data Dependencies

This document provides a comprehensive list of all dependencies required to run the full Meta2Data workflow.

## Quick Start

The recommended way to install all dependencies is via Conda:

```bash
conda env create -f env.yml
conda activate Meta2Data
pip install -e .
```

## Python (>= 3.10)

### Core Python Packages (`pyproject.toml`)

| Package | Version | Purpose |
|---|---|---|
| biopython | >= 1.84 | NCBI Entrez API access, sequence parsing |
| pandas | >= 2.2.2 | Metadata manipulation and CSV/Excel handling |
| openpyxl | >= 0.0 | Excel file read/write support |
| requests | >= 2.32.3 | HTTP requests (CNCB/GSA API, ENA) |
| numpy | == 1.26.4 | Numerical computation |
| biom-format | >= 2.1.15 | BIOM feature table format support |

## Bioinformatics Tools

### QIIME2 2024.10 (Amplicon Distribution)

Core workflow orchestration framework. The following plugins are required:

| Plugin | Purpose |
|---|---|
| q2-demux | Demultiplexing |
| q2-dada2 | DADA2 denoising |
| q2-deblur | Deblur denoising |
| q2-cutadapt | Primer trimming |
| q2-quality-filter | Quality filtering |
| q2-feature-table | Feature table operations |
| q2-feature-classifier | Taxonomy classification |
| q2-alignment | Sequence alignment (MAFFT) |
| q2-phylogeny | Phylogenetic tree building |
| q2-fragment-insertion | SEPP fragment insertion |
| q2-diversity / q2-diversity-lib | Alpha/beta diversity |
| q2-taxa | Taxonomy visualization |
| q2-vsearch | Clustering and chimera detection |
| q2-composition | Compositional data analysis |
| q2-emperor | 3D PCoA visualization |
| q2-longitudinal | Longitudinal analysis |
| q2-metadata | Metadata handling |
| q2-quality-control | Quality control |
| q2-sample-classifier | Machine learning classifiers |
| q2-stats | Statistical tests |
| q2-vizard | Visualization |
| q2cli | QIIME2 command-line interface |
| q2templates / q2galaxy | Template and Galaxy support |

### Sequence Quality Control & Processing

| Tool | Version | Purpose |
|---|---|---|
| fastp | any | Adapter removal and quality control |
| cutadapt | >= 4.9 | Primer trimming |
| vsearch | >= 2.22.1 | Clustering, chimera detection, dereplication |
| SeqKit | any | Sequence manipulation toolkit |

### Sequence Alignment & Taxonomy

| Tool | Version | Purpose |
|---|---|---|
| BLAST | >= 2.16.0 | Sequence alignment / similarity search |
| MAFFT | >= 7.526 | Multiple sequence alignment |
| MUSCLE | >= 3.8.1551 | Sequence alignment |
| ClustalW | >= 2.1 | Multiple sequence alignment |

### Phylogenetic Tree Building

| Tool | Version | Purpose |
|---|---|---|
| FastTree | >= 2.1.11 | Fast phylogenetic tree inference |
| RAxML | >= 8.2.13 | Maximum likelihood phylogenetics |
| IQ-TREE | >= 2.3.6 | Phylogenetic tree building |
| SEPP | >= 4.5.5 | Fragment insertion phylogenetic placement |

### Data Download & SRA Processing

| Tool | Version | Purpose |
|---|---|---|
| sra-tools | any | NCBI SRA data download (fasterq-dump, prefetch) |
| entrez-direct | >= 22.4 | NCBI Entrez command-line utilities (esearch, efetch) |

### Other Bioinformatics Tools

| Tool | Version | Purpose |
|---|---|---|
| Bowtie2 | >= 2.5.4 | Sequence alignment |
| SAMtools | >= 1.21 | BAM/SAM file processing |
| HMMER | >= 3.4 | HMM protein/nucleotide sequence analysis |
| SortMeRNA | >= 2.0 | rRNA sequence filtering |
| Deblur | >= 1.1.1 | Denoising algorithm |
| UniFrac | >= 1.3 | Phylogenetic distance metrics |

## R / Bioconductor (via Conda)

| Package | Version | Purpose |
|---|---|---|
| r-base | >= 4.3.3 | R runtime |
| bioconductor-dada2 | >= 1.30.0 | DADA2 denoising engine |
| bioconductor-phyloseq | >= 1.46.0 | Microbiome data analysis |
| bioconductor-decipher | >= 2.30.0 | Sequence alignment |
| bioconductor-decontam | >= 1.22.0 | Contamination removal |
| bioconductor-ancombc | >= 2.4.0 | Differential abundance analysis |
| bioconductor-mia | >= 1.10.0 | Microbiome analysis |
| bioconductor-biostrings | >= 2.70.1 | Biological sequence classes |

See `env.yml` for the full list of ~200+ R/Bioconductor packages included.

## System Libraries

| Library | Version | Purpose |
|---|---|---|
| libffi | >= 3.4.2 | Foreign function interface |
| libcurl | >= 8.8.0 | URL retrieval |
| libxml2 | >= 2.12.7 | XML parsing |
| libxslt | >= 1.1.39 | XSLT processing |
| zlib | >= 1.2.13 | Compression |
| bzip2 | >= 1.0.8 | Compression |
| HDF5 | >= 1.14.3 | Hierarchical data format |
| OpenSSL | >= 3.3.2 | Cryptography / SSL |
| SQLite | >= 3.46.0 | Lightweight database |

## System Utilities

| Tool | Purpose |
|---|---|
| GCC / GFortran / Make | Build toolchain |
| wget / curl | File and data download |
| git | Version control |
| jq | JSON query processor |
| pigz / pbzip2 | Parallel compression |
| sed / Perl | Text processing |

## External Databases & APIs

### Databases (downloaded at runtime)

| Database | File | Purpose |
|---|---|---|
| GreenGenes2 (GG2) | `2024.09.backbone.full-length.nb.sklearn-1.4.2.qza` | Pre-trained Naive Bayes taxonomy classifier |
| GG2 SEPP reference | `sepp-refs-gg-13-8.qza` | SEPP reference tree for fragment insertion |

Use `Meta2Data ggCOMBO --dl` to download these automatically.

### External APIs

| API | Purpose |
|---|---|
| NCBI Entrez | BioProject/BioSample metadata retrieval |
| CNCB/GSA (ngdc.cncb.ac.cn) | Chinese sequence data retrieval |
| ENA (European Nucleotide Archive) | Fallback FASTQ download |
| NCBI SRA FTP | SRA file download |

## Dependencies by Subcommand

| Subcommand | Key Dependencies |
|---|---|
| **MetaDL** | biopython, requests, pandas, openpyxl, NCBI Entrez API, CNCB API |
| **AmpliconPIP** | QIIME2 full suite, fastp, cutadapt, vsearch, sra-tools, entrez-direct, SeqKit, DADA2/Deblur |
| **ggCOMBO** | QIIME2, q2-feature-classifier, q2-fragment-insertion, SEPP, GreenGenes2 database |

## System Requirements

| Resource | Minimum | Recommended |
|---|---|---|
| OS | Linux (Ubuntu/CentOS) | Linux |
| Python | 3.10 | 3.10+ |
| RAM | 8 GB | 16 GB+ |
| Storage | 50 GB free | Varies by dataset |
| CPU | Any | Multi-core (parallelism via `-t` flag) |
