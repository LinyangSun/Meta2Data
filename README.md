<p align="center">
  <img src="docs/llloooooog.png" alt="Meta2Data logo" width="100%">
</p>

# Meta2Data

**Automated Bioinformatics Pipeline for microbiome Sequencing Data Processing from public Databases**

Meta2Data is a command-line tool for downloading, processing, and analyzing metabarcoding data (maybe also include metagenome in future) from public databases (INSDC, CNCB/GSA). It integrates metadata retrieval, SRA data download, quality control, and QIIME2-based analysis into a single, automated workflow.

## 3 module
1. MetaDL >> For metadata preprocessing.
2. AmpliconPIP >> For sequencing data processing.
3. AmpliconTAXA >> For taxonomy annotation.

## Features

- **Metadata Download and Pre-clean**: (MetaDL module) Search, download, and pre-clean metadata from INSDC and CNCB databases by keywords, BioProject ID, or BioSample ID. Auto-fetches BioProject descriptions and standardizes column names.
- **Multi-Platform Support**: (AmpliconPIP module) Automatic detection and processing of Illumina, PacBio, Ion Torrent, 454, and Oxford Nanopore (ONT) sequencing platforms. The platform is detected automatically from INSDC/CNCB for each downloaded dataset (or set explicitly with `--platform` in local mode).
- **DADA2 / vsearch workflows**: (AmpliconPIP module) Choose the denoising strategy with a required `--dada2` / `--vsearch` flag. `--dada2` runs DADA2 for single-base ASV resolution (Illumina / Ion Torrent / PacBio CCS); `--vsearch` runs vsearch for 97%-identity OTUs and is robust on all 5 platforms, including degraded / binned-quality and ONT data.
- **Local Mode**: (AmpliconPIP module) Process FASTQ files you already have with `--local` — no download, no INSDC lookup. One run handles one dataset/one platform; you supply `--platform` and, optionally, explicit primers (`--primer-fwd`/`--primer-rev`, otherwise auto-detected). Your original files are never modified.
- **Smart Primer Detection**: (AmpliconPIP module) Automatic entropy-based primer detection and trimming for amplicon data (no need to provide primer details); explicit primers can be given in local mode.
- **Per-Dataset Summary & Region Detection**: (AmpliconPIP module) After processing, a `per_dataset_summary.tsv` reports each dataset's platform, quality status, and the amplified 16S V-region (e.g. `V3-V4`, `V1-V9`), inferred by aligning representative sequences to the E. coli 16S reference. Summaries and the unified status log are re-run safe (upserted / append-only).
- **QIIME2 Integration**: (AmpliconPIP module) Integration with QIIME2 2024.10 for downstream analysis.
- **Taxonomy Assignment**: (AmpliconTAXA module) Taxonomy classification (GreenGenes2 and SILVA supported) and phylogenetic tree generation, in ASV or OTU mode (must match the mode used in AmpliconPIP).
- **OS**: Only for linux.
- **Others**: Parallel task supported for AmpliconPIP.

## Notice

Please ensure you allocate sufficient time for your task (1-2 days). The data download and phylogenetic tree generation steps can be very time-consuming.

## Updates

### 2026-09 revision

- Permanent per-step read counts extend `summary.csv`, with `dada2_` / `vsearch_` branch columns, a BioProject total table, and TAXA retained/lost abundance tables. See [read-count outputs and units](docs/read_counts.md).
- Optional `M2D_PROFILE_DIR` records command time, CPU, memory observations, I/O and nested invocation history without changing `summary.csv`. Actual FASTQ bases/bytes/layout are retained with the raw-count audit. See [resource profiling](docs/resource_profile.md).

- Workflow flags and result names now use `dada2` and `vsearch`.
- Primer defaults follow benchmark b1: first 20 bp, database first, strict `fold < 16` fallback; unknown primers are trimmed by 20 bp unless `--skip-unknown-primers` is set.
- Local pairing, sample IDs, raw counts and manifests use one naming rule, including `_1_001` / `_2_001` and multi-lane files.
- TAXA discovers complete result pairs recursively, deduplicates artifacts, and selects the region workflow independently with `--singleV` or `--notree` (default: multi-region SEPP).
- `--parameter` loads validated JSON overrides. Effective settings and input provenance are saved, and changed inputs/settings invalidate affected checkpoints.

### 2026-06-12

- **Fix (AmpliconPIP `--local`): force `--max-parallel 1` in local mode.** A `--local` run always resolves to exactly one dataset, but `--max-parallel` defaulted to `2`, so the per-dataset thread budget was `threads ÷ 2` and half of the requested `-t` threads sat idle. Local runs now give the single dataset all `-t` threads (if a different `--max-parallel` is passed it is overridden, with a notice).
- **Fix (AmpliconPIP ASV / DADA2): run DADA2 denoising multithreaded.** All four `qiime dada2` denoise methods (`denoise-paired`, `denoise-pyro`, `denoise-single`, `denoise-ccs`) were invoked without `--p-n-threads`, so they ran on QIIME2's default of a single thread — the serial bottleneck of every ASV run. They now pass `--p-n-threads` with the per-dataset thread budget. Affects both download and `--local` ASV runs; ASV outputs are unchanged (DADA2 is deterministic across thread counts), only faster.

## Installation

Meta2Data can be installed in a local folder to avoid contaminating your QIIME2 environment and to make updates easier. QIIME2 is only required if you plan to run AmpliconPIP or AmpliconTAXA; MetaDL runs on any Python 3 interpreter.

### Step 1: Install QIIME2 and associated software (only for AmpliconPIP / AmpliconTAXA)

Create the conda environment from the provided `env.yml` (replace `<env-name>` with a name of your choice):

```bash
conda env create -n <env-name> -f env.yml
conda activate <env-name>
```

### Step 2: Clone the repo and add it to PATH

```bash
git clone https://github.com/LinyangSun/Meta2Data.git
echo 'export PATH="'"$PWD"'/Meta2Data/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```
> If you also use conda, the order of lines in your rc file matters.
> Place the export PATH=$HOME/Meta2Data/bin:$PATH line before any conda activate command.


### Step 3: Verify installation

```bash
Meta2Data --help
```

## Requirements

### System Requirements
- **OS**: Linux (tested on Ubuntu/CentOS)
- **Computation Resources**: For AmpliconPIP, 10 GB RAM, 4 CPU are recommended per parallel task. For AmpliconTAXA, 40-60GB RAM and 10CPU are recommonded



## Usage

Meta2Data provides several subcommands for different stages of the workflow:

```bash
Meta2Data <command> [options]

Available commands:
    MetaDL         Search keywords combination in INSDC and CNCB. Download and preclean metadata.
    AmpliconPIP    Download and process amplicon sequencing data based on user provided metadata.
    AmpliconTAXA        Merge amplicon datasets (--dada2 | --vsearch) and assign taxonomy using GreenGenes2 or SILVA.
    ShortreadsPIP  (In development)
```

---

### MetaDL: Metadata Download

Download metadata from INSDC and CNCB databases with parallel processing and checkpoint/resume capability. Automatically fetches BioProject descriptions and standardizes column names (CamelCase normalization, synonym merging via dictionary).

**Two modes:**

| Mode | Required Options | Description |
|------|-----------------|-------------|
| ID Input | `-i`, `-o` | Provide a directory of txt files containing BioProject IDs (PRJ*) and/or BioSample IDs (SAM*) — they can be mixed |
| Keyword Search | `-o`, `--keywords`, `--field`, `--organism` | Search INSDC + CNCB by keywords, then download metadata for matched BioProjects |

```
Required:
    -o, --output DIR              Output directory

ID Input Mode:
    -i, --input DIR               Directory with provided ID txt files (Only for BioProject and/or BioSample)

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

Main outputs (in the `-o` output directory):
- `all_metadata_merged.csv` — Final merged dataset with auto-fetched BioProject descriptions and SRA experiment design descriptions. Column names are standardized (CamelCase normalization + synonym merging via dictionary).
- `status.tsv` — Processing status for each input ID (`has_data` / `no_data` / `no_run_info` / `download_error`)
- `column_description.tsv` — Per-column statistics: fill rate, number of datasets covered, top 5 values, and column type (`core` / `cncb` / `rare`)
- `bioproject_absdesc.tsv` — One row per BioProject with publication info (`PMID`, `PMC`, `DOI`, `ArticleTitle`, `ArticleAbstract`, `PubSource`). Always generated; failures are non-fatal.
- `RecordWithoutRUNinfo.csv` — Records without SRA Run info (only generated when such records exist)

Keyword search mode only (`searched_keywords/`):
- `searched_keywords/combined_results.csv` — BioProjects matched by keyword search
- `searched_keywords/bioproject_ids.txt` — Matched BioProject accessions, one per line
- `searched_keywords/search_summary.txt` — Query summary

Internal / resume (`tmp/`, safe to delete after the run):
- `tmp/checkpoints/download_state.json` — Resume state for interrupted runs
- `tmp/<BioProject>.processed.csv` — Per-BioProject processed result (doubles as the resume checkpoint)
- `tmp/<group>_biosample.txt`, `tmp/<group>_sra_runinfo.csv` — Raw BioSample / SRA RunInfo fetched from the databases
- `tmp/<accession>.temp.csv`, `tmp/BIOSAMPLE_INPUT.processed.csv`, `tmp/SRA_INPUT.processed.csv` — scratch / direct SAM*/SRR*-input intermediates

> **Tip 1 — AI-assisted metadata screening**
>
> MetaDL automatically fetches a `Description` column for every BioProject, which summarizes each study's purpose, target organism, and experimental design. This makes the merged CSV well-suited for AI-based screening. Upload the CSV to any AI tool with a Team/collaborative workspace (Claude Team, ChatGPT Team, Gemini, etc.) and ask it to filter based on your criteria. For example:
>
> *"Here is my metadata CSV. Based on the Description and other columns, keep only gut microbiome samples from healthy human adults sequenced on Illumina with 16S amplicons. I do not want include any datasets that have illness. You need to label each item as include exclude and NotSure. You need to assign 3 workers to screen the full datasets independently, and a leader to give a final decision. Then output with xxx xxx xxx files. Ask me anything unclear to you before starting."*

> **Tip 2 — Broader keyword search with genus-level terms**
>
> Keyword search results depend heavily on how authors annotate their BioProjects. To maximize coverage for a taxonomic group of interest, don't rely solely on high-level terms (e.g., "bee"). Instead, collect genus names from a published phylogeny or species tree for your clade (for example from tree of life), and include them as `--organism` terms. For example:
>
> ```bash
> # Instead of just "bee", also search by genus names from the Apoidea phylogeny
> organism=("bee" "Apis" "Bombus" "Megachile" "Osmia" "Andrena" "Halictus")
> methods=("16S rRNA" "amplicon")
> Meta2Data MetaDL \
>     -o metadata/ \
>     --keywords \
>     --field "${methods[@]}" \
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
    --dada2 | --vsearch                 Denoising mode (REQUIRED: specify exactly one; no default)
                                  --dada2 : DADA2, single-base ASV resolution.
                                          Platforms Illumina / Ion Torrent / PacBio CCS;
                                          454, ONT and degraded/binned quality are skipped.
                                  --vsearch : vsearch, 97%-identity OTU.
                                          All 5 platforms; robust on mixed / degraded data.

Optional:
    -o, --output DIR              Output directory
                                  Default: current directory (--test) / metadata dir (normal)
    -t, --threads INT             Total CPU threads (default: 12)
                                  Auto-split: per-dataset threads = threads / max-parallel
    --max-parallel INT            Datasets to process in parallel (default: 2)
    --test                        Run in test mode
                                  Without -m: use built-in test data
                                  With -m: subset metadata (2 SRA per BioProject)
    -h, --help                    Show help

Local mode (process existing FASTQ instead of downloading):
    --skip-unknown-primers        Skip the dataset when an unmatched primer has fold <16.
                                  Default: trim the first 20 bp and continue.
    --parameter FILE             JSON overrides (docs/parameters.default.json).
    --local                       Read FASTQ straight from a folder; no download,
                                  no platform detection. -m is the INPUT FOLDER
                                  (one dataset = that folder; every FASTQ directly
                                  inside it is a sample; id = folder name).
    --platform PLATFORM           REQUIRED with --local. One of:
                                  ILLUMINA | LS454 | ION_TORRENT | PACBIO_SMRT | OXFORD_NANOPORE
                                  (applies to the whole run — one run = one platform).
    --primer-fwd SEQ              Optional. Forward primer to trim with cutadapt.
    --primer-rev SEQ              Optional. Reverse primer (paired-end). If no primer
                                  is given, the entropy auto-detector is used (as in
                                  download mode). --col-* are not needed in --local.
```

> **Note:** A denoising mode (`--dada2` or `--vsearch`) is mandatory and the two are mutually exclusive — there is no default. This replaces the previous auto-mix-by-platform behaviour, which silently combined ASV and OTU results within a single run (a hidden batch-effect risk). AmpliconTAXA automatically identifies the method when its input contains only one method. If both are present, select one explicitly.

> **Local mode:** `--local` skips download + automatic platform detection, so you must pass `--platform`. One `--local` run handles exactly **one dataset / one platform** (no sub-folder recursion) — for mixed-platform data, run each folder separately. Your original FASTQ files are never modified (they are symlinked read-only into the working directory). Example:
> ```bash
> Meta2Data AmpliconPIP --local --platform ILLUMINA --vsearch \
>     -m /path/to/my_fastq_folder -o /path/to/output -t 8
> # with explicit primers (cutadapt):
> Meta2Data AmpliconPIP --local --platform ILLUMINA --vsearch \
>     --primer-fwd GTGYCAGCMGCCGCGGTAA --primer-rev GGACTACNVGGGTWTCTAAT \
>     -m /path/to/my_fastq_folder -o /path/to/output -t 8
> ```
> ```bash
> Meta2Data AmpliconPIP --local --platform ILLUMINA --dada2 \
>     -m /path/to/my_fastq_folder -o /path/to/output -t 8
> # with explicit primers (cutadapt):
> Meta2Data AmpliconPIP --local --platform ILLUMINA --dada2 \
>     --primer-fwd GTGYCAGCMGCCGCGGTAA --primer-rev GGACTACNVGGGTWTCTAAT \
>     -m /path/to/my_fastq_folder -o /path/to/output -t 8
> ```

**Metadata CSV format** (column names customizable via `--col-*`):
```csv
Bioproject,Run
PRJNA12345,SRR123456
PRJNA12345,SRR123457
PRJNA67890,SRR234567
```

**Processing pipeline:**
1. Download SRA data (via FTP) — or, with `--local`, read FASTQ straight from a folder (no download)
2. Detect sequencing platform (Illumina, PacBio, Ion Torrent, 454, ONT) and layout (single/paired-end) automatically from INSDC/CNCB — or use `--platform` in local mode
3. Quality control
4. Identify and trim primers (entropy auto-detection, or explicit primers in local mode)
5. Mode-specific denoising:
   - `--dada2`: DADA2 (Illumina / Ion Torrent / PacBio CCS); 454, ONT and degraded/binned-quality data are skipped
   - `--vsearch`: vsearch 97% OTU clustering (all 5 platforms, robust on degraded/binned-quality and ONT data)
6. Generate QIIME2 artifacts (`.qza` files)
7. Per-dataset summary: write `per_dataset_summary.tsv` (platform + quality status + amplified 16S region)

**Output structure** (`<mode>` = `dada2` | `vsearch`):
```
<output_dir>/
├── datasets_ID.txt                            # Generated dataset list
├── <dataset_ID>/                              # One directory per dataset
│   ├── <dataset_ID>_sra.txt                  # SRA accession list
│   ├── ori_fastq/                             # Downloaded FASTQ files
│   ├── <dataset_ID>-<mode>-final-rep-seqs.qza
│   └── <dataset_ID>-<mode>-final-table.qza
├── logs/                                      # Per-dataset verbose logs (<dataset_ID>.log)
├── datasets.log                               # Unified status log (append-only): SUCCESS|FAILED|SKIPPED|LOW_QUALITY, one dated header per run
├── summary.csv                                # Ordered per-sample stage counts; DADA2/VSEARCH prefixes; keyed by dataset/method/sample
├── pip_dataset_read_counts.csv                # Per-BioProject stage totals (including pooled VSEARCH abundances)
└── per_dataset_summary.tsv                    # Per-dataset: platform + quality status + amplified 16S region (V-region via E. coli alignment); upserted by Bioproject
```

The amplified region in `per_dataset_summary.tsv` is inferred by aligning each dataset's
representative sequences to the E. coli 16S reference (`docs/ecoli_16S_J01859.fasta`) and
mapping the median E. coli span to the V1–V9 regions (e.g. `V3-V4`, `V1-V9` for full-length),
with a confidence = fraction of rep-seqs agreeing.

---

### AmpliconTAXA: Merge & Taxonomy Assignment

```text
Required:
    --db DIR                     Database directory
    -i, --input DIR              Search this directory recursively for PIP final results

Optional:
    --dada2 | --vsearch           Select a method; inferred when only one is present
    --singleV                    Single-region alignment and de novo tree
                                 Default: multi-region SEPP reference-tree insertion
    --notree, --no-tree           Skip tree construction and tree-dependent filtering
                                 Mutually exclusive with --singleV
    --parameter FILE             JSON parameter overrides
    --db-type TYPE               greengenes (default) or silva
    --confidence FLOAT           Classification confidence (default: 0.7)
    -o, --output DIR             Output directory (default: input directory)
    -t, --threads INT            Threads (default: 4)
    --dl                         Download missing database files
```

Dataset directory names are unrestricted. TAXA finds matching
`<id>-<method>-final-table.qza` and `<id>-<method>-final-rep-seqs.qza`
in the same directory, validates the artifact types, and removes copied or linked
duplicate artifact pairs. Incomplete/invalid pairs are recorded and skipped;
conflicting dataset names or artifacts are reported as errors. Temporary
directories and TAXA aggregate directories marked by `taxa-run-state.json` are
excluded. Local dataset names beginning with `final-` remain valid. Both methods
in the same input require an explicit `--dada2` or `--vsearch`; they are never
silently combined. Sample IDs must be unique across included datasets; overlapping
sample IDs cause a merge error. Rename local sample files before PIP processing
if separate datasets reuse the same sample names.

New vsearch results use sequence-derived SHA-256 feature IDs so independently
numbered centroids cannot collide across datasets. Identical sequences share IDs;
this does not change clustering parameters or sample counts.

All TAXA modes merge feature tables and sequences, orient sequences using
the GreenGenes2 backbone, filter the table to oriented features, and classify with
the chosen GreenGenes2 or SILVA classifier. `--singleV` then builds a de novo tree
with MAFFT/masking/FastTree. By default, SEPP inserts features into the
existing Greengenes 13_8 reference (`sepp-refs-gg-13-8.qza`) and filters tables and
sequences to placed features. `--notree` (alias `--no-tree`) skips tree construction
and tree-dependent filtering, and requires no SEPP reference, including with `--dl`.
The oriented table and representative sequences are its final outputs.
`--singleV` and `--notree` are mutually exclusive, including when selected through
the parameter file. The feature method does not choose the tree workflow.
Different region sequences remain distinct features; taxonomy classification does
not automatically collapse them into taxon-level abundance counts.

Outputs are isolated by method and region workflow:

```text
final-<dada2|vsearch>-<singleV|multiV|notree>/
    collection.json                 # included, duplicate and skipped result pairs
    parameters.json                 # resolved settings
    taxa-run-state.json             # cache dependencies
    <gg2|silva>Taxonomy.qza
    # singleV:
    orientedTable.qza, orientedRepSeqs.qza, denovoRootedTree.qza
    # multiV:
    treeFilteredTable.qza, treeFilteredRepSeqs.qza, seppTree.qza
    # notree:
    orientedTable.qza, orientedRepSeqs.qza
    tmp/                            # intermediate artifacts and diagnostics
```

Input, reference and classification-parameter changes invalidate their dependent
artifacts. A database change reuses the merge and orientation steps when their
inputs are unchanged. GG2 and SILVA classifications are cached independently,
so switching classifiers preserves each valid result. Required tree insertion
and feature filtering failures return a nonzero exit status and retain
intermediate artifacts for diagnosis.

### Custom parameters and primer decisions

PIP and TAXA accept `--parameter my_parameters.json`. The complete template is
[docs/parameters.default.json](docs/parameters.default.json); supply only the
fields you want to override. Explicit CLI flags take precedence over the JSON
file, which takes precedence over built-in defaults. Unknown fields, incorrect
types and out-of-range values fail before data processing.

```json
{
  "primer": {"fold_threshold": 16, "skip_unknown": true},
  "vsearch": {"maxee": 1.0, "cluster_identity": 0.97},
  "taxa": {"confidence": 0.8, "singleV": false, "notree": false}
}
```

The b1 detector uses the first 20 bases of all quality-filtered reads from the
selected sample (forward and reverse analyzed separately). Bases supported by
at least 10% of reads define C/D/T/Q states, contributing fold factors 1/2/3/4.
Database matching takes precedence, with identity >=0.85 and informative-position
fraction >=0.50. The bundled benchmark database contains 10 forward and 13 reverse
primers. Only unmatched sequences use the strict fold rule:

| Decision | Default behavior | With `--skip-unknown-primers` |
|---|---|---|
| Database match | Trim to the matched endpoint | Same |
| No match and fold <16 | Trim 20 bp | Skip the whole dataset |
| No match and fold >=16 | Leave reads unchanged | Same |
| No valid detection reads | Record failure | Same |

`<dataset>-<method>-primer_info.json` records layout, decisions, thresholds and
actual trimming outside temporary directories. In mixed-orientation PE data,
actual R1 and R2 reads are evaluated separately within each orientation group.
Unknown primers on any required end/orientation trigger the skip policy.
Explicit cutadapt runs report variable trim lengths as `null`; detection-only
runs record zero trimming. Additional fixed 5-prime trimming after primer removal
is zero by default for degraded-quality and Ion workflows; the relevant JSON
settings can request additional trimming explicitly.

Both PacBio workflows retain the full-length eligibility check: more than half
of the first 1,000 adapter-removed reads must exceed 1,400 bp. The configurable
PacBio minimum/maximum lengths apply to subsequent read filtering and do not
change this platform eligibility check.

PacBio vsearch applies automatic primer trimming, including mixed orientations.
DADA2 CCS detects primers without pre-trimming because `denoise-ccs` requires a
known or explicitly supplied forward primer to orient and trim reads itself.
If that prerequisite cannot be met, the dataset is skipped with an explanation;
use explicit primers or the vsearch workflow for that dataset.

---

## Examples

### Case 1: Different way to run metaDL

```bash
conda activate <env-name>

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
conda activate <env-name>

# Custom column names matching your CSV headers
Meta2Data AmpliconPIP \
    -m my_samples.csv \ # your metadata.csv
    --col-bioproject "ProjectID" \ # column name for bioproject
    --col-sra "SRA_Accession" \ # column name for sra run (the sra normally start with SRR, ERR, DRR OR CRR)
    --vsearch \
    -o amplicon_output/ \
    -t 8
```

### Case 3: Taxonomy Assignment Only (AmpliconPIP Already Complete)

Run AmpliconTAXA independently on existing AmpliconPIP results, e.g., to compare databases.

```bash
conda activate <env-name>

# vsearch results, multi-region workflow, GreenGenes2 classification
Meta2Data AmpliconTAXA --vsearch \
    --db path/to/your/metafile/databases/ \ # Prepare an empty path, the pip will download database in it
    --db-type greengenes \
    --dl \
    --confidence 0.7 \
    -i path/to/your/metafile/amplicon_output/ \
    -t 16

# DADA2 results from one region, SILVA 138.99 classification
Meta2Data AmpliconTAXA --dada2 --singleV \
    --db path/to/your/metafile/databases/ \
    --db-type silva \
    --dl \
    --confidence 0.7 \
    -i path/to/your/metafile/amplicon_output/ \
    -t 16

# Compare databases on the SAME mode: the second run reuses the merge/orient/tree
# steps and only recomputes the taxonomy artifact.
Meta2Data AmpliconTAXA --vsearch --db-type greengenes --db databases/ -i amplicon_output/ -t 16
Meta2Data AmpliconTAXA --vsearch --db-type silva      --db databases/ -i amplicon_output/ -t 16
```



### Case 4: Test Mode — Quick Validation

Verify the pipeline works before running on your full dataset. A denoising
mode (`--dada2` / `--vsearch`) is required even in test mode.

`--test` has two forms:
- **without `-m`** — runs the built-in test data (fastest way to verify the install).
- **with `-m`** — subsets *your* metadata to 2 SRA runs per BioProject and runs the
  real pipeline on that small slice, so you can confirm your CSV format, `--col-*`
  names, and accessions all work before launching the full run.

```bash

# Use built-in test data (fastest way to verify installation)
Meta2Data AmpliconPIP --test --vsearch -t 8

# Or smoke test on your own data: it confirms your CSV format, your --col-* names, and that your specific accessions download and process — but in a few minutes on 2 runs/project instead of hours/days on the full set
Meta2Data AmpliconPIP --test \
    -m path/to/your/metafile/metadata.csv \
    --col-bioproject Bioproject \
    --col-sra Run \
    --vsearch \
    -o test_output/ \
    -t 8
```



### Case 5: Process Local FASTQ (No Download)

Already have the FASTQ files? Use `--local` to process a folder directly — no download, no NCBI lookup. The folder **is** one dataset (its name becomes the dataset id) and every FASTQ file directly inside it is a sample or a direction of sample. You must supply `--platform` and an explicit output directory (`-o`); your original files are never modified (they are symlinked read-only).

```bash

# Auto-detect primers (b1), vsearch method, Illumina data
Meta2Data AmpliconPIP --local --platform ILLUMINA --vsearch \
    -m path/to/my_fastq_folder/ \   # input folder = one dataset (id = folder name)
    -o local_output/ \
    -t 8

# Provide explicit primers (trimmed with cutadapt) instead of auto-detection
Meta2Data AmpliconPIP --local --platform ILLUMINA --vsearch \
    --primer-fwd GTGYCAGCMGCCGCGGTAA \
    --primer-rev GGACTACNVGGGTWTCTAAT \
    -m path/to/my_fastq_folder/ \
    -o local_output/ \
    -t 8
```

> - `--platform` must be one of `ILLUMINA | LS454 | ION_TORRENT | PACBIO_SMRT | OXFORD_NANOPORE`.
> - One `--local` run = one dataset / one platform; for mixed-platform data, run each folder separately.
> - `--max-parallel` is forced to `1` in `--local` mode (a single dataset has nothing to parallelize across), so the single dataset always gets all `-t` threads. To process several folders concurrently, launch one `--local` run per folder.
> - Paired-end is detected from `_1`/`_2` or `_R1`/`_R2` filename suffixes; `.fq`/`.fq.gz` are accepted (normalized to `.fastq`).
> - The same outputs as download mode are produced (`datasets.log`, `summary.csv`, `per_dataset_summary.tsv`).

Local paired filenames may use `_R1/_R2`, `_1/_2`, `_R1_001/_R2_001`,
or `_1_001/_2_001`, with `.fastq`, `.fq`, or gzip-compressed extensions. For example:

```text
PRM190709-001_S205_L001_1_001.fastq.gz
PRM190709-001_S205_L001_2_001.fastq.gz
```

The direction is `1/2`; the last `001` is a chunk number. Matching lane/chunk pairs
are combined in a stable order under the same sample ID. Original files are left
unchanged, and `read_layout.json` records their mapping. Unpaired or ambiguous
names produce a filename error with a renaming example; there is no user manifest
option. Plain single-end reads should use names such as `sample.fastq.gz`.
Mixed PE/SE datasets should be separated into different folders.

One local folder remains one dataset, named after the folder. To process a subset
of downloaded datasets, supply a metadata CSV containing just those datasets.
For TAXA, choose an input directory containing the desired results. No dataset-ID,
manifest, or dataset-selection flags are provided.

For validation coverage and execution limits, see
[docs/revision-validation.md](docs/revision-validation.md).

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


## Workflow

<p align="center">
  <img src="docs/meta2data workflow.png" alt="Meta2Data logo" width="100%">
</p>




### Per-sample protection against false inferred adapters

AmpliconPIP accepts `--adapter-guard --db /path/to/gg2 --dl` (or
`--adapter-guard --adapter-ref reference.qza`). This optional guard checks
fastp's inferred adapter sequences against the 16S reference, then reruns only
samples with a strong biological match from their original FASTQs. Existing
FASTQ names and `summary.csv` stay compatible. See [adapter guard](docs/adapter_guard.md)
for thresholds, provenance, fallback behavior and limitations.
