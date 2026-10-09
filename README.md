<p align="center">
  <img src="docs/llloooooog.png" alt="Meta2Data logo" width="100%">
</p>

# Meta2Data

Meta2Data is a Linux toolkit for amplicon analysis, combining public metadata retrieval, sequencing data processing, taxonomic classification, and phylogenetic tree construction. It also supports local FASTQ datasets.

## Workflow

```text
MetaDL: retrieve and organize metadata
   ↓ Select projects and samples for analysis
AmpliconPIP: download/read FASTQ → adapter and primer processing → quality control → dada2 or vsearch
   ↓ Feature tables + representative sequences
AmpliconTAXA: collect results recursively → merge and orient sequences → classify → optionally build a tree
```

| Module | Input | Main output |
|---|---|---|
| `MetaDL` | Project/sample/Run IDs, or keywords | Merged metadata CSV |
| `AmpliconPIP` | Public metadata CSV, local dataset CSV, or both | Feature table, representative sequences, and processing statistics for each dataset |
| `AmpliconTAXA` | A common parent directory containing one or more PIP results | Merged feature table, representative sequences, taxonomy, and an optional phylogenetic tree |

Select a processing method explicitly in PIP: `--dada2` supports Illumina, Ion Torrent, and full-length PacBio CCS; `--vsearch` also supports 454, Oxford Nanopore, and data with degraded or binned quality scores. In this pipeline, dada2 skips unsupported data rather than switching methods automatically.

## Installation

Choose either Conda or SIF.

### Conda

```bash
export M2D_ROOT=/absolute/path/to/Meta2Data  # Source directory
conda env create -n meta2data -f "$M2D_ROOT/env.yml"
conda activate meta2data
export PATH="$M2D_ROOT/bin:$PATH"
Meta2Data --help
```

`env.yml` includes QIIME 2 2024.10 and the pipeline dependencies.

### SIF / Apptainer

Use a SIF version that provides the interface documented here. Meta2Data and its dependencies run directly from the image. After installing Apptainer, choose a version from [GitHub Packages](https://github.com/LinyangSun/Meta2Data/pkgs/container/meta2data) and replace `TAG` below with its tag:

```bash
apptainer pull Meta2Data.sif \
  oras://ghcr.io/linyangsun/meta2data:TAG

apptainer exec --cleanenv Meta2Data.sif Meta2Data --help
```

For convenience, define a shortcut:

```bash
export M2D_SIF="$(realpath Meta2Data.sif)"

m2d() {
  apptainer exec --cleanenv "$M2D_SIF" Meta2Data "$@"
}

m2d --help
```

Replace `Meta2Data` with `m2d` in the examples below. Run commands from your working directory; PIP and TAXA write to its `results/` directory. Set the image path and define the `m2d` function again in each new terminal.

See the [container guide](docs/container.md) for more examples.

## Preparing your analysis

Run the examples below from the same working directory. PIP and TAXA share fixed input/output locations, while MetaDL uses a separately specified output directory:

```text
working_directory/
├── id_files/                # Input directory for ID-file mode; location is configurable
│   ├── projects.txt
│   └── metadata/            # MetaDL output for ID-file mode
├── metadata/
│   └── bee/                 # Example MetaDL output for keyword mode
├── metadata.csv             # Selected public projects and Runs for PIP
├── local_metadata.csv       # Local datasets, paths, and platforms
├── local_data/
│   ├── project_A/
│   └── project_B/
└── results/
    ├── pip/
    │   ├── PRJNA1402755/     # Public project
    │   ├── project_A/        # Local dataset
    │   └── project_B/
    ├── db/                  # Shared reference resources
    └── taxa/                # Merged tables, taxonomy, and trees
```

Specify MetaDL output with `-o`. The examples place ID-file output in a new subdirectory of the input directory, and keyword output in `metadata/<task_name>/`. After filtering the metadata, save the projects and Runs selected for analysis as `metadata.csv` for PIP.

All PIP modes—public batch, selected projects, local, mixed, and test—write to `results/pip/` under the directory where the command starts. TAXA collects results recursively from there by default, or from an existing PIP directory specified with `-i`; its output is always `results/taxa/` under the launch directory. Reference resources are cached in `results/db/`. Use a separate working directory for each analysis.

Use the same processing method for datasets merged in one TAXA run, and keep sample IDs unique across datasets. Run independent PIP commands sequentially when they share an output root; a second command reports that the directory is in use if another PIP run is active. After rerunning a project, retain only the intended version for merging so that TAXA does not collect both old and new results.

## 1. MetaDL: retrieve and organize metadata

MetaDL supports INSDC and CNCB/GSA and retrieves project, sample, and Run information. Filter the metadata according to your research question before analysis.

Choose one of two independent input modes: **keyword mode** searches by terms; **ID-file mode** reads BioProject, BioSample, or SRA accessions supplied by the user. ID-file mode does not require a previous keyword search.

### Keyword mode

Enable `--keywords` and supply search terms directly.

| Parameter | Description |
|---|---|
| `--keywords` | Enable keyword mode |
| `--field TERM...` | Required; method-related or similar search terms |
| `--organism TERM...` | Required; terms describing the organisms of interest |
| `--opt TERM...` | Optional additional search terms |

For example, use Bash arrays for bee-related terms and methods:

```bash
bee=("Apis" "Bombus" "Megachile" "Osmia" "Andrena" "Halictus")
methods=("16S rRNA" "amplicon")

Meta2Data MetaDL \
  --keywords \
  --field "${methods[@]}" \
  --organism "${bee[@]}" \
  -o metadata/bee
```

### ID-file mode

**Without `--keywords`, specify an input directory using `-i, --input DIR`.** Put the BioProject or BioSample IDs to retrieve in one or more `.txt` files, one ID per line without a header, and place these files in that directory. SRA accessions are also accepted, including Run IDs starting with `SRR`, `ERR`, or `DRR`.

For example, prepare a directory containing:

```text
id_files/
├── projects.txt      # One BioProject ID per line
└── samples.txt       # One BioSample ID per line
```

Either file type alone is sufficient; directory and file names are configurable. Example contents of `projects.txt`:

```text
PRJNA1402755
PRJNA863317
```

Pass the directory containing the files:

```bash
Meta2Data MetaDL \
  -i /path/to/id_files \
  -o /path/to/id_files/metadata
```

The program creates the `metadata/` output directory under `/path/to/id_files/`:

```text
id_files/
├── projects.txt
├── samples.txt
└── metadata/
    ├── all_metadata_merged.csv
    ├── status.tsv
    └── tmp/
```

`-i` takes a **directory path, not an individual `.txt` file**. The directory can be anywhere accessible and is independent of keyword-mode output. MetaDL reads `.txt` files directly within it, without scanning subdirectories.

### Parameters shared by both modes

| Parameter | Description |
|---|---|
| `-o, --output DIR` | Required metadata output directory; the program creates the subdirectories shown in the examples |
| `-k, --api-key KEY` | Optional NCBI API key |
| `-w, --max-workers N` | Parallel workers; defaults to `3` without an API key or `8` with one |
| `-h, --help` | Show help |

### Main output

| File or directory | Contents |
|---|---|
| `all_metadata_merged.csv` | Merged metadata |
| `status.tsv` | Download status |
| `column_description.tsv` | Field statistics |
| `bioproject_absdesc.tsv` | Publication information associated with projects |
| `searched_keywords/` | Search results generated in keyword mode only |

Resume an interrupted task by rerunning it with the same output directory. Use separate output directories for different ID lists or keyword searches.

## 2. AmpliconPIP: sequencing data processing

### Automatic primer detection and trimming

Primers are detected automatically for each dataset unless supplied manually. This applies to standard public-data batch processing and to local datasets whose primer columns are not mapped or whose mapped primer fields are both empty. Manually supplied primers are still trimmed, but automatic detection and fallback are disabled.

The default rules are:

| Detection result | Action |
|---|---|
| Matches the known primer database | Trim through the end of the match |
| No match, and fold in the first 20 bp is strictly less than 16 | Treat as an unknown primer and trim the first 20 bp by default |
| No match, and fold is at least 16 | Treat this end as primer-free and do not trim it |

`--skip-unknown-primers`: **skip the entire dataset** if automatic detection finds an unknown primer. This option does not bypass primer processing and continue the analysis.

**PacBio dada2 has an additional requirement:** after passing the near-full-length check, a known forward primer must still be detected or explicitly supplied. The preliminary step detects primers without trimming them; `denoise-ccs` then handles orientation and primer removal. Datasets without a forward primer are skipped. Full-length reads are not necessarily primer-free, and the output retains near-full-length representative sequences. vsearch does not require detection of a known forward primer, but still follows the detection and trimming rules above.

The following settings can be configured in a JSON file supplied with `--parameter`:

| JSON parameter | Default | Description |
|---|---:|---|
| `primer.window` | `20` | Detection window at the start of a read, in bp |
| `primer.fold_threshold` | `16` | When no database primer matches, fold strictly below this threshold indicates an unknown primer |
| `primer.support_frequency` | `0.10` | Minimum frequency for a base to contribute support at each position; affects the consensus and fold |
| `primer.database_identity` | `0.85` | Minimum proportion of IUPAC-compatible non-N positions when matching the primer database |
| `primer.informative_fraction` | `0.50` | Minimum ratio of non-N positions in the matched segment to the full length of the database primer |
| `primer.unknown_trim_length` | `20` | Bases removed from the start when an unknown primer is found and the dataset is not skipped |
| `primer.skip_unknown` | `false` | Whether to skip the entire dataset when an unknown primer is detected; corresponds to `--skip-unknown-primers` |

The following parameters select **reads used for detection only**. They do not directly filter reads across the full dataset:

| JSON parameter | Default | Condition for inclusion in detection |
|---|---:|---|
| `primer.min_length` | `50` | Read length of at least 50 bp |
| `primer.min_average_quality` | `20` | Mean Phred score of at least 20 across the read; this filter is skipped when constant placeholder quality scores are detected |
| `primer.min_complexity` | `0.3` | Number of distinct 2-mers divided by 16, at least 0.3 |
| `primer.min_entropy` | `1.0` | Shannon entropy of the A/C/G/T composition of at least 1.0 |

Detection uses the first sample after sorting the dataset's samples by name, and the resulting decision applies to the entire dataset. If no valid detection reads remain, processing fails; this is not interpreted as “no primer.” Manually supplied primers do not use the automatic detection parameters above.
See the [JSON parameter guide](docs/parameters.md) for accepted ranges, platform applicability, and configuration examples.

### Adapter protection

Adapter protection is enabled by default. It checks whether an “adapter” inferred by fastp could actually be a biological 16S sequence. If a match is found, the sample is reprocessed from the original FASTQ: adapter trimming is disabled for single-end reads; for paired-end reads, automatic adapter inference is disabled while overlap-based adapter trimming remains enabled. `--no-adapter-guard` disables only this additional check. Standard adapter removal, primer trimming, and downstream quality control still run.

### Processing modes

| Mode | Input and purpose |
|---|---|
| Public and local data together | Supply both a public metadata CSV and a local dataset CSV to process all data in one run |
| Rerun selected public projects | Select projects from public metadata and supply known primers for each project |
| Public data only | Download data using public metadata and automatically identify the platform and primers |
| Local data only | Read FASTQ files listed in a local CSV, with a platform and optional primers specified for each row |

All modes write to `results/pip/` under the directory where the command is launched. Select either `--dada2` or `--vsearch`.

### Processing public and local data together

Prepare two CSV files. The public metadata file, `metadata.csv`, contains one Run per row:

```csv
Bioproject,Run
PRJNA1402755,SRR36832824
PRJNA863317,SRR20818414
```

The local metadata file, `local_metadata.csv`, contains one dataset per row:

```csv
datasets,path,platform
project_A,local_data/project_A,ILLUMINA
project_B,local_data/project_B,ION_TORRENT
```

Each local path points to a directory that directly contains FASTQ files. Relative paths are resolved against **the directory containing the local CSV**. Dataset names determine the output directory names. `datasets`, `path`, and `platform` are example column names; all three column-mapping parameters must be supplied explicitly.

```bash
Meta2Data AmpliconPIP \
  --public-m metadata.csv \
  --public-bioproject-colNAME Bioproject --public-sra-colNAME Run \
  --local-m local_metadata.csv \
  --local-datasets-colNAME datasets \
  --local-path-colNAME path --local-platform-colNAME platform \
  --dada2 -t 8
```

Public and local results are written to `results/pip/<project-or-dataset-name>/`. To read local primer columns, also supply `--local-primer-f-colNAME` and, optionally, `--local-primer-r-colNAME`; see the local examples below. Mixed input cannot be combined with `--public-bioprojectIDs` or `--test`.

### Main parameters

**Public data input**

| Parameter | Description |
|---|---|
| `--public-m FILE` | Required for public or mixed input; public metadata CSV |
| `--public-bioproject-colNAME NAME` | Required; project ID column in the public CSV |
| `--public-sra-colNAME NAME` | Required; Run ID column in the public CSV |

**Local data input**

| Parameter | Description |
|---|---|
| `--local-m FILE` | Required for local or mixed input; local dataset CSV |
| `--local-datasets-colNAME NAME` | Required; dataset name column, with no default |
| `--local-path-colNAME NAME` | Required; FASTQ directory column, with no default |
| `--local-platform-colNAME NAME` | Required; platform column, with no default |
| `--local-primer-f-colNAME NAME` | Optional; forward primer column |
| `--local-primer-r-colNAME NAME` | Optional; reverse primer column; requires the forward primer column to be specified |

**General parameters**

| Parameter | Description / default |
|---|---|
| `--dada2` / `--vsearch` | Choose exactly one; no default |
| `-t, --threads N` | Total threads; default `4` |
| `--max-parallel N` | Concurrent datasets; default `2`, also applies to multiple local datasets |
| `--skip-unknown-primers` | Skip the entire dataset if an unknown primer is detected; see the start of this section |
| `--parameter FILE` | JSON parameter overrides; see the advanced parameters section |
| `--adapter-guard` | Explicitly enable adapter protection; already enabled by default |
| `--no-adapter-guard` | Disable the additional adapter protection check |
| `--test` | Public-data test list; see the test dataset section |
| `-h, --help` | Show help |

Threads are divided evenly among concurrent datasets. For example, `-t 8 --max-parallel 2` assigns 4 threads per dataset. A single dataset uses all requested threads.

### Rerunning selected public projects with known primers

`--public-bioprojectIDs` is disabled by default. If automatic processing produces unexpected results, select all Runs belonging to specific projects from the original public CSV and rerun them with known primers for each project. This mode cannot be combined with local input.

| Parameter | Description |
|---|---|
| `--public-bioprojectIDs ID...` | Public projects to process, separated by spaces |
| `--public-primer-fwd SEQ...` | Required; one sequence per project, in the same order |
| `--public-primer-rev SEQ...` | May be omitted entirely; if supplied, the number and order must also match the projects |

Public primer parameters take **lists of nucleotide sequences** directly; local primer parameters take **CSV column names**. The following project IDs and primers only illustrate the pairing. Replace them using your public metadata and original experimental records:

```bash
projects=("PRJNA123456" "PRJNA654321")
forward=("GTGYCAGCMGCCGCGGTAA" "CCTACGGGNGGCWGCAG")
reverse=("GGACTACNVGGGTWTCTAAT" "GACTACHVGGGTATCTAATCC")
Meta2Data AmpliconPIP \
  --public-m metadata.csv \
  --public-bioproject-colNAME Bioproject --public-sra-colNAME Run \
  --public-bioprojectIDs "${projects[@]}" \
  --public-primer-fwd "${forward[@]}" --public-primer-rev "${reverse[@]}" \
  --dada2 -t 8
```

Each selected project uses the supplied primers, with no automatic detection or fallback. Standard public batch processing does not accept manually supplied primers. The selected metadata is saved as `selected_metadata.csv`, and the project-to-primer mapping as `project_primers.json`. The original CSV is not modified.

### Processing public data only

Save the filtered public metadata as `metadata.csv` using the format above, and supply only the public input parameters:

```bash
Meta2Data AmpliconPIP \
  --public-m metadata.csv \
  --public-bioproject-colNAME Bioproject --public-sra-colNAME Run \
  --dada2 -t 8
```

All projects in the CSV are processed, with automatic platform and primer detection. Additional platform or quality-description columns in the public CSV do not directly control processing.

### Processing local data only

Each row in the local CSV corresponds to a directory that directly contains FASTQ files, for example:

```text
local_data/
├── project_A/
│   ├── A01_R1.fastq.gz
│   └── A01_R2.fastq.gz
└── project_B/
    └── B01.fastq.gz
```

Use `local_metadata.csv` from above to detect primers automatically for each dataset:

```bash
Meta2Data AmpliconPIP --local-m local_metadata.csv \
  --local-datasets-colNAME datasets \
  --local-path-colNAME path --local-platform-colNAME platform \
  --dada2 -t 8
```

Results are written to `results/pip/project_A/` and `results/pip/project_B/`.

**Primer columns are read only when explicitly specified.** For example, add two columns to the local CSV:

```csv
datasets,path,platform,primer_f,primer_r
project_A,local_data/project_A,ILLUMINA,GTGYCAGCMGCCGCGGTAA,GGACTACNVGGGTWTCTAAT
project_B,local_data/project_B,ION_TORRENT,,
```

```bash
Meta2Data AmpliconPIP --local-m local_metadata.csv \
  --local-datasets-colNAME datasets \
  --local-path-colNAME path --local-platform-colNAME platform \
  --local-primer-f-colNAME primer_f --local-primer-r-colNAME primer_r \
  --dada2 -t 8
```

If both mapped primer fields are empty, primers are detected automatically. If a forward primer is supplied, the sequences in that row are used, with no fallback to automatic detection. The reverse primer may be empty, but cannot be supplied on its own. Known primers at both ends of long single-end reads are trimmed; shorter reads that do not reach the reverse primer are retained. Replace the example sequences using the original experimental records.

Different local datasets may use different platforms: `ILLUMINA`, `LS454`, `ION_TORRENT`, `PACBIO_SMRT`, or `OXFORD_NANOPORE`. A single dataset cannot mix platforms or single-end and paired-end reads.

Supported extensions are `.fastq`, `.fq`, and their `.gz` equivalents. Paired-end files support the following naming patterns; extensions are omitted in the table:

| Forward filename | Reverse filename |
|---|---|
| `sample_R1` | `sample_R2` |
| `sample_1` | `sample_2` |
| `sample_R1_001` | `sample_R2_001` |
| `sample_1_001` | `sample_2_001` |
| `sample_L001_R1_001` | `sample_L001_R2_001` |

Local single-end files should use names without a pairing suffix, such as `sample.fastq.gz`. An unpaired `sample_R1.fastq.gz` or `sample_1.fastq.gz` causes a missing-mate error. Single-end files named `sample_L001.fastq` and `sample_L002.fastq` are treated as two separate samples and are not merged automatically.

For paired-end data, multiple lanes or chunks from the same sample are sorted by name, then merged separately for R1 and R2. Every group must contain both mates. The paired-end sample ID is the filename prefix after removing the extension and read, lane, and chunk suffixes. For example, `A_S1_L001_R1_001.fastq.gz` becomes `A_S1`. Single-end sample IDs are filenames with the extension removed.

Dataset names must be unique and must not conflict with public project IDs in the same run. The same FASTQ directory cannot be registered twice, and datasets with identical names but different source directories are not overwritten. Duplicate sample IDs within a local or mixed run cause a preflight error; they are not renamed automatically. Sample names must also remain unique across separate runs, and TAXA checks them during merging. Input sources and output directories must not overlap. Keep metadata CSV files outside the results directory. Original FASTQ files are not modified.

### Main outputs

```text
results/pip/
├── <dataset>/
│   ├── <dataset>-<method>-final-table.qza
│   ├── <dataset>-<method>-final-rep-seqs.qza
│   └── read_counts/                 # Per-stage counts and quality-control reports
├── summary.csv                      # Per-sample read counts at each processing stage
├── pip_dataset_read_counts.csv      # Dataset-level summary
├── per_dataset_summary.tsv          # Platform, quality, amplicon region, etc.
├── effective-parameters-<method>.json
├── datasets.log                     # Success, failure, skipped, and other statuses
└── logs/<dataset>.log                # Detailed log for each dataset
```

`<method>` is `dada2` or `vsearch`. After successful processing, working FASTQ files and temporary files are removed; final results and statistics are retained. Completed steps can be reused when inputs and parameters are unchanged. Changes to inputs or parameters trigger reprocessing of the affected steps.

## 3. AmpliconTAXA: merge, classify, and build trees

By default, TAXA recursively collects complete public and local results from `results/pip/`, merges feature tables and representative sequences, orients sequences, assigns taxonomy, and optionally builds a tree. Output is always written to `results/taxa/` under the launch directory.

```bash
Meta2Data AmpliconTAXA --dada2 --classifier greengenes -t 8
```

| Parameter | Description / default |
|---|---|
| `-i, --input DIR` | Recursively search this directory; defaults to `results/pip` under the launch directory |
| `--dada2` / `--vsearch` | Mutually exclusive; must match PIP. Detected automatically if only one method is present; select explicitly when both are present |
| `--classifier greengenes / silva` | Classifier; defaults to `greengenes` and downloads automatically if missing |
| `--singleV` | For data from the same amplified region: align sequences and build a tree de novo |
| `--notree` / `--no-tree` | Merge, orient, and classify without building a tree |
| `--confidence FLOAT` | Classification confidence; default `0.7`, range `0–1` |
| `-t, --threads N` | Threads; default `4` |
| `--parameter FILE` | Override parameters using JSON |
| `-h, --help` | Show help |

All modes orient sequences against the GG2 sequence reference and remove features that cannot be oriented from both the representative sequences and feature table. The default SEPP workflow also removes features that cannot be inserted into the tree. `--singleV` and `--notree` are mutually exclusive, and neither requires the SEPP reference. `--notree` skips tree construction and tree-insertion filtering only. Classifier selection is independent of the tree method; resources are prepared automatically in `results/db/`.

If the recorded state identifies results from an input source that has since been replaced, TAXA reports an error and preserves the files. Rerun PIP with that method first. External results without a state record can still be imported.

To classify with SILVA without building a tree:

```bash
Meta2Data AmpliconTAXA --dada2 --classifier silva --notree -t 8
```

Results are saved in `results/taxa/final-<method>-<multiV|singleV|notree>/`:

| Workflow | Main output |
|---|---|
| Default, multiple regions | `treeFilteredTable.qza`, `treeFilteredRepSeqs.qza`, `seppTree.qza` |
| `--singleV` | `orientedTable.qza`, `orientedRepSeqs.qza`, `denovoRootedTree.qza` |
| `--notree` | `orientedTable.qza`, `orientedRepSeqs.qza` |
| All workflows | `gg2Taxonomy.qza` or `silvaTaxonomy.qza`, `taxa_read_counts.csv`, `taxa_read_losses.csv`, configuration and input records |

## Test datasets

Test files contain real accessions and require an internet connection to download complete Runs. Runtime and disk usage depend on the data size. Use a separate working directory; tests also use the fixed `results/pip`, `results/taxa`, and `results/db` paths.

```bash
# Built-in test: 6 projects and 6 Runs across multiple platforms
Meta2Data AmpliconPIP --test --vsearch -t 8

# Classify the test results without building a tree
Meta2Data AmpliconTAXA --vsearch --notree -t 8
```

The extended list, [newpipe.csv](test/newpipe.csv), contains 12 projects and 24 Runs. Save it in your working directory and run:

```bash
Meta2Data AmpliconPIP --test --vsearch \
  --public-m newpipe.csv \
  --public-bioproject-colNAME Bioproject --public-sra-colNAME Run -t 8
```

`--test` is for public data only and cannot be combined with `--local-m`. Without `--public-m`, it uses the built-in `test/ampliconpiptest.csv`; with a CSV, it selects the first 2 rows per project. When combined with `--public-bioprojectIDs`, projects are selected first, then the first 2 rows are taken; matching primers are still required. A processing method must be selected. Skipping unsupported platforms is expected when using dada2.

MetaDL has no `--test` option. To test it, prepare a small number of IDs as described in ID-file mode and pass their directory with `-i`. Test TAXA using results already generated by PIP.

## Advanced parameters and diagnostics

For 454, after primer processing and tail trimming, the pipeline calculates the median read length `L` from all reads in each Run/local sample. It removes reads shorter than `ceil(0.5 × L)`, then removes reads containing more than 1 `N`/`n` in total. Retained reads are not truncated to a common length. Set the fraction with JSON parameter `ls454.length_fraction` (default `0.5`) and the maximum N count with `ls454.max_n` (default `1`). This screens for length and ambiguous bases without relying on quality scores; it does not guarantee error-free sequences.

Filtering is followed by exact dereplication, 99% preclustering, vsearch denoising, chimera removal, 97% clustering, and mapping reads back to the representatives. Sequences with abundance 1 are retained before 99% preclustering; denoising still uses `vsearch.minsize`. Mapping uses filtered reads, excluding members assigned to identified chimeras. There is no additional removal of final features with count 1. Per-sample thresholds and losses are recorded in `ls454_quality-vsearch.json`; member exclusion statistics are in `ls454_members-vsearch.json`.

For Ion Torrent, both vsearch and dada2 automatically set `EE = L × 10^(-19/10)` by default, where `L` is the dataset's median read length after primer removal. The actual threshold is recorded in `ion_quality-<method>.json` in the dataset directory.

For Oxford Nanopore, reads and representative sequences are ordered deterministically before clustering and consensus polishing. Sequence identifiers are stable, and repeated reads retain their full contribution to abundance.

PIP and TAXA both support `--parameter settings.json`. Include only the fields to override; explicit CLI options take precedence. See the [JSON parameter guide](docs/parameters.md) for definitions, units, ranges, applicable platforms, and effects of parameter changes. The [default configuration](docs/parameters.default.json) is available as a template.

For example, save the following as `settings.json`:

```json
{
  "primer": {"skip_unknown": true},
  "vsearch": {"maxee": 1.0, "cluster_identity": 0.97},
  "taxa": {"confidence": 0.8}
}
```

For local input, the dataset name, path, and platform column parameters remain required when using a configuration file:

```bash
Meta2Data AmpliconPIP --local-m local_metadata.csv \
  --local-datasets-colNAME datasets \
  --local-path-colNAME path --local-platform-colNAME platform \
  --vsearch --parameter settings.json -t 8
Meta2Data AmpliconTAXA --vsearch --parameter settings.json --notree -t 8
```

AmpliconPIP adapter protection automatically uses the GG2 sequence reference. The BLAST index and query cache are stored in `results/db/adapter_guard/` under the launch directory.

If a dataset fails, check `datasets.log` and its dataset log first, resolve the issue, and rerun the same command. In `summary.csv`, `0` means that no reads were retained at that stage; `NA` means the stage was not run or was not applicable. Initial counts for paired-end data are reported as read pairs.

Further documentation: [Processing workflow](docs/amplicon_pipeline.md) · [Adapter protection](docs/adapter_guard.md) · [Read counts](docs/read_counts.md) · [Resource statistics](docs/resource_profile.md).
