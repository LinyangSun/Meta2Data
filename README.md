<p align="center">
  <img src="docs/llloooooog.png" alt="Meta2Data logo" width="100%">
</p>

# Meta2Data

> [!IMPORTANT]
> **Help expand Meta2Data | 欢迎共同拓展 Meta2Data**
>
> Meta2Data's analysis workflow currently supports only 16S rRNA amplicon data. We aim to extend it to other marker genes and short- and long-read metagenomic data. **Contributions and collaborations are welcome!**
>
> Meta2Data 的分析流程目前仅支持 16S rRNA 扩增子数据。未来希望扩展至其他标记基因，以及短读长和长读长宏基因组数据。**欢迎参与开发与合作！**

Meta2Data is a Linux toolkit for amplicon analysis, combining public metadata retrieval, sequencing data processing, taxonomic classification, and phylogenetic tree construction. It also supports local FASTQ datasets.

## Citation

If you use Meta2Data in your research, please cite:

**Li J, Duplouy A, Sun L. (2026).** Meta2Data: Automated data processing for large scale microbiome data re-analysis. *Bioinformatics*, btag828.  
[https://doi.org/10.1093/bioinformatics/btag828](https://doi.org/10.1093/bioinformatics/btag828)

**Please also cite the software and reference databases used in your analysis.**

| Workflow | Core software to cite |
|---|---|
| **VSEARCH workflow** (`--vsearch`) | [fastp](https://doi.org/10.1002/imt2.70078), [VSEARCH](https://doi.org/10.7717/peerj.2584), and [QIIME 2](https://doi.org/10.1038/s41587-019-0209-9) |
| **DADA2 workflow** (`--dada2`) | [fastp](https://doi.org/10.1002/imt2.70078), [DADA2](https://doi.org/10.1038/nmeth.3869), and [QIIME 2](https://doi.org/10.1038/s41587-019-0209-9) |

Additional citations apply to AmpliconTAXA, platform-specific processing, and the optional steps used in your workflow.

<details>
<summary><strong>Additional citations by analysis step</strong></summary>

Use the references for the steps and database versions used in your run. The core software references are also listed in full below.

| Step or resource used | Additional references |
|---|---|
| Oxford Nanopore processing with `--vsearch` | [Chopper / NanoPack2](https://doi.org/10.1093/bioinformatics/btad311), [minimap2](https://doi.org/10.1093/bioinformatics/bty191), and [Racon](https://doi.org/10.1101/gr.214270.116) |
| Primer matching or trimming with Cutadapt | [Cutadapt](https://doi.org/10.14806/ej.17.1.200) |
| Adapter protection using BLAST and the GG2 reference | [BLAST+](https://doi.org/10.1186/1471-2105-10-421) and [Greengenes2](https://doi.org/10.1038/s41587-023-01845-1) |
| AmpliconTAXA sequence orientation and taxonomic classification, with either denoising method | [RESCRIPt](https://doi.org/10.1371/journal.pcbi.1009581), [VSEARCH](https://doi.org/10.7717/peerj.2584), [q2-feature-classifier](https://doi.org/10.1186/s40168-018-0470-z), and [Greengenes2](https://doi.org/10.1038/s41587-023-01845-1) |
| SILVA classifier (`--classifier silva-138.1` or `--classifier silva-144`) | [SILVA (Chuvochina et al., 2026)](https://doi.org/10.1093/nar/gkaf1247) |
| Default SEPP tree construction | [SEPP](https://doi.org/10.1142/9789814366496_0024), [q2-fragment-insertion](https://doi.org/10.1128/mSystems.00021-18), and the Greengenes 13_8 reference: [DeSantis et al. (2006)](https://doi.org/10.1128/AEM.03006-05) and [McDonald et al. (2012)](https://doi.org/10.1038/ismej.2011.139) |
| De novo tree construction with `--singleV` | [MAFFT](https://doi.org/10.1093/molbev/mst010) and [FastTree](https://doi.org/10.1371/journal.pone.0009490) |
| Full-length PacBio CCS processing with DADA2 `denoise-ccs` | [Callahan et al. (2019)](https://doi.org/10.1093/nar/gkz569), in addition to the DADA2 reference above |
| VSEARCH denoising (`--cluster_unoise`) | [UNOISE (Edgar, 2016)](https://doi.org/10.1101/081257), in addition to VSEARCH |
| VSEARCH chimera removal (`--uchime3_denovo`) | [UCHIME2 (Edgar, 2016)](https://doi.org/10.1101/074252), in addition to VSEARCH |

AmpliconTAXA uses RESCRIPt and VSEARCH to orient sequences against Greengenes2 even when the SILVA classifier is selected. The Greengenes 13_8 reference used by SEPP is a separate resource from Greengenes2. With `--notree`, no tree-construction references are needed.

The SILVA reference follows the database's [current citation guidance](https://www.arb-silva.de/contact/citation). Report the database release selected in your run. The legacy `silva` option refers to SILVA 138.1: the original official file was named `silva-138-99-nb-classifier.qza`, but its embedded training provenance identifies release 138.1.

### Core software references

- **fastp:** Chen S. (2025). fastp 1.0: An ultra-fast all-round tool for FASTQ data quality control and preprocessing. *iMeta 4(5):e70078*. [doi:10.1002/imt2.70078](https://doi.org/10.1002/imt2.70078)
- **VSEARCH:** Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016). VSEARCH: a versatile open source tool for metagenomics. *PeerJ 4:e2584*. [doi:10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)
- **DADA2:** Callahan BJ et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods 13:581–583*. [doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
- **QIIME 2:** Bolyen E et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology 37:852–857*. [doi:10.1038/s41587-019-0209-9](https://doi.org/10.1038/s41587-019-0209-9)

### References for additional steps

- **Chopper / NanoPack2:** De Coster W, Rademakers R. (2023). NanoPack2: population-scale evaluation of long-read sequencing data. *Bioinformatics 39(5):btad311*. [doi:10.1093/bioinformatics/btad311](https://doi.org/10.1093/bioinformatics/btad311)
- **minimap2:** Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics 34(18):3094–3100*. [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)
- **Racon:** Vaser R, Sović I, Nagarajan N, Šikić M. (2017). Fast and accurate de novo genome assembly from long uncorrected reads. *Genome Research 27(5):737–746*. [doi:10.1101/gr.214270.116](https://doi.org/10.1101/gr.214270.116)
- **Cutadapt:** Martin M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal 17(1):10–12*. [doi:10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)
- **BLAST+:** Camacho C et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics 10:421*. [doi:10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)
- **RESCRIPt:** Robeson MS II et al. (2021). RESCRIPt: Reproducible sequence taxonomy reference database management. *PLoS Computational Biology 17(11):e1009581*. [doi:10.1371/journal.pcbi.1009581](https://doi.org/10.1371/journal.pcbi.1009581)
- **q2-feature-classifier:** Bokulich NA et al. (2018). Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2’s q2-feature-classifier plugin. *Microbiome 6:90*. [doi:10.1186/s40168-018-0470-z](https://doi.org/10.1186/s40168-018-0470-z)
- **Greengenes2:** McDonald D et al. (2024). Greengenes2 unifies microbial data in a single reference tree. *Nature Biotechnology 42:715–718*. [doi:10.1038/s41587-023-01845-1](https://doi.org/10.1038/s41587-023-01845-1)
- **SILVA:** Chuvochina M et al. (2026). SILVA in 2026: a global core biodata resource for rRNA within the DSMZ digital diversity. *Nucleic Acids Research 54(D1):D334–D341*. [doi:10.1093/nar/gkaf1247](https://doi.org/10.1093/nar/gkaf1247)
- **SEPP:** Mirarab S, Nguyen N, Warnow T. (2012). SEPP: SATé-enabled phylogenetic placement. *Pacific Symposium on Biocomputing:247–258*. [doi:10.1142/9789814366496_0024](https://doi.org/10.1142/9789814366496_0024)
- **q2-fragment-insertion:** Janssen S et al. (2018). Phylogenetic Placement of Exact Amplicon Sequences Improves Associations with Clinical Information. *mSystems 3(3):e00021-18*. [doi:10.1128/mSystems.00021-18](https://doi.org/10.1128/mSystems.00021-18)
- **Greengenes database:** DeSantis TZ et al. (2006). Greengenes, a Chimera-Checked 16S rRNA Gene Database and Workbench Compatible with ARB. *Applied and Environmental Microbiology 72(7):5069–5072*. [doi:10.1128/AEM.03006-05](https://doi.org/10.1128/AEM.03006-05)
- **Greengenes taxonomy:** McDonald D et al. (2012). An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea. *The ISME Journal 6:610–618*. [doi:10.1038/ismej.2011.139](https://doi.org/10.1038/ismej.2011.139)
- **MAFFT:** Katoh K, Standley DM. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution 30(4):772–780*. [doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)
- **FastTree:** Price MN, Dehal PS, Arkin AP. (2010). FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments. *PLoS ONE 5(3):e9490*. [doi:10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490)
- **DADA2 for PacBio CCS:** Callahan BJ et al. (2019). High-throughput amplicon sequencing of the full-length 16S rRNA gene with single-nucleotide resolution. *Nucleic Acids Research 47(18):e103*. [doi:10.1093/nar/gkz569](https://doi.org/10.1093/nar/gkz569)
- **UNOISE:** Edgar RC. (2016). UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing. *bioRxiv 081257 (preprint)*. [doi:10.1101/081257](https://doi.org/10.1101/081257)
- **UCHIME2:** Edgar RC. (2016). UCHIME2: improved chimera prediction for amplicon sequencing. *bioRxiv 074252 (preprint)*. [doi:10.1101/074252](https://doi.org/10.1101/074252)

### Earlier fastp versions

The main fastp link follows the current citation for fastp 1.x. For an older environment, choose the paper appropriate to the version and features used:

- **fastp, original paper:** Chen S, Zhou Y, Chen Y, Gu J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics 34(17):i884–i890*. [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
- **fastp, subsequent improvements:** Chen S. (2023). Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. *iMeta 2(2):e107*. [doi:10.1002/imt2.107](https://doi.org/10.1002/imt2.107)

</details>

## New features

### 20 November 2026: SILVA 144 support

- **SILVA 144:** AmpliconTAXA now supports the official full-length SILVA 144 classifier. Select it with `--classifier silva-144`.

- **Updated environment and SIF v1.1:** The runtime has been upgraded to **QIIME 2 2026.7.0** and **scikit-learn 1.7.1**, available through both Conda installation and the new **SIF v1.1** image.

- **Compatible GG2 and SILVA 138.1 classifiers:** Full-length Greengenes2 2024.09 and SILVA 138.1 classifiers have been retrained for the updated runtime and are downloaded automatically when needed. The existing `greengenes` and `silva` aliases continue to select GG2 2024.09 and SILVA 138.1, respectively.

See [Installation](#installation) for setup instructions and [Built-in classifiers](#built-in-classifiers) for database options.

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

The compatible Greengenes2 2024.09 and SILVA 138.1 replacements have passed
their model audits and are published. Fresh native validation passed 317 unit
tests with no failures, errors or skips, 14 PIP cases and five TAXA cases.
The final SIF passed two Illumina PE PIP cases (DADA2 and VSEARCH) and five
TAXA cases. See the [final validation record](docs/validation/2026-11-19/release-gates.md)
for the tested source, exact scope and evidence.

### Conda

```bash
export M2D_ROOT=/absolute/path/to/Meta2Data  # Source directory
export M2D_ENV=/absolute/path/to/meta2data-2026.7  # New environment directory
python3 "$M2D_ROOT/scripts/create_conda_environment.py" --prefix "$M2D_ENV"
conda activate "$M2D_ENV"
export PATH="$M2D_ROOT/bin:$PATH"
Meta2Data --help
```

`env.yml` lists the direct requirements for Linux x86_64 with glibc 2.28 or
newer: QIIME 2 2026.7.0, the eight plugins Meta2Data uses, scikit-learn 1.7.1,
and the pipeline's Python modules and command-line tools.
`environment-constraints-linux-64.txt` retains the official release dependency
versions and their resolved builds; `environment-linux-64.lock`
records their archives and checksums. Unused plugins such as `q2-fondue` are
excluded. The same files are used for Conda and the fresh SIF.
Create a separate environment when upgrading. Classifier models built with
scikit-learn 1.4.2 cannot be reused in this environment.

The installer requires Conda with the libmamba solver. It checks that the lock
matches `env.yml` and every release constraint, installs the exact package archives with SHA256 checksums,
and excludes configured extra channels and default packages. This also avoids the
[Conda 26.7–26.9 `nodefaults` regression](https://github.com/conda/conda/issues/16816).
It does not change your Conda configuration or replace an existing environment.
It then installs the checksum-pinned official NCBI SRA Toolkit 3.4.1 under the
new environment and exposes `prefetch` and `fasterq-dump` there. A provenance
record is retained with that component. Use `--conda /path/to/conda` if Conda is
not on PATH; `--dry-run` checks the lock and Conda command without downloading
or installing either component. A dry run does not verify file ownership or
runtime compatibility.

### SIF / Apptainer

The validated `v1.1` image is published in
[GitHub Packages](https://github.com/LinyangSun/Meta2Data/pkgs/container/meta2data).
Install Apptainer, download the image and define the `m2d` shortcut below.
Meta2Data and its dependencies run directly from the image.

```bash
apptainer pull Meta2Data-v1.1.sif oras://ghcr.io/linyangsun/meta2data:v1.1

export M2D_SIF="$(pwd)/Meta2Data-v1.1.sif"

m2d() {
    apptainer exec --cleanenv "$M2D_SIF" Meta2Data "$@"
}

m2d --help
```

```bash
m2d AmpliconPIP \
    --public-m metadata.csv \
    --public-bioproject-colNAME Bioproject \
    --public-sra-colNAME Run \
    --vsearch \
    --threads 8 \
    --max-parallel 2
```

Run commands from your working directory; PIP and TAXA write to its `results/` directory. Set the image path and define the `m2d` function again in each new terminal.

For SIF, replace the leading `Meta2Data` with `m2d` in the examples below.


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
| `--classifier NAME` | Built-in classifier: `gg2-2024.09`, `silva-138.1`, or `silva-144`; legacy aliases `greengenes` and `silva` are retained. Default: GG2 2024.09 |
| `--singleV` | For data from the same amplified region: align sequences and build a tree de novo |
| `--notree` / `--no-tree` | Merge, orient, and classify without building a tree |
| `--confidence FLOAT` | Classification confidence; default `0.7`, range `0–1` |
| `-t, --threads N` | Threads; default `4` |
| `--parameter FILE` | Override parameters using JSON |
| `-h, --help` | Show help |

All modes orient sequences against the GG2 sequence reference and remove features that cannot be oriented from both the representative sequences and feature table. The default SEPP workflow also removes features that cannot be inserted into the tree. `--singleV` and `--notree` are mutually exclusive, and neither requires the SEPP reference. `--notree` skips tree construction and tree-insertion filtering only. Classifier selection is independent of the tree method; resources are prepared automatically in `results/db/`.

### Built-in classifiers

| Database release | Explicit selection | Legacy alias | Taxonomy output |
|---|---|---|---|
| Greengenes2 2024.09 | `--classifier gg2-2024.09` (default) | `--classifier greengenes` | `gg2Taxonomy.qza` |
| SILVA 138.1 | `--classifier silva-138.1` | `--classifier silva` | `silvaTaxonomy.qza` |
| SILVA 144 | `--classifier silva-144` | None | `silva144Taxonomy.qza` |

All three choices use full-length reference classifiers built with scikit-learn 1.7.1.
The selected model is prepared automatically in `results/db/`, reusing a valid
cache or bundle and downloading from a fixed release URL when needed.
Conda and SIF builds use the same resource definitions.

Input classifier files are cached in `results/db/`. These model artifacts are
separate from the taxonomy output artifacts listed above:

| Canonical classifier | Input model filename in `results/db/` |
|---|---|
| `gg2-2024.09` | `gg2-2024.09-full-length-sklearn-1.7.1.qza` |
| `silva-138.1` | `silva-138.1-full-length-sklearn-1.7.1.qza` |
| `silva-144` | `silva-144-full-length-sklearn-1.7.1.qza` |

The `silva` alias remains tied to 138.1; it does not select the latest SILVA
release. Classifier choice is a command-line option, not a JSON parameter or a
path to a user-supplied `.qza` file. The selected database identity is recorded in
the TAXA state. Classifier builds have separate cache filenames, so older
scikit-learn 1.4.2 files cannot be mistaken for the replacement models.
Before reusing a taxonomy result, Meta2Data checks its recorded classifier UUID,
input-sequence UUID, and confidence. A valid archive with mismatched provenance
is recomputed even when its filename and saved run state appear correct.

Taxonomy strings are preserved as produced by the selected classifier. Meta2Data
does not infer ranks from fixed field positions. When using downstream tools,
check the selected classifier's rank scheme before choosing a numeric level.

If the recorded state identifies results from an input source that has since been replaced, TAXA reports an error and preserves the files. Rerun PIP with that method first. External results without a state record can still be imported.

To classify with SILVA 144 without building a tree:

```bash
Meta2Data AmpliconTAXA --dada2 --classifier silva-144 --notree -t 8
```

Results are saved in `results/taxa/final-<method>-<multiV|singleV|notree>/`:

| Workflow | Main output |
|---|---|
| Default, multiple regions | `treeFilteredTable.qza`, `treeFilteredRepSeqs.qza`, `seppTree.qza` |
| `--singleV` | `orientedTable.qza`, `orientedRepSeqs.qza`, `denovoRootedTree.qza` |
| `--notree` | `orientedTable.qza`, `orientedRepSeqs.qza` |
| All workflows | Selected `gg2Taxonomy.qza`, `silvaTaxonomy.qza`, or `silva144Taxonomy.qza`; `taxa_read_counts.csv`, `taxa_read_losses.csv`, configuration and input records |

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

For normal-quality Illumina paired-end data in `--vsearch`, an authorized R1
5-prime primer cut is also applied once to R2. Unknown-primer skip rules take
precedence. Each mate then loses only its contiguous Q<20 tail. All samples
finish merging before one project-wide route is chosen: >=50% pooled merge
success uses merged reads; otherwise all usable R1 reads are retained for
filtering, including successfully merged pairs and R1 with unusable R2. The
selected input is filtered at EE<=1.5 (`vsearch.illumina_pe_maxee`) without EE
truncation. Other platforms, DADA2, original single-end input and degraded-quality
branches are unchanged. Detailed counts and logs survive temporary-file cleanup;
see [the pipeline description](docs/amplicon_pipeline.md#normal-quality-illumina-paired-end-vsearch).

For Ion Torrent, both vsearch and dada2 automatically set `EE = L × 10^(-19/10)` by default, where `L` is the dataset's median read length after primer removal. The actual threshold is recorded in `ion_quality-<method>.json` in the dataset directory.

For Oxford Nanopore, reads and representative sequences are ordered deterministically before clustering and consensus polishing. Sequence identifiers are stable, and repeated reads retain their full contribution to abundance.

PIP and TAXA both support `--parameter settings.json`. Include only the fields to override; explicit CLI options take precedence. See the [JSON parameter guide](docs/parameters.md) for definitions, units, ranges, applicable platforms, and effects of parameter changes. The [default configuration](docs/parameters.default.json) is available as a template.

For example, save the following as `settings.json`:

```json
{
  "primer": {"skip_unknown": true},
  "vsearch": {"maxee": 1.0, "illumina_pe_maxee": 1.5, "cluster_identity": 0.97},
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
