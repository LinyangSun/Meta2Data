# AmpliconFunction.sh - Unified Amplicon Processing Functions

## Overview

`AmpliconFunction.sh` consolidates all amplicon sequencing processing functions from individual scripts into a single unified library.

## File Statistics

- **Total Lines**: 646
- **Total Functions**: 13 implemented + 3 ONT placeholders
- **Platforms Supported**: 454, Illumina, PacBio, ONT (partial)

## Function List

### Common Functions (3)
Functions shared across all sequencing platforms:

1. `Amplicon_Common_MakeManifestFileForQiime2()` - Generate QIIME2 manifest files
2. `Amplicon_Common_ImportFastqToQiime2()` - Import FASTQ to QIIME2 format
3. `Amplicon_Common_FinalFilesCleaning()` - Clean up and organize final outputs

### 454 Platform Functions (5)
Roche 454 pyrosequencing pipeline:

4. `Amplicon_454_ImportToQiime2()` - Import 454 data to QIIME2
5. `Amplicon_454_QualityControl()` - Quality filtering for 454 reads
6. `Amplicon_454_Deduplication()` - Remove duplicate sequences
7. `Amplicon_454_ChimerasRemoval()` - Detect and remove chimeric sequences
8. `Amplicon_454_ClusterDenovo()` - De novo OTU clustering

### Illumina Platform Functions (3)
Illumina short-read sequencing pipeline:

9. `Amplicon_Illumina_PrimerDetectionAndQualityControl()` - Primer trimming and QC
10. `Amplicon_Illumina_QualityControlForQZA()` - Quality control for imported data
11. `Amplicon_Illumina_DenosingDeblur()` - Denoising with Deblur algorithm

### PacBio Platform Functions (3)
PacBio HiFi/CCS long-read sequencing pipeline:

12. `Amplicon_Pacbio_PrimerDetectionAndQualityControl()` - Primer detection for long reads
13. `Amplicon_Pacbio_QualityControlForQZA()` - Quality control for PacBio data
14. `Amplicon_Pacbio_DenosingDada2()` - Denoising with DADA2

### ONT Platform Functions (3) - TODO
Oxford Nanopore long-read sequencing pipeline (placeholders):

- `Amplicon_ONT_PrimerDetectionAndQualityControl()` - TODO: Not implemented
- `Amplicon_ONT_QualityControlForQZA()` - TODO: Not implemented  
- `Amplicon_ONT_DenosingMethod()` - TODO: Not implemented

## Usage

### In Pipeline Scripts

Source this file at the beginning of your pipeline script:

```bash
#!/bin/bash
source scripts/AmpliconFunction.sh

# Now you can call any function
Amplicon_Common_MakeManifestFileForQiime2
Amplicon_Illumina_PrimerDetectionAndQualityControl
# ... etc
```

### Required Environment Variables

Most functions expect these variables to be set:

- `dataset_path` - Path to the dataset directory
- `dataset_name` - Name of the dataset
- `sequence_type` - Either "single" or "paired"
- `F_primer` - Forward primer sequence (platform-specific)
- `R_primer` - Reverse primer sequence (platform-specific)
- `cpu` - Number of CPU threads to use

### Example Usage

```bash
#!/bin/bash
set -e

# Source the unified functions
source scripts/AmpliconFunction.sh

# Set required variables
export dataset_path="/path/to/dataset/"
export dataset_name="PRJNA12345"
export sequence_type="paired"
export cpu=8

# Run Illumina pipeline
Amplicon_Illumina_PrimerDetectionAndQualityControl
Amplicon_Common_MakeManifestFileForQiime2
Amplicon_Common_ImportFastqToQiime2
Amplicon_Illumina_QualityControlForQZA
Amplicon_Illumina_DenosingDeblur
Amplicon_Common_FinalFilesCleaning
```

## Original Scripts Consolidated

This file consolidates the following original scripts:

### Common
- `Amplicon_Common_MakeManifestFileForQiime2.sh`
- `Amplicon_Common_ImportFastqToQiime2.sh`
- `Amplicon_Common_FinalFilesCleaning.sh`

### 454
- `Amplicon_454_ImportToQiime2.sh`
- `Amplicon_454_QualityControl.sh`
- `Amplicon_454_Deduplication.sh`
- `Amplicon_454_ChimerasRemoval.sh`
- `Amplicon_454_ClusterDenovo.sh`

### Illumina
- `Amplicon_Illumina_PrimerDetectionAndQualityControl.sh`
- `Amplicon_Illumina_QualityControlForQZA.sh`
- `Amplicon_Illumina_DenosingDeblur.sh`

### PacBio
- `Amplicon_Pacbio_PrimerDetectionAndQualityControl.sh`
- `Amplicon_Pacbio_QualityControlForQZA.sh`
- `Amplicon_Pacbio_DenosingDada2.sh`

## Maintenance

To update functions:
1. Edit `AmpliconFunction.sh` directly, or
2. Edit individual scripts and regenerate using the merge script

## TODO

- [ ] Implement ONT platform functions
- [ ] Add multi-platform merging functions (vsearch-based)
- [ ] Add multi-software integration functions
- [ ] Add comprehensive error handling
- [ ] Add logging functionality

## Notes

- All functions maintain their original behavior and dependencies
- Functions use `set -e` for error handling
- External dependencies: QIIME2, vsearch, fastp, py_16s.py

