# MetaDL - Multi-Database Metadata Download Tool

Automated tool for downloading and merging BioProject metadata from NCBI and CNCB/GSA databases.

**🆕 NEW: Unified pipeline now supports both NCBI (PRJNA/PRJEB/PRJDB) and CNCB/GSA (PRJCA) databases!**

For comprehensive documentation on the unified pipeline, see [UNIFIED-METADL-README.md](UNIFIED-METADL-README.md).

## Quick Start

### Installation

```bash
cd /Users/a1-6/Project-helsinki/Meta2Data

# Install Python dependencies
pip install -r requirements_metadl.txt
```

### Basic Usage

```bash
# Download metadata for BioProjects
./bin/Meta2Data MetaDL \
    -i /path/to/bioproject_ids/ \
    -o /path/to/output/ \
    -e your@email.com
```

## Prerequisites

- Python 3.7+
- biopython >= 1.80 (NCBI API access)
- pandas >= 2.0.0 (data processing)
- selenium >= 4.17.0 (CNCB/GSA downloads)
- openpyxl >= 3.0.0 (Excel parsing)
- ChromeDriver (for CNCB downloads: `brew install --cask chromedriver`)

## Input Format

### BioProject ID Files

Create `.txt` files containing BioProject IDs (one per line):

**Example: `my_projects.txt`**
```
PRJNA123456
PRJEB234567
PRJCA345678
```

**Supported formats:**
- `PRJNA*` - NCBI (US)
- `PRJEB*` - ENA (Europe)
- `PRJDB*` - DDBJ (Japan)
- `PRJCA*` - CNCB/GSA (China) ← NEW!

Multiple `.txt` files can be placed in the input directory - all will be processed.

## Workflow

The unified tool automatically performs these steps:

1. **Read and classify BioProject IDs** (NCBI vs CNCB/GSA)
2. **Download NCBI data** (BioSample + SRA) for PRJNA/PRJEB/PRJDB
3. **Download CNCB data** (Excel metadata) for PRJCA
4. **Parse and merge NCBI data** (SRA + BioSample)
5. **Parse and convert CNCB data** (Excel → CSV)
6. **Final merge** (NCBI + CNCB with normalized columns)

## Output Files

### Individual Project Files

- `PRJNA*_biosample.txt` - Raw BioSample data for each project
- `PRJNA*_sra_runinfo.csv` - Raw SRA RunInfo for each project

### Combined Files

- `all_biosamples_combined.csv` - All BioSample metadata combined
- `all_sra_runs_combined.csv` - All SRA metadata combined

### Merged Files (Primary Output)

**NCBI-only files:**
- `ncbi_merged_*.csv` - Merged NCBI data (SRA + BioSample)

**CNCB-only files:**
- `cncb_combined.csv` - Combined CNCB data (converted from Excel)

**Final output:**
- **`all_metadata_merged.csv`** - **Main output**: Combined NCBI + CNCB data with normalized columns

## Usage Examples

### Example 1: Single Keyword Search

```bash
# 1. Create input directory
mkdir -p bioproject_ids

# 2. Create BioProject ID file
cat > bioproject_ids/microbiome_projects.txt << EOF
PRJNA123456
PRJNA234567
EOF

# 3. Run MetaDL
./bin/Meta2Data MetaDL \
    -i bioproject_ids/ \
    -o metadata_output/ \
    -e researcher@university.edu
```

### Example 2: Multiple Keyword Searches

```bash
# Organize by keywords
mkdir -p bioproject_ids

# Gut microbiome projects
cat > bioproject_ids/gut_microbiome.txt << EOF
PRJNA123456
PRJNA234567
EOF

# Soil microbiome projects
cat > bioproject_ids/soil_microbiome.txt << EOF
PRJNA345678
PRJNA456789
EOF

# Run MetaDL (will process all txt files)
./bin/Meta2Data MetaDL \
    -i bioproject_ids/ \
    -o metadata_output/ \
    -e researcher@university.edu
```

### Example 3: Large Dataset

```bash
# For many projects, process in batches
./bin/Meta2Data MetaDL \
    -i large_dataset/batch1/ \
    -o output/batch1/ \
    -e researcher@university.edu

./bin/Meta2Data MetaDL \
    -i large_dataset/batch2/ \
    -o output/batch2/ \
    -e researcher@university.edu
```

## Getting BioProject IDs from NCBI

### Method 1: NCBI Website (Recommended)

1. Go to [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/)
2. Search with your keywords (e.g., "human gut microbiome")
3. Click "Send to" → "File" → Format: "Accession List"
4. Download the file
5. Place in your input directory

### Method 2: Manual Collection

Create a text file with one BioProject ID per line:

```bash
cat > my_projects.txt << EOF
PRJNA123456
PRJNA234567
PRJNA345678
EOF
```

## Output Data Structure

### BioSample Fields (example)

- `BioSample` - BioSample accession (e.g., SAMN12345678)
- `BioProject_ID` - Associated BioProject
- `Organism` - Source organism
- `Sample_Name` - Sample identifier
- `host` - Host organism (for host-associated samples)
- `isolation_source` - Sample origin
- `collection_date` - When sample was collected
- `geo_loc_name` - Geographic location
- Custom attributes specific to the study

### SRA Fields (example)

- `Run` - SRA run accession (e.g., SRR12345678)
- `BioSample` - Associated BioSample
- `BioProject` - Associated BioProject
- `LibraryStrategy` - Sequencing strategy (AMPLICON, WGS, etc.)
- `LibraryLayout` - SINGLE or PAIRED
- `Platform` - Sequencing platform (ILLUMINA, PACBIO, etc.)
- `Model` - Sequencer model
- `spots` - Number of spots
- `bases` - Total bases
- `avgLength` - Average read length
- `download_path` - FTP path to FASTQ files

## Performance

- **Small dataset** (1-10 BioProjects): 1-5 minutes
- **Medium dataset** (10-100 BioProjects): 10-60 minutes
- **Large dataset** (100+ BioProjects): 1+ hours

Download time depends on:
- Number of BioProjects
- Number of samples per BioProject
- NCBI server response time
- Network speed

## Troubleshooting

### Error: "No .txt files found"

**Solution**: Ensure your input directory contains `.txt` files with BioProject IDs

```bash
ls -l bioproject_ids/
# Should show .txt files
```

### Error: "No BioProject IDs found"

**Solution**: Check that your txt files contain valid BioProject IDs

```bash
cat bioproject_ids/my_projects.txt
# Should show PRJNA*, PRJEB*, PRJDB*, or PRJCA* IDs
```

### Slow Downloads

**Solution**:
- Process in smaller batches
- Check network connection
- Try during off-peak hours

### NCBI API Errors

**Solution**:
- Verify email address is correct
- Check NCBI website is accessible
- Wait and retry (NCBI may be experiencing issues)

## Integration with AmpliconPIP

The merged metadata can be converted to AmpliconPIP format:

```bash
# 1. Download metadata with MetaDL
./bin/Meta2Data MetaDL -i bioproject_ids/ -o metadata/ -e your@email.com

# 2. Convert matched_sra_biosample_*.csv to AmpliconPIP format
# (Manual step: select relevant columns and rename to match AmpliconPIP requirements)

# 3. Run AmpliconPIP
./bin/Meta2Data AmpliconPIP -m converted_metadata.csv -o results/ -t 8
```

## Advanced Usage

### Python API

The core downloader can be used as a Python module:

```python
from scripts.MetaDL.ncbi_metadata_downloader import run_complete_pipeline

results = run_complete_pipeline(
    input_folder='bioproject_ids/',
    output_folder='output/',
    email='your@email.com'
)

# Access results
bioprojects = results['bioproject_ids']
biosamples_df = results['biosamples']
sra_df = results['sra_runs']
merged_df = results['merge_results']['matched']
```

### Custom Processing

Modify `scripts/MetaDL/ncbi_metadata_downloader.py` to:
- Add custom parsers for specific attributes
- Implement additional quality filters
- Export to different formats

## Files Structure

```
Meta2Data/
├── bin/
│   └── MetaDL                          # Main command wrapper
├── scripts/
│   └── MetaDL/
│       ├── ncbi_metadata_downloader.py # Core Python script
│       ├── getmeta.ipynb               # Original Jupyter notebook
│       └── TestData/                   # Test files
└── requirements_metadl.txt             # Python dependencies
```

## Citation

If you use MetaDL in your research, please cite:
- Meta2Data pipeline
- NCBI databases (BioProject, BioSample, SRA)

## Support

For issues or questions:
1. Check this README
2. Review error messages
3. Verify input format
4. Check NCBI website accessibility

---

**Last Updated**: 2025-12-09
