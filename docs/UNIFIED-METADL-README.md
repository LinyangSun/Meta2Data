# Unified MetaDL - Multi-Database Metadata Download Tool

Automated tool for downloading and merging BioProject metadata from **NCBI** (PRJNA/PRJEB/PRJDB) and **CNCB/GSA** (PRJCA) databases.

## Features

- **Multi-database support**: NCBI, ENA, DDBJ, and CNCB/GSA
- **Automatic routing**: Classifies BioProjects by prefix and routes to appropriate downloader
- **Intelligent merging**: Three-stage merge strategy for optimal data integration
- **Column normalization**: Standardizes field names across databases
- **Comprehensive output**: Individual project files + combined datasets

## Quick Start

### Installation

```bash
cd /Users/a1-6/Project-helsinki/Meta2Data

# Install dependencies (includes NCBI and CNCB tools)
pip install -r requirements_metadl.txt

# Install ChromeDriver for CNCB downloads (macOS)
brew install --cask chromedriver
```

### Basic Usage

```bash
# Download metadata for mixed NCBI and CNCB BioProjects
./bin/Meta2Data MetaDL \
    -i /path/to/bioproject_ids/ \
    -o /path/to/output/ \
    -e your@email.com
```

## Input Format

### BioProject ID Files

Create `.txt` files containing BioProject IDs (one per line):

**Example: `mixed_projects.txt`**
```
PRJNA123456
PRJEB234567
PRJCA345678
PRJDB456789
```

**Supported formats:**
- `PRJNA*` - NCBI (US)
- `PRJEB*` - ENA (Europe)
- `PRJDB*` - DDBJ (Japan)
- `PRJCA*` - CNCB/GSA (China)

Multiple `.txt` files can be placed in the input directory.

## Three-Stage Merge Workflow

The pipeline performs intelligent data integration in three stages:

### Stage 1: NCBI Data Merge (SRA + BioSample)

For NCBI BioProjects (PRJNA/PRJEB/PRJDB):

1. **Download** BioSample metadata via Entrez API
2. **Download** SRA RunInfo via Entrez API
3. **Merge** using two strategies:
   - Strategy A: Match by BioSample ID only
   - Strategy B: Match by BioProject + BioSample composite key
   - Auto-selects best strategy based on match rate

**Output:** `ncbi_merged_*.csv`

### Stage 2: CNCB Data Conversion (Excel → CSV)

For CNCB/GSA BioProjects (PRJCA):

1. **Download** Excel metadata via Selenium browser automation
2. **Parse** Excel files and extract all sheets
3. **Convert** to CSV format
4. **Add** source tracking fields

**Output:** `cncb_combined.csv`

### Stage 3: Final Merge (NCBI + CNCB)

Combine data from both sources:

1. **Normalize** column names across databases
   - `Accession` (CNCB) → `BioSample` (standard)
   - `Sample name` (CNCB) → `Sample_Name` (standard)
   - `Project accession` (CNCB) → `BioProject` (standard)
2. **Identify** common columns (BioProject, BioSample, Organism, etc.)
3. **Concatenate** vertically (combine rows)
4. **Add** source tracking (`Source_Database` field)

**Output:** `all_metadata_merged.csv`

## Output Files

### Individual Project Files

**NCBI Projects:**
- `PRJNA*_biosample.txt` - Raw BioSample data
- `PRJNA*_sra_runinfo.csv` - Raw SRA RunInfo

**CNCB Projects:**
- `CRA*.xlsx` - Raw Excel metadata (auto-converted from PRJCA)

### Combined Files

- `ncbi_merged_biosample_id.csv` - Merged NCBI data (SRA + BioSample)
- `cncb_combined.csv` - Combined CNCB data (converted from Excel)

### Final Merged File (Primary Output)

- **`all_metadata_merged.csv`** - Complete dataset with NCBI + CNCB data
  - Contains all records from both databases
  - Normalized column names
  - Source tracking field

## Data Structure

### Common Fields (Normalized)

Fields available in final merged output:

| Field | Description | Source |
|-------|-------------|--------|
| `BioProject` | BioProject accession | Both |
| `BioSample` | BioSample accession | Both |
| `Organism` | Source organism | Both |
| `Sample_Name` | Sample identifier | Both |
| `collection_date` | Collection date | Both |
| `geo_loc_name` | Geographic location | Both |
| `isolation_source` | Sample origin | Both |
| `Source_Database` | NCBI or CNCB/GSA | Auto-added |

### NCBI-Specific Fields

- `Run` - SRA run accession (e.g., SRR12345678)
- `LibraryStrategy` - Sequencing strategy (AMPLICON, WGS, etc.)
- `LibraryLayout` - SINGLE or PAIRED
- `Platform` - Sequencing platform (ILLUMINA, PACBIO, etc.)
- `spots`, `bases`, `avgLength` - Sequencing metrics
- `download_path` - FTP path to FASTQ files

### CNCB-Specific Fields

- `ID` - Internal sample ID
- `Accession` - BioSample accession (normalized to `BioSample`)
- `Sample title`, `Public description` - Descriptive fields
- Additional metadata fields vary by project

## Usage Examples

### Example 1: Mixed Database Projects

```bash
# Create input file with mixed BioProjects
cat > bioproject_ids/international_study.txt << EOF
PRJNA1335420
PRJCA004523
PRJNA1335413
EOF

# Run MetaDL
./bin/Meta2Data MetaDL \
    -i bioproject_ids/ \
    -o metadata_output/ \
    -e researcher@university.edu

# Check results
ls -lh metadata_output/
# all_metadata_merged.csv  <- Final output
# ncbi_merged_*.csv        <- NCBI data
# cncb_combined.csv        <- CNCB data
```

### Example 2: NCBI-Only Projects

```bash
# US + European projects
cat > bioproject_ids/ncbi_projects.txt << EOF
PRJNA123456
PRJEB234567
EOF

./bin/Meta2Data MetaDL \
    -i bioproject_ids/ \
    -o ncbi_output/ \
    -e your@email.com
```

### Example 3: CNCB-Only Projects

```bash
# Chinese database projects
cat > bioproject_ids/cncb_projects.txt << EOF
PRJCA004523
PRJCA012345
EOF

./bin/Meta2Data MetaDL \
    -i bioproject_ids/ \
    -o cncb_output/ \
    -e your@email.com
```

### Example 4: Large-Scale Download

```bash
# Process 100+ BioProjects from multiple databases
./bin/Meta2Data MetaDL \
    -i large_dataset/bioproject_list/ \
    -o large_dataset/metadata/ \
    -e researcher@university.edu

# Monitor progress (output is streamed)
tail -f large_dataset/metadata/download.log
```

## Performance

### Expected Download Times

| Dataset Size | NCBI Projects | CNCB Projects | Total Time |
|--------------|---------------|---------------|------------|
| Small (1-10) | 1-5 min | 2-10 min | ~5-15 min |
| Medium (10-50) | 10-30 min | 20-60 min | ~30-90 min |
| Large (50-100) | 30-60 min | 60-120 min | ~2-3 hours |

**Factors affecting speed:**
- Number of samples per BioProject
- NCBI/CNCB server response time
- Network bandwidth
- Rate limiting delays (0.4s between NCBI requests)

## Troubleshooting

### Error: "No .txt files found"

**Solution**: Ensure input directory contains `.txt` files

```bash
ls -l bioproject_ids/
# Should show .txt files
```

### Error: "Could not find download button" (CNCB)

**Cause**: CNCB website structure changed or BioProject doesn't exist

**Solutions**:
1. Verify BioProject ID on [CNCB GSA website](https://ngdc.cncb.ac.cn/gsa/)
2. Check if project is public (some are restricted)
3. Try manual download first to confirm accessibility

### Error: "No SRA data found"

**Cause**: BioProject has no associated SRA data (BioSample-only project)

**Solution**: This is expected for some projects. The pipeline will continue with BioSample data only.

### Slow CNCB Downloads

**Cause**: Selenium browser automation is slower than API

**Solutions**:
- Process CNCB projects in smaller batches
- Run during off-peak hours
- Use headless mode (set in code)

### Column Name Mismatches

**Issue**: Cannot find expected columns in merged file

**Solution**: Check column normalization mapping in code:

```python
# CNCB → Standard mapping
'Accession' → 'BioSample'
'Sample name' → 'Sample_Name'
'Project accession' → 'BioProject'
```

## Integration with AmpliconPIP

Use the merged metadata for downstream analysis:

```bash
# 1. Download metadata
./bin/Meta2Data MetaDL -i bioproject_ids/ -o metadata/ -e your@email.com

# 2. Use merged file with AmpliconPIP
./bin/Meta2Data AmpliconPIP \
    -m metadata/all_metadata_merged.csv \
    -o analysis_results/ \
    -t 8
```

**Note**: You may need to filter or format the merged CSV to match AmpliconPIP's expected columns.

## Advanced Usage

### Python API

Use the pipeline programmatically:

```python
from scripts.MetaDL.unified_metadata_downloader import run_unified_pipeline

results = run_unified_pipeline(
    input_folder='bioproject_ids/',
    output_folder='output/',
    email='your@email.com'
)

# Access results
bioproject_ids = results['bioproject_ids']
classified = results['classified']
ncbi_df = results['ncbi_df']
cncb_df = results['cncb_df']
final_df = results['final_df']

# Export custom format
final_df.to_excel('custom_output.xlsx', index=False)
```

### Custom Column Normalization

Modify normalization rules in `unified_metadata_downloader.py`:

```python
def normalize_column_names(df, database):
    if database == 'CNCB':
        rename_map = {
            'Accession': 'BioSample',
            'Sample name': 'Sample_Name',
            # Add your custom mappings here
            'Your_Field': 'Standard_Field',
        }
    return df
```

### Filtering Results

Filter merged data by criteria:

```bash
# Filter by organism
python3 -c "
import pandas as pd
df = pd.read_csv('metadata/all_metadata_merged.csv')
filtered = df[df['Organism'].str.contains('human', case=False)]
filtered.to_csv('human_samples.csv', index=False)
"
```

## File Structure

```
Meta2Data/
├── bin/
│   └── MetaDL                          # CLI entry point
├── scripts/
│   └── MetaDL/
│       ├── unified_metadata_downloader.py  # Main pipeline (NEW)
│       ├── ncbi_metadata_downloader.py     # NCBI-only (legacy)
│       └── dl_metadata_from_CNCB.py        # CNCB-only (standalone)
├── requirements_metadl.txt             # Combined dependencies
└── docs/
    └── UNIFIED-METADL-README.md        # This file
```

## Technical Details

### BioProject Classification

```python
PRJNA* → NCBI (USA)
PRJEB* → ENA (Europe) → Routed to NCBI Entrez
PRJDB* → DDBJ (Japan) → Routed to NCBI Entrez
PRJCA* → CNCB/GSA (China) → Selenium downloader
```

### Data Sources

- **NCBI Entrez API**: BioSample and SRA metadata
- **CNCB/GSA Website**: Excel-based metadata exports

### Merge Strategy Selection

The NCBI merge automatically selects the best strategy:

```python
if matched_biosample_only >= matched_composite:
    use_strategy = "BioSample ID only"
else:
    use_strategy = "BioProject + BioSample composite"
```

## Citation

If you use Unified MetaDL in your research, please cite:
- Meta2Data pipeline
- NCBI databases (BioProject, BioSample, SRA)
- CNCB/GSA database

## Support

For issues or questions:
1. Check this README
2. Review error messages in output
3. Verify input format and BioProject IDs
4. Check database website accessibility
5. Report issues with detailed error logs

---

**Last Updated**: 2025-12-10
**Version**: 2.0 (Unified NCBI + CNCB)
