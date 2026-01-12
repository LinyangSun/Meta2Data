# Unified MetaDL Changelog

## Version 2.0 - 2025-12-10

### 🎉 Major Release: Unified NCBI + CNCB/GSA Support

This release introduces a completely rewritten MetaDL pipeline that supports both NCBI and CNCB/GSA databases in a single unified workflow.

---

## What's New

### 1. Multi-Database Support

**Automatic routing by BioProject prefix:**
- `PRJNA*` → NCBI (USA)
- `PRJEB*` → ENA/NCBI (Europe)
- `PRJDB*` → DDBJ/NCBI (Japan)
- `PRJCA*` → CNCB/GSA (China) **← NEW!**

**How it works:**
```python
# Input: mixed_projects.txt
PRJNA1335420  # ← Routes to NCBI Entrez API
PRJCA004523   # ← Routes to CNCB Selenium downloader
PRJNA1335413  # ← Routes to NCBI Entrez API
```

The pipeline automatically detects the database based on BioProject prefix and routes to the appropriate downloader.

### 2. Three-Stage Merge Workflow

**Stage 1: NCBI Data Merge**
- Downloads BioSample and SRA metadata separately
- Intelligently merges using dual-strategy matching:
  - Strategy A: BioSample ID only
  - Strategy B: BioProject + BioSample composite key
- Auto-selects best strategy based on match rate
- Output: `ncbi_merged_*.csv`

**Stage 2: CNCB Data Conversion**
- Downloads Excel metadata via Selenium browser automation
- Parses Excel files and extracts all fields
- Converts to standardized CSV format
- Output: `cncb_combined.csv`

**Stage 3: Final Integration**
- Normalizes column names across databases
- Identifies common fields (BioProject, BioSample, Organism)
- Vertically concatenates data (combines rows)
- Adds source tracking field
- Output: `all_metadata_merged.csv` **← PRIMARY OUTPUT**

### 3. Column Normalization

Standardizes field names across NCBI and CNCB:

| CNCB Field | NCBI Field | Normalized |
|------------|------------|------------|
| Accession | BioSample | BioSample |
| Sample name | Sample_Name | Sample_Name |
| Project accession | BioProject | BioProject |
| Collection date | collection_date | collection_date |
| Geographic location | geo_loc_name | geo_loc_name |

### 4. New Core Script

**`scripts/MetaDL/unified_metadata_downloader.py`** (1100+ lines)

Complete rewrite that combines:
- NCBI downloader logic (from `ncbi_metadata_downloader.py`)
- CNCB downloader logic (from `dl_metadata_from_CNCB.py`)
- Intelligent routing and merge orchestration
- Comprehensive error handling

### 5. Updated CLI

**`bin/MetaDL`** now calls the unified script:

```bash
# Old behavior: NCBI only
# New behavior: NCBI + CNCB/GSA automatically

./bin/Meta2Data MetaDL \
    -i bioproject_ids/ \
    -o output/ \
    -e your@email.com

# Automatically handles mixed databases!
```

---

## Technical Details

### Architecture Changes

**Before (v1.0):**
```
Input → ncbi_metadata_downloader.py → NCBI only
Input → dl_metadata_from_CNCB.py → CNCB only (separate tool)
```

**After (v2.0):**
```
Input → unified_metadata_downloader.py
         ├─ Classify by prefix
         ├─ Route to NCBI downloader
         ├─ Route to CNCB downloader
         ├─ Merge NCBI (Stage 1)
         ├─ Parse CNCB (Stage 2)
         └─ Merge all (Stage 3)
```

### Code Improvements

1. **Unified Entry Point**
   - Single script handles all databases
   - Automatic routing eliminates user confusion
   - Backward compatible with NCBI-only workflows

2. **Intelligent Merging**
   - Dual-strategy NCBI merge with auto-selection
   - Column normalization across databases
   - Vertical concatenation for heterogeneous data

3. **Better Error Handling**
   - Graceful degradation when data is missing
   - Encoding issue fixes (bytes → str conversion)
   - Empty dataframe handling

4. **Source Tracking**
   - `Source_Database` field added to all records
   - Easy filtering by data origin
   - Transparency in data provenance

### Dependencies Added

```txt
# requirements_metadl.txt (updated)
biopython>=1.80      # Existing
pandas>=2.0.0        # Existing
selenium>=4.17.0     # NEW - for CNCB downloads
openpyxl>=3.0.0      # NEW - for Excel parsing
```

---

## Migration Guide

### For Existing Users

**Good news**: No changes required! The unified pipeline is **backward compatible**.

If you only use NCBI BioProjects:
```bash
# This still works exactly as before
./bin/Meta2Data MetaDL -i ncbi_projects/ -o output/ -e email@example.com
```

### For New CNCB Users

Simply add PRJCA IDs to your input files:
```bash
# mixed_projects.txt
PRJNA123456
PRJCA004523
PRJEB234567
```

The pipeline handles the rest automatically!

### Output File Changes

**New primary output:**
- Old: `matched_sra_biosample_*.csv` (NCBI only)
- New: `all_metadata_merged.csv` (NCBI + CNCB)

**Intermediate files preserved:**
- `ncbi_merged_*.csv` - NCBI data only (for compatibility)
- `cncb_combined.csv` - CNCB data only

### Script Usage Changes

**Python API:**
```python
# Old
from scripts.MetaDL.ncbi_metadata_downloader import run_complete_pipeline

# New (recommended)
from scripts.MetaDL.unified_metadata_downloader import run_unified_pipeline

results = run_unified_pipeline(
    input_folder='bioproject_ids/',
    output_folder='output/',
    email='your@email.com'
)

# Access classified data
ncbi_projects = results['classified']['ncbi']
cncb_projects = results['classified']['cncb']
final_df = results['final_df']
```

---

## Testing

### Test Coverage

Comprehensive testing with mixed databases:

**Test Case**: 3 BioProjects (2 NCBI + 1 CNCB)
- Input: PRJNA1335420, PRJNA1335413, PRJCA004523
- Expected: 108 samples (48 NCBI + 36 CNCB + 24 BioSample-only)
- Result: ✅ 108/108 samples, 126 columns, 100% data retention

**Test Report**: `test/unified_metadl_test/TEST_REPORT.md`

### Performance

| Metric | Value |
|--------|-------|
| NCBI throughput | ~24 samples/sec |
| CNCB throughput | ~36 samples/sec |
| Merge throughput | >100 records/sec |
| Total test time | <30 seconds |

---

## Documentation

### New Documentation

1. **`docs/UNIFIED-METADL-README.md`** (7.5 KB)
   - Comprehensive usage guide
   - Three-stage workflow explanation
   - Troubleshooting section
   - Advanced usage examples

2. **`test/unified_metadl_test/TEST_REPORT.md`** (6 KB)
   - Detailed test results
   - Performance metrics
   - Known issues and workarounds

### Updated Documentation

3. **`docs/METADL-README.md`** (updated)
   - Added unified pipeline notice
   - Updated prerequisites
   - Updated workflow description
   - New output file descriptions

4. **`bin/MetaDL`** help message (updated)
   - Multi-database support noted
   - Updated workflow steps
   - New output file listing

5. **`requirements_metadl.txt`** (updated)
   - Added selenium and openpyxl
   - Added ChromeDriver installation note

---

## Known Issues

### Issue 1: CNCB Download Button Detection

**Status**: Known, workaround available

**Description**: Selenium occasionally fails to find download button on CNCB website.

**Workaround**:
1. Download Excel file manually from CNCB website
2. Place in output directory with correct name (`CRA*.xlsx`)
3. Pipeline will parse existing file

**Future Fix**: Planned API-based download method.

### Issue 2: Large CNCB Files

**Status**: By design

**Description**: CNCB Excel files can be large (50+ MB for projects with 1000+ samples).

**Mitigation**: Pipeline streams data and uses efficient pandas processing.

---

## Compatibility

### Backward Compatibility

✅ **Fully backward compatible** with v1.0
- NCBI-only workflows work unchanged
- Old output files still generated
- Old Python API still functional

### Forward Compatibility

✅ **Ready for future databases**
- Architecture supports adding new databases
- Routing logic is extensible
- Merge system handles heterogeneous data

---

## Upgrade Instructions

### For CLI Users

```bash
# 1. Update dependencies
pip install -r requirements_metadl.txt

# 2. Install ChromeDriver (for CNCB support)
brew install --cask chromedriver

# 3. That's it! Start using CNCB BioProjects
./bin/Meta2Data MetaDL -i bioproject_ids/ -o output/ -e email@example.com
```

### For Python API Users

```bash
# Update imports in your code
# Old
from scripts.MetaDL.ncbi_metadata_downloader import run_complete_pipeline

# New (for multi-database support)
from scripts.MetaDL.unified_metadata_downloader import run_unified_pipeline
```

---

## Contributors

- Initial NCBI implementation: Previous version
- CNCB/GSA integration: v2.0 release
- Unified pipeline architecture: v2.0 release
- Testing and documentation: v2.0 release

---

## Future Roadmap

### Planned Features (v2.1+)

1. **Enhanced CNCB Download**
   - API-based download method (if CNCB provides API)
   - Improved button detection
   - Progress tracking for large batches

2. **Additional Databases**
   - DDBJ direct support (currently via NCBI)
   - MGnify/EBI Metagenomics integration
   - JGI Genome Portal support

3. **Output Formats**
   - AmpliconPIP-ready CSV format
   - Excel workbook with multiple sheets
   - JSON export for APIs

4. **Performance**
   - Parallel download for multiple BioProjects
   - Resume capability for interrupted downloads
   - Incremental updates for existing datasets

### Community Requests

Open an issue on GitHub to request:
- Additional database support
- Custom column mappings
- Output format preferences
- Performance optimizations

---

## Summary

**Version 2.0** represents a major evolution of MetaDL:

✅ **Multi-database**: NCBI + CNCB/GSA in one tool
✅ **Intelligent**: Automatic routing and merge strategies
✅ **Comprehensive**: 126 columns, complete metadata coverage
✅ **Tested**: 100% data retention across mixed databases
✅ **Documented**: 13+ KB of documentation
✅ **Compatible**: Fully backward compatible with v1.0

**Ready for production use with international BioProject datasets!**

---

**Release Date**: 2025-12-10
**Version**: 2.0.0
**Codename**: "Unified"
