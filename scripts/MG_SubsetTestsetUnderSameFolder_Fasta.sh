#!/usr/bin/env bash
# ############################################################################
# ##                                                                        ##
# ##     TRANSFERABLE TEST MODE MODULE FOR FASTA FILES                     ##
# ##                                                                        ##
# ##  Can be used as:                                                      ##
# ##    1. Inline function (copy-paste into your script)                   ##
# ##    2. Standalone script (source it or run separately)                 ##
# ##                                                                        ##
# ############################################################################

# ╔════════════════════════════════════════════════════════════════════════╗
# ║                      TEST SAMPLE EXTRACTION FUNCTION                   ║
# ╚════════════════════════════════════════════════════════════════════════╝

MG_SubsetTestsetUnderSameFolder_Fasta() {
  # Function that extracts test samples from FASTA files
  # Returns the path to test samples directory (via stdout)
  # 
  # Usage: 
  #   NEW_DIR=$(prepare_test_samples -i "$INPUT_DIR")
  #   NEW_DIR=$(prepare_test_samples -i "$INPUT_DIR" -s 5 -r 50)
  # 
  # Options:
  #   -i INPUT_DIR    Input directory containing FASTA files (required)
  #   -s NUM_SAMPLES  Number of samples to extract (default: 3)
  #   -r NUM_SEQS     Number of sequences per sample (default: 20)
  # 
  # NOTE: All messages are sent to stderr (>&2) so that only the directory
  #       path is returned via stdout. This allows clean variable capture.
  #       Test samples are created in parent directory of INPUT_DIR
  
  # Default values
  local INPUT_DIR=""
  local TEST_SAMPLE_COUNT=3           # Default: 3 samples
  local TEST_SEQ_PER_SAMPLE=20        # Default: 20 sequences per sample
  
  # Parse arguments
  local OPTIND
  while getopts "i:s:r:" opt; do
    case $opt in
      i) INPUT_DIR="$OPTARG" ;;
      s) TEST_SAMPLE_COUNT="$OPTARG" ;;
      r) TEST_SEQ_PER_SAMPLE="$OPTARG" ;;
      *)
        echo "ERROR: Invalid option" >&2
        echo "Usage: prepare_test_samples -i INPUT_DIR [-s NUM_SAMPLES] [-r NUM_SEQS]" >&2
        return 1
        ;;
    esac
  done
  
  # Validate required argument
  if [[ -z "$INPUT_DIR" ]]; then
    echo "ERROR: -i INPUT_DIR is required" >&2
    echo "Usage: prepare_test_samples -i INPUT_DIR [-s NUM_SAMPLES] [-r NUM_SEQS]" >&2
    return 1
  fi
  
  if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR" >&2
    return 1
  fi
  
  # Get parent directory of INPUT_DIR and create test_samples there
  local PARENT_DIR=$(dirname "$INPUT_DIR")
  local TEST_DIR="${PARENT_DIR}/test_samples"
  
  echo "" >&2
  echo "╔════════════════════════════════════════════════════════════════╗" >&2
  echo "║                                                                ║" >&2
  echo "║                    ⚠️  TEST MODE ENABLED ⚠️                     ║" >&2
  echo "║                                                                ║" >&2
  echo "║  Test Configuration:                                          ║" >&2
  echo "║    • Samples: ${TEST_SAMPLE_COUNT} sample(s)                                   ║" >&2
  echo "║    • Sequences per sample: ${TEST_SEQ_PER_SAMPLE} sequences                    ║" >&2
  echo "║                                                                ║" >&2
  echo "╚════════════════════════════════════════════════════════════════╝" >&2
  echo "" >&2
  
  # Create test directory
  mkdir -p "$TEST_DIR"
  
  echo "Preparing test samples..." >&2
  echo "  Input directory: $INPUT_DIR" >&2
  echo "  Output directory: $TEST_DIR" >&2
  echo "─────────────────────────────────────────────────────────────────" >&2
  
  local test_count=0
  
  for input_fasta in "${INPUT_DIR}"/*.fa "${INPUT_DIR}"/*.fasta; do
    [[ -e "$input_fasta" ]] || continue
    
    if [[ $test_count -lt $TEST_SAMPLE_COUNT ]]; then
      local filename=$(basename "$input_fasta")
      local output_fasta="${TEST_DIR}/${filename}"
      
      echo "" >&2
      echo "  [Sample $((test_count+1))/${TEST_SAMPLE_COUNT}] Processing: $filename" >&2
      
      # Count original sequences
      local original_count=$(grep -c "^>" "$input_fasta")
      echo "    Original sequences: ${original_count}" >&2
      
      # Extract first N sequences using awk
      awk -v n="$TEST_SEQ_PER_SAMPLE" '
        BEGIN { seq_count=0; printing=0 }
        /^>/ {
          if (seq_count < n) {
            printing=1
            seq_count++
            print
          } else {
            printing=0
          }
          next
        }
        printing { print }
      ' "$input_fasta" > "$output_fasta"
      
      local extracted_count=$(grep -c "^>" "$output_fasta")
      echo "    Extracted sequences: ${extracted_count}" >&2
      
      ((test_count++))
    else
      break
    fi
  done
  
  echo "" >&2
  echo "─────────────────────────────────────────────────────────────────" >&2
  
  if [[ $test_count -eq 0 ]]; then
    echo "❌ ERROR: No FASTA files found in ${INPUT_DIR}" >&2
    return 1
  fi
  
  echo "✓ Test sample preparation completed" >&2
  echo "  Total samples prepared: ${test_count}" >&2
  echo "  Test samples directory: ${TEST_DIR}" >&2
  echo "" >&2
  
  # CRITICAL: Return the test directory path to stdout
  # This allows the caller to reassign their variable
  echo "$TEST_DIR"
}

# ############################################################################
# ##                   END OF FUNCTION DEFINITION                         ##
# ############################################################################


# ============================================================================
# USAGE INSTRUCTIONS:
# ============================================================================
# 
# METHOD 1: INLINE USAGE (Copy function into your script)
# --------------------------------------------------------
# 
# #!/usr/bin/env bash
# 
# # 1. Define your input directory (can be ANY variable name)
# GENOME_DIR="/path/to/genomes"
# 
# # 2. Copy the prepare_test_samples function here
# 
# # 3. Set test mode
# TEST_MODE="true"
# 
# # 4. Conditionally replace your input directory
# if [[ "$TEST_MODE" == "true" ]]; then
#   # Default: 3 samples, 20 sequences each
#   GENOME_DIR=$(prepare_test_samples -i "$GENOME_DIR")
#   
#   # OR with custom values:
#   # GENOME_DIR=$(prepare_test_samples -i "$GENOME_DIR" -s 5 -r 50)
# fi
# 
# # 5. Use GENOME_DIR normally in your analysis
# for file in "${GENOME_DIR}"/*.fasta; do
#   echo "Processing: $file"
#   # Your analysis here
# done
# 
# ============================================================================
# 
# METHOD 2: STANDALONE SCRIPT (Source this file)
# -----------------------------------------------
# 
# #!/usr/bin/env bash
# 
# # Source the test mode module
# source /path/to/test_mode_module.sh
# 
# # Define your paths
# BIN_DIR="/data/genomes"
# 
# # Use the function with defaults
# TEST_MODE="true"
# if [[ "$TEST_MODE" == "true" ]]; then
#   BIN_DIR=$(prepare_test_samples -i "$BIN_DIR")
# fi
# 
# # Or with custom values
# # BIN_DIR=$(prepare_test_samples -i "$BIN_DIR" -s 10 -r 100)
# 
# # Your analysis
# for genome in "${BIN_DIR}"/*.fasta; do
#   run_analysis "$genome"
# done
# 
# ============================================================================
# 
# METHOD 3: STANDALONE EXECUTABLE
# --------------------------------
# 
# # Make this script executable
# chmod +x test_mode_module.sh
# 
# # Run it with default values (3 samples, 20 sequences)
# NEW_DIR=$(./test_mode_module.sh -i "/path/to/genomes")
# 
# # Run it with custom values
# NEW_DIR=$(./test_mode_module.sh -i "/path/to/genomes" -s 5 -r 50)
# 
# # Use the returned directory
# for file in "${NEW_DIR}"/*.fasta; do
#   run_analysis "$file"
# done
# 
# ============================================================================
# 
# IMPORTANT NOTES:
# ============================================================================
# 
# • Test samples are created in the PARENT directory of INPUT_DIR
#   Example: If INPUT_DIR="/data/genomes"
#            Test dir will be "/data/test_samples"
# 
# • Default values:
#   - Number of samples: 3
#   - Sequences per sample: 20
# 
# • Options:
#   -i : Input directory (REQUIRED)
#   -s : Number of samples (optional, default=3)
#   -r : Sequences per sample (optional, default=20)
# 
# ============================================================================

# If this script is executed directly (not sourced), run the function
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  if [[ $# -eq 0 ]]; then
    echo "Usage: $0 -i <input_directory> [-s <num_samples>] [-r <num_sequences>]" >&2
    echo "" >&2
    echo "This script prepares test samples from FASTA files:" >&2
    echo "" >&2
    echo "Options:" >&2
    echo "  -i INPUT_DIR    Input directory containing FASTA files (required)" >&2
    echo "  -s NUM_SAMPLES  Number of samples to extract (default: 3)" >&2
    echo "  -r NUM_SEQS     Number of sequences per sample (default: 20)" >&2
    echo "" >&2
    echo "Examples:" >&2
    echo "  $0 -i /data/genomes" >&2
    echo "  $0 -i /data/genomes -s 5 -r 50" >&2
    echo "" >&2
    echo "Returns: Path to test samples directory" >&2
    echo "Note: Test samples are created in parent directory of INPUT_DIR" >&2
    exit 1
  fi
  
  prepare_test_samples "$@"
fi