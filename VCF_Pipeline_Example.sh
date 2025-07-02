#!/bin/bash

# VCF Processing Pipeline for Genomic Variant Analysis
# Author: Lavanyaa Gupta
# Description: Automated preprocessing of VCF files for cancer genomics analysis

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${SCRIPT_DIR}/input_vcf"
OUTPUT_DIR="${SCRIPT_DIR}/processed_vcf"
LOG_FILE="${SCRIPT_DIR}/processing.log"
TEMP_DIR="${SCRIPT_DIR}/temp"

# Create directories as necessary
mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

# Logging function to create timestamped log entries
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Check dependencies
check_dependencies() {
    log "Checking dependencies..."
    
    local deps=("bcftools" "tabix" "bgzip")
    for dep in "${deps[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            log "ERROR: $dep is not installed or not in PATH"
            exit 1
        fi
    done
    
    log "All dependencies found"
}

# Function to validate VCF file
validate_vcf() {
    local vcf_file="$1"
    log "Validating VCF file: $(basename "$vcf_file")"
    
    if [[ ! -f "$vcf_file" ]]; then
        log "ERROR: VCF file not found: $vcf_file"
        return 1
    fi
    
    # Check if file is compressed
    if [[ "$vcf_file" == *.gz ]]; then
        if ! bgzip -t "$vcf_file" 2>/dev/null; then
            log "ERROR: Invalid compressed VCF file: $vcf_file"
            return 1
        fi
    fi
    
    # Basic VCF format validation
    local header_check
    if [[ "$vcf_file" == *.gz ]]; then
        header_check=$(zcat "$vcf_file" | head -1)
    else
        header_check=$(head -1 "$vcf_file")
    fi
    
    if [[ ! "$header_check" =~ ^##fileformat=VCF ]]; then
        log "ERROR: Invalid VCF header in file: $vcf_file"
        return 1
    fi
    
    log "VCF validation passed: $(basename "$vcf_file")"
    return 0
}

# Function to perform quality filtering
quality_filter() {
    local input_vcf="$1"
    local output_vcf="$2"
    local min_qual="${3:-30}"
    local min_depth="${4:-10}"
    
    log "Applying quality filters (QUAL >= $min_qual, DP >= $min_depth)"
    
    bcftools filter \
        -i "QUAL >= $min_qual && INFO/DP >= $min_depth" \
        -O z \
        -o "$output_vcf" \
        "$input_vcf"
    
    # Index the filtered VCF
    tabix -p vcf "$output_vcf"
    
    # Count variants before and after filtering
    local original_count
    local filtered_count
    
    if [[ "$input_vcf" == *.gz ]]; then
        original_count=$(zcat "$input_vcf" | grep -v "^#" | wc -l)
    else
        original_count=$(grep -v "^#" "$input_vcf" | wc -l)
    fi
    
    filtered_count=$(zcat "$output_vcf" | grep -v "^#" | wc -l)
    
    log "Quality filtering complete: $original_count -> $filtered_count variants"
}

# Function to normalize variants
normalize_variants() {
    local input_vcf="$1"
    local output_vcf="$2"
    local reference_fasta="$3"
    
    log "Normalizing variants with reference: $(basename "$reference_fasta")"
    
    bcftools norm \
        -f "$reference_fasta" \
        -m -both \
        -O z \
        -o "$output_vcf" \
        "$input_vcf"
    
    tabix -p vcf "$output_vcf"
    log "Variant normalization complete"
}

# Function to annotate variants
annotate_variants() {
    local input_vcf="$1"
    local output_vcf="$2"
    
    log "Adding basic annotations"
    
    # Add INFO fields for analysis
    bcftools annotate \
        -a "$input_vcf" \
        -h <(echo '##INFO=<ID=VARIANT_TYPE,Number=1,Type=String,Description="Type of variant (SNP, INDEL, etc.)">') \
        -O z \
        -o "$output_vcf" \
        "$input_vcf"
    
    tabix -p vcf "$output_vcf"
    log "Variant annotation complete"
}

# Function to generate summary statistics
generate_stats() {
    local vcf_file="$1"
    local stats_file="$2"
    
    log "Generating summary statistics"
    
    bcftools stats "$vcf_file" > "$stats_file"
    
    # Extract key statistics
    local total_variants
    local snps
    local indels
    
    total_variants=$(grep "^SN" "$stats_file" | grep "number of records:" | cut -f4)
    snps=$(grep "^SN" "$stats_file" | grep "number of SNPs:" | cut -f4)
    indels=$(grep "^SN" "$stats_file" | grep "number of indels:" | cut -f4)
    
    log "Statistics generated - Total: $total_variants, SNPs: $snps, INDELs: $indels"
}

# Function to process a single VCF file
process_vcf() {
    local input_vcf="$1"
    local base_name
    base_name=$(basename "$input_vcf" .vcf)
    base_name=$(basename "$base_name" .vcf.gz)
    
    log "Processing VCF file: $(basename "$input_vcf")"
    
    # Validate input
    if ! validate_vcf "$input_vcf"; then
        log "ERROR: Validation failed for $input_vcf"
        return 1
    fi
    
    # Define intermediate and output files
    local filtered_vcf="${TEMP_DIR}/${base_name}_filtered.vcf.gz"
    local normalized_vcf="${TEMP_DIR}/${base_name}_normalized.vcf.gz"
    local final_vcf="${OUTPUT_DIR}/${base_name}_processed.vcf.gz"
    local stats_file="${OUTPUT_DIR}/${base_name}_stats.txt"
    
    # Step 1: quality filtering
    quality_filter "$input_vcf" "$filtered_vcf"
    
    # Step 2: normalization (if reference is available)
    if [[ -f "${SCRIPT_DIR}/reference.fa" ]]; then
        normalize_variants "$filtered_vcf" "$normalized_vcf" "${SCRIPT_DIR}/reference.fa"
        local processing_vcf="$normalized_vcf"
    else
        log "WARNING: Reference genome not found, skipping normalization"
        local processing_vcf="$filtered_vcf"
    fi
    
    # Step 3: annotation
    annotate_variants "$processing_vcf" "$final_vcf"
    
    # Step 4: generate statistics
    generate_stats "$final_vcf" "$stats_file"
    
    log "Successfully processed: $(basename "$input_vcf") -> $(basename "$final_vcf")"
}

# Function to process all VCF files in input directory
process_all_vcfs() {
    log "Starting batch processing of VCF files"
    
    local vcf_files=("$INPUT_DIR"/*.vcf "$INPUT_DIR"/*.vcf.gz)
    local processed_count=0
    local failed_count=0
    
    for vcf_file in "${vcf_files[@]}"; do
        if [[ -f "$vcf_file" ]]; then
            if process_vcf "$vcf_file"; then
                ((processed_count++))
            else
                ((failed_count++))
                log "ERROR: Failed to process $vcf_file"
            fi
        fi
    done
    
    log "Batch processing complete: $processed_count successful, $failed_count failed"
}

# Function to generate processing report
generate_report() {
    local report_file="${OUTPUT_DIR}/processing_report.txt"
    
    log "Generating processing report"
    
    cat > "$report_file" << EOF


VCF PROCESSING PIPELINE REPORT
==============================
Processing Date: $(date)
Script Version: 1.0
Author: Lavanyaa Gupta

CONFIGURATION:
- Input Directory: $INPUT_DIR
- Output Directory: $OUTPUT_DIR
- Quality Threshold: >= 30
- Depth Threshold: >= 10

PROCESSED FILES:
EOF
    
    # List all processed files with their statistics
    for stats_file in "$OUTPUT_DIR"/*_stats.txt; do
        if [[ -f "$stats_file" ]]; then
            local base_name
            base_name=$(basename "$stats_file" _stats.txt)
            
            echo "" >> "$report_file"
            echo "File: ${base_name}" >> "$report_file"
            echo "----------------------------------------" >> "$report_file"
            
            # Extract key statistics
            grep "^SN" "$stats_file" | head -10 >> "$report_file"
        fi
    done
    
    cat >> "$report_file" << EOF

QUALITY METRICS:
- All variants passed minimum quality thresholds
- Normalization applied where reference available
- Basic annotations added for downstream analysis

NEXT STEPS:
1. Review individual file statistics
2. Proceed with Python analysis pipeline
3. Generate visualizations and clinical reports

For questions or issues, contact: lgupta8@asu.edu
EOF
    
    log "Processing report generated: $report_file"
}

# Function to clean up temporary files
cleanup() {
    log "Cleaning up temporary files"
    rm -rf "$TEMP_DIR"
    log "Cleanup complete"
}

# Main execution function
main() {
    log "Starting VCF Processing Pipeline"
    log "Script Directory: $SCRIPT_DIR"
    
    # Check if running with proper permissions
    if [[ ! -w "$SCRIPT_DIR" ]]; then
        log "ERROR: No write permission in script directory"
        exit 1
    fi
    
    # Check dependencies
    check_dependencies
    
    # Process command line arguments
    case "${1:-all}" in
        "single")
            if [[ -z "${2:-}" ]]; then
                log "ERROR: Please provide VCF file path for single file processing"
                exit 1
            fi
            process_vcf "$2"
            ;;
        "all"|"batch")
            process_all_vcfs
            ;;
        "help"|"-h"|"--help")
            cat << EOF

            