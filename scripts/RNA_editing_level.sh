#!/bin/bash

# RNA A-to-I Editing Detection Pipeline for PacBio Long Reads
# Requires: samtools, bedtools, Python (for downstream analysis)

# Input files
BED_FILE="/private5/Projects/raismor/dsRID_project/dsRID/data/dsRID_novel_dsRNA_merged.bed"  # BED file with regions of interest
BAM_FILE="/private5/Projects/raismor/dsRID_project/raw_data/GSE219452/input/ENCFF260AWP_sorted.bam"    # Aligned PacBio long reads (sorted and indexed)
REFERENCE="/private/dropbox/Genomes/Human/hg38/hg38.fa"  # Reference genome (hg38)
OUTPUT_DIR="/private5/Projects/raismor/dsRID_project/raw_data/GSE219452/results"         # Output directory

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

echo "Step 1: Extract reads overlapping BED regions..."
bedtools intersect -a $BAM_FILE -b $BED_FILE > $OUTPUT_DIR/filtered_reads.bam

# Index the filtered BAM file
echo "Indexing filtered BAM file..."
samtools index $OUTPUT_DIR/filtered_reads.bam

echo "Step 2: Generate pileup file..."
samtools mpileup -f $REFERENCE $OUTPUT_DIR/filtered_reads.bam > $OUTPUT_DIR/pileup.txt

# Python script to identify A-to-G mismatches and calculate editing levels
PYTHON_SCRIPT=$OUTPUT_DIR/identify_a_to_g_edits.py
cat << 'EOF' > $PYTHON_SCRIPT
import pandas as pd

# Load pileup file
pileup = pd.read_csv("/private5/Projects/raismor/dsRID_project/raw_data/GSE219452/results/pileup.txt", sep='\t', header=None, 
                     names=['chrom', 'pos', 'ref_base', 'depth', 'read_bases'])

# Filter for A-to-G mismatches
a_to_g_edits = pileup[
    (pileup['ref_base'] == 'A') &  # Reference is A
    (pileup['read_bases'].str.contains('G'))  # Reads contain G
]

# Calculate editing levels
def calculate_editing_level(read_bases):
    total_reads = len(read_bases)
    g_reads = read_bases.count('G')
    return g_reads / total_reads if total_reads > 0 else 0

# Apply the function to calculate editing levels
a_to_g_edits['editing_level'] = a_to_g_edits['read_bases'].apply(calculate_editing_level)

# Save results
a_to_g_edits.to_csv("/private5/Projects/raismor/dsRID_project/raw_data/GSE219452/results/a_to_g_edits.csv", index=False)
EOF

# Run the Python script
echo "Step 3: Process pileup file to identify A-to-G edits..."
python3 $PYTHON_SCRIPT

# Final output files
echo "Results saved to:"
echo "$OUTPUT_DIR/filtered_reads.bam"
echo "$OUTPUT_DIR/pileup.txt"
echo "$OUTPUT_DIR/a_to_g_edits.csv"

# Optional: Cleanup intermediate files
# Uncomment the line below if you want to remove intermediate files after analysis
# rm $OUTPUT_DIR/filtered_reads.bam $OUTPUT_DIR/filtered_reads.bam.bai $OUTPUT_DIR/pileup.txt

echo "Pipeline completed successfully!"
