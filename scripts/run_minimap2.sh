#!/bin/bash

# List of folders
folder_list=("GSE219452" "GSE219545" "GSE219574" "GSE219641" "GSE219645" "GSE219655" "GSE219675" "GSE219875")

# Path to the reference genome
reference_genome="/private/dropbox/Genomes/Human/hg38/hg38.fa"

# Loop through each folder
for folder in "${folder_list[@]}"; do
  # Define the input directory path
  input_directory="${folder}/input"

  # Find all .fastq files in the input directory
  find "$input_directory" -type f -name "*.fastq" | while read -r fastq_file; do
    # Extract the base name (without directory and extension)
    base_name=$(basename "$fastq_file" .fastq)

    # Define the output SAM file name
    output_sam="${folder}/input/${base_name}.sam"

    # Run the minimap2 command in the background
    echo "Running: minimap2 -ax splice -t 8 --cs -L $reference_genome $fastq_file > $output_sam &"
    nohup minimap2 -ax splice -t 8 --cs -L "$reference_genome" "$fastq_file" > "$output_sam" &
  done
done
