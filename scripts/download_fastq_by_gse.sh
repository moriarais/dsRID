#!/bin/bash

# List of ENCODE experiment IDs
experiment_ids=(
    # "ENCSR462COR"
    # "ENCSR169YNI"
    # "ENCSR257YUB"
    # "ENCSR690QHM"
    # "ENCSR316ZTD"
    "ENCSR697ASE"
    "ENCSR094NFM"
    "ENCSR463IDK"
    "ENCSR205QMF"
)

# Base URL for ENCODE experiments
base_url="https://www.encodeproject.org/experiments"

# Loop through each experiment ID
for experiment_id in "${experiment_ids[@]}"; do
    echo "Processing $experiment_id..."

    # Fetch the experiment page and extract the GSE number
    gse_number=$(curl -s "${base_url}/${experiment_id}/" | grep -oP 'GEO:GSE[0-9]+' | head -1 | sed 's/GEO://')
    
    if [ -z "$gse_number" ]; then
        echo "No GSE number found for $experiment_id. Skipping..."
        continue
    fi

    echo "Found GSE number: $gse_number"

    # Create a directory for the GSE number
    mkdir -p "$gse_number"

    # Fetch FASTQ file links specifically under "Raw sequencing data" section
    fastq_links=$(curl -s "${base_url}/${experiment_id}/" | grep -oP '/ENCFF[0-9A-Z]+/@@download/[^\"]+\.fastq\.gz' | sed 's#^#https://www.encodeproject.org/files#')
    # https://www.encodeproject.org/files/ENCFF708BOP/@@download/ENCFF708BOP.fastq.gz

    echo "$fastq_links"

    # Download each FASTQ file into the GSE directory
    if [ -z "$fastq_links" ]; then
        echo "No FASTQ files found for $experiment_id. Skipping..."
        continue
    fi

    for link in $fastq_links; do
        echo "Downloading $link to raw_data/$gse_number/input/"
        nohup wget -P "$gse_number" "$link" &
    done
done

echo "Download complete. Files are organized by GSE numbers."
