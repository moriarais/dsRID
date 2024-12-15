#!/bin/bash

# Function to recursively unzip all .fastq.gz files in a given directory
unzip_fastq_files() {
  local base_directory=$1

  # Find all .fastq.gz files under the specified directory
  find "$base_directory" -type f -name "*.fastq.gz" | while read -r file; do
    echo "Decompressing: $file"
    nohup gzip -dk -f "$file" &
  done
}

# Check if a directory argument is provided
if [ $# -eq 0 ]; then
  echo "Usage: $0 <base_directory>"
  exit 1
fi

# Call the function with the provided base directory
unzip_fastq_files "$1"