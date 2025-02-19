# !bin/bash

# Save arguments into variables
input_file=$1
output_file=$2

# Check if both arguments are provided
if [ -z "$input_file" ] || [ -z "$output_file" ]; then
  echo "Usage: bash dsRID.sh <input_file> <output_directory>"
  exit 1
fi

# Check if the BAM file exists
if [ ! -f "$input_file" ]; then
  echo "Error: BAM file '$input_file' does not exist."
  exit 1
fi

# Check if the output directory exists, and create it if it doesn't
if [ ! -d "$output_file" ]; then
  echo "Output directory '$output_file' does not exist. Creating it now..."
  # mkdir -p "$output_file"  # -p creates any missing parent directories
fi

# Define the output file for all errors and warnings
ERROR_LOG="../errors_and_warnings.txt"

# Clear previous content in the error log (optional)
> "$ERROR_LOG"

echo "----- Errors and Warnings for extract_train.py -----" >> "$ERROR_LOG"
python ./src/extract_train.py -b "$input_file" -o "$output_file"/dsRID_train.tsv 2>> "$ERROR_LOG"

echo "----- Errors and Warnings for extract_null.py -----" >> "$ERROR_LOG"
python ./src/extract_null.py -b "$input_file" -o "$output_file"/dsRID_null.tsv 2>> "$ERROR_LOG"

echo "----- Errors and Warnings for concat_train_null.py -----" >> "$ERROR_LOG"
python ./src/concat_train_null.py -t "$output_file"/dsRID_train.tsv \
 -n "$output_file"/dsRID_null.tsv -o "$output_file"/dsRID_data.tsv 2>> "$ERROR_LOG"

echo "----- Errors and Warnings for extract_whole.py -----" >> "$ERROR_LOG"
python ./src/extract_whole.py -b "$input_file" -o "$output_file"/dsRID_whole.tsv 2>> "$ERROR_LOG"

echo "----- Errors and Warnings for model_predict.py -----" >> "$ERROR_LOG"
python ./src/model_predict.py "$output_file"/dsRID_data.tsv "$output_file" 2>> "$ERROR_LOG"
