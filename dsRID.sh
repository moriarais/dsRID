# !bin/bash

# Save arguments into variables
input_file=$1
output_file=$2

python ./src/extract_train.py -b "$input_file" -o "$output_file"/dsRID_train.tsv
python ./src/extract_null.py -b "$input_file" -o "$output_file"/dsRID_null.tsv

python ./src/concat_train_null.py -t "$output_file"/dsRID_train.tsv \
 -n "$output_file"/dsRID_null.tsv -o "$output_file"/dsRID_data.tsv

python ./src/extract_whole.py -b "$input_file" -o "$output_file"/dsRID_whole.tsv

python ./src/model_predict.py "$output_file"/dsRID_data.tsv "$output_file"
