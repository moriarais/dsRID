 #!/bin/bash

# List of folders
# folder_list=( "GSE219641" "GSE219645" "GSE219871")
folder_list=( "GSE219574" "GSE219641" "GSE219645" "GSE219655" "GSE219871" "GSE219875")

# https://www.encodeproject.org/files/ENCFF708BOP/@@download/ENCFF708BOP.fastq.gz

# Path to the reference genome
reference_genome="/private/dropbox/Genomes/Human/hg38/hg38.fa"
echo "$folder_list"
# Loop through each folder
for folder in "${folder_list[@]}"; do
  # Define the input directory path
  input_directory="/private11/Projects/raismor/dsRID_project/raw_data/${folder}/input"
  echo "$input_directory"
  
  #   # Find all .fastq files in the input directory
  # find "$input_directory" -type f -name "*.fastq" | while read -r fastq_file; do
  #   # Extract the base name (without directory and extension)
  #   base_name=$(basename "$fastq_file" .fastq)

  #   # Define the output SAM file name
  #   output_sam="/private11/Projects/raismor/dsRID_project/raw_data/${folder}/input/${base_name}.sam"

  #   # Run the minimap2 command in the background
  #   echo "Running: minimap2 -ax splice -t 8 --cs -L $reference_genome $fastq_file > $output_sam &"
  #   minimap2 -ax splice --cs -t 8 -L  "$reference_genome" "$fastq_file" > "$output_sam" &
  # done

  # Find all .sam files in the input directory
  # find "$input_directory" -type f -name "*.sam" | while read -r sam_file; do
  #   # Extract the base name (without directory and extension)
  #   base_name=$(basename "$sam_file" .sam)
  #   echo "$base_name"
  #   # Define the output BAM file name
  #   output_bam="/private11/Projects/raismor/dsRID_project/raw_data/${folder}/input/${base_name}.bam"

  #   # Run the samtools command to convert SAM to BAM
  #   echo "Running: samtools view -b $sam_file > $output_bam"
  #   nohup samtools view -b "$sam_file" > "$output_bam" &
  # done

  # Find all .bam files in the input directory
  # find "$input_directory" -type f -name "*.bam" | while read -r bam_file; do
  #   # Extract the base name (without directory and extension)
  #   base_name=$(basename "$bam_file" .bam)
  #   echo "$base_name"
  #   # Define the output sorted BAM file name
  #   output_sorted_bam="/private11/Projects/raismor/dsRID_project/raw_data/${folder}/input/${base_name}_sorted.bam"

  #   # Run the samtools sort command
  #   echo "Running: samtools sort $bam_file -o $output_sorted_bam"
  #   nohup samtools sort "$bam_file" -o "$output_sorted_bam" &
  # done

  # find "$input_directory" -type f -name "*_sorted.bam" -print | while read -r sorted_bam_file; do

  #   # Extract the base name (without directory and extension)
  #   base_name=$(basename "$sorted_bam_file" _sorted.bam)

  #   echo "$base_name"

  #   # Define the output sorted BAM file name
  #   output_sorted_bam="${sorted_bam_file}"

  #   # Run the samtools index command
  #   echo "Running: samtools index $output_sorted_bam"
  #   nohup samtools index "$output_sorted_bam" &
  # done


    find "$input_directory" -type f -name "*_sorted.bam" | while read -r sorted_bam_file; do
      # Extract the base name (without directory and extension)
      base_name=$(basename "$sorted_bam_file" _sorted.bam)
      output_directory="/private11/Projects/raismor/dsRID_project/raw_data/${folder}/output"

      # Run the samtools index command
      echo "Running: bash dsRID.sh $sorted_bam_file"
      nohup bash /private11/Projects/raismor/dsRID_project/dsRID/dsRID.sh ${sorted_bam_file} ${output_directory} &

    done
done

