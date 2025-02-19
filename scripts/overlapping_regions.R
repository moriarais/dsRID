library(ggplot2)
library(dplyr)

# Load BED file (adjust column names as needed)
overlap_data <- read.table("/private5/Projects/raismor/dsRID_project/external_data/overlapping_regions_alu_novel_dsRID.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Assign column names based on BED format
colnames(overlap_data) <- c("chr_A", "start_A", "end_A", "chr_B", "start_B", "end_B")

# Optional: Add an index column for visualization
overlap_data <- overlap_data %>% mutate(id = row_number())

# Print first few rows
head(overlap_data)

ggplot(overlap_data) +
  geom_segment(aes(x = start_A, xend = end_A, y = id, yend = id), color = "blue", size = 2) +  # Regions from file1
  geom_segment(aes(x = start_B, xend = end_B, y = id + 0.3, yend = id + 0.3), color = "red", size = 2) +  # Regions from file2
  labs(title = "Overlapping Genomic Regions", x = "Genomic Position", y = "Region ID") +
  theme_minimal()

