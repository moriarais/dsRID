import matplotlib.pyplot as plt  
import pandas as pd  

# Load BED file correctly  
bed = pd.read_csv("/private5/Projects/raismor/dsRID_project/external_data/overlapping_regions_alu_novel_dsRID_nooverlap.bed",
                   delim_whitespace=True, header=None, 
                  names=["chr1", "start1", "end1", "chr2", "start2", "end2", "overlap"])

# Convert start and end positions to integers
bed[["start1", "end1", "start2", "end2"]] = bed[["start1", "end1", "start2", "end2"]].astype(int)

# Create plot
fig, ax = plt.subplots(figsize=(10, 5))

for i, row in bed.iterrows():
    ax.plot([row["start1"], row["end1"]], [i, i], lw=5, color='blue', label="alu" if i == 0 else "")
    ax.plot([row["start2"], row["end2"]], [i, i], lw=5, color='red', label="novel_dsRID" if i == 0 else "")

ax.set_xlabel("Genomic Position")
ax.set_ylabel("Overlapping Regions")
ax.set_yticks([])
ax.legend()

# Save the plot
output_path = "/private5/Projects/raismor/dsRID_project/plots/overlapping_regions_final.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")  # Save as high-resolution image
plt.close()  # Close the figure to free memory

print(f"Plot saved as {output_path}")
