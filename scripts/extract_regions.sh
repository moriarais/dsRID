region_bed=$1
outdir=${2:="/private9/Projects/dsRNAProject/BamLinksWithIndexPooled"}
name=$(basename ${region_bed/bed/})
name=${name/txt/}
fixed_region_bed=${region_bed/.bed/.fixed.bed}
fixed_region_bed=${fixed_region_bed/txt/fixed.bed}
# fix if needed
cat $region_bed | tr "[:\-]" "\t" > $fixed_region_bed
# execute
time samtools-1.19 merge --threads 40 -b /private9/Projects/dsRNAProject/BamLinksWithIndexPooled/AllTissues_AlignedHyperEditedReadsPE_0.05_0.6_30_0.6_0.1_0.8_0.2.sortedByCoord.out.files.txt -L $fixed_region_bed -o "${outdir}/AllTissues_${name}_AlignedHyperEditedReadsPE_0.05_0.6_30_0.6_0.1_0.8_0.2.sortedByCoord.out.bam##idx##${outdir}/AllTissues_${name}_AlignedHyperEditedReadsPE_0.05_0.6_30_0.6_0.1_0.8_0.2.sortedByCoord.out.bam.bai" --write-index
# remove
rm $fixed_region_bed
