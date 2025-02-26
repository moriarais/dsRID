import sys

from scipy.fft import skip_backend

import re
import pysam
import argparse as ap
import pandas as pd
import numpy as np
from utils import *
import random

def main(args):
    print("extract_null.py init")
    
    # Argument parser for command-line options
    argp = ap.ArgumentParser(description="extract features for dsrna prediction from randomly sampled region",
                             formatter_class=ap.ArgumentDefaultsHelpFormatter)
    
    # Input BAM file argument
    argp.add_argument(
        "-b", "--bam_file",
        help="input bam file, with cs tags, sorted and indexed",
        type=str,
        default=""
    )

    # Output file argument
    argp.add_argument(
        "-o", "--out_file",
        help="output file location",
        type=str,
        default="./data/null.tsv"
    )

    # Read threshold argument
    argp.add_argument(
        "-r", "--read_threshold",
        help="number of reads required for prediction",
        type=int,
        default=6
    )
    
    # GTF file for splicing annotation
    argp.add_argument(
        "-s", "--splice_anno",
        help="gtf file for splicing region",
        type=str,
        default="./data/gencode.v34.annotation.sorted.gtf.gz"
    )

    # Chromosomes to be analyzed
    argp.add_argument(
        "-c", "--chr",
        help = "chromosomes to be analyzed",
        nargs='*',
        type=str,
        default = 'chr1'
    )

    args = argp.parse_args(args)

    mapq_thr = 20  # Minimum mapping quality threshold
    window = 2500  # Window size for region selection
    
    # Define chromosomes to analyze
    chrs = ["chrX", "chrY"]
    chrs.extend(["chr" + str(i) for i in range(1, 23)])
    print(chrs)
    
    sample_num = 10000  # Number of random samples to extract

    # Open BAM file for reading
    sam = pysam.AlignmentFile(args.bam_file, 'rb')
    feat_lst = list()

    # Get splicing annotation dictionary from GTF file
    sp_dic = get_gtf_splice_pos(args.splice_anno, chrs)

    while (len(feat_lst) < sample_num):
        # Randomly select a chromosome and exon site
        chr = random.choice(chrs)
        exonsite_ind = random.choice(range(len(sp_dic[chr])))
        
        try:
            int_start = sp_dic[chr][exonsite_ind][1]
            int_end = sp_dic[chr][exonsite_ind + 1][0]
        except:
            continue
        
        if int_start + 10 - window > int_end - 10 + window:
            continue
        
        # Select a random start position within the exon-intron boundary
        start = np.random.randint(int_start + 10 - window, int_end - 10 + window)
        end = start + window
        print(chr, start, end)
        
        # Count reads covering the selected region
        coverage = sam.count(chr, start, end)
        read_list = list()
        
        if coverage < args.read_threshold:
            continue
        
        # Iterate over reads in the selected region
        for read in sam.fetch(chr, start, end):
            if read.mapq < mapq_thr:
                continue
            elif read.is_secondary:  # Skip secondary alignments
                continue
            
            pos = read.reference_start  # Read start position
            read_start = read.reference_start
            read_end = read.reference_end
            
            # Convert CIGAR string to DataFrame
            cs = CIGAR_to_df(read.cigartuples, pos)
            cs_splice = cs.loc[cs['ope'] == 3]  # Select spliced reads
            
            if len(cs_splice) == 0:
                continue
            
            for ri, row in cs_splice.iterrows():
                low = int(row['low'])
                high = int(row['high'])
                
                # Extract sequences at splice junctions
                bp_start = read.query_sequence[low-read_start - 2 : low - read_start]
                bp_end = read.query_sequence[high-read_start : high - read_start + 2]
                length = int(row['val'])
                skipped_bases = read.query_sequence[low - read_start : high - read_start]
                
                # Store read splice site data
                read_list.append(
                    [read.query_name,
                    chr, low, high, length,
                    bp_start, bp_end, read_start, read_end, skipped_bases]
                )
        
        # Convert collected read data into a DataFrame
        mapped_splices = pd.DataFrame(read_list,
        columns=["read", "chr", "pos_start",
        "pos_end", "pos_len", "bp_start", "bp_end",
        "read_start", "read_end", "skipped_bases"])

        # Filter splices within the selected region
        sp_sites = mapped_splices.loc[(mapped_splices['pos_start'] > start - 10) &
      (mapped_splices['pos_end'] < end + 10) ]
        
        if len(sp_sites) == 0:
            continue

        print(sp_sites['pos_start'][sp_sites["pos_start"].isna()])
        
        # Convert read data to feature matrix
        feat_mat = reads_to_feature(mapped_splices, start, end, coverage)
        feat_mat['name'] = chr + ":" + str(start) + "-" + str(end)
        feat_mat['chr'] = chr
        feat_mat['start'] = start
        feat_mat['end'] = end
        print(feat_mat)
        
        feat_lst.append(feat_mat)
        
        # Periodically save results to file
        if len(feat_lst) % 1000 == 0:
            feat_total = pd.concat(feat_lst)
            feat_total.to_csv(args.out_file, sep='\t')
        print(len(feat_lst))
    
    # Save final collected features to output file
    feat_total = pd.concat(feat_lst)
    feat_total.to_csv(args.out_file, sep='\t')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
