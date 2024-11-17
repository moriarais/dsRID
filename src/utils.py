import re
from statistics import mean
import pysam
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
import warnings
warnings.simplefilter(action="ignore", category=pd.errors.SettingWithCopyWarning)

# Split the 'cs' string into a list of tuples containing operation and value parts
def split_cs_string(cs_string):
    return list(
        zip(
            re.sub('[0-9a-z]', ' ', cs_string).split(),  # Replace numbers and lowercase letters with spaces
            re.sub('[:*\-+~]', ' ', cs_string).split()   # Replace operation characters with spaces
        )
    )

# Define functions to interpret lengths based on 'cs' string operations
cslenfuncs = {
    ':': int,                       # Matches
    '*': lambda x: 1,               # Mismatches
    '+': lambda x: 0,               # Insertions
    '-': len,                       # Deletions
    '~': lambda x: int(re.sub('[a-z]', '', x))  # Large skips, extract numeric value
}

# Define functions to interpret lengths based on 'CIGAR' operations
cigarlenfuncs = {
    0 : lambda x: x,   # Match/Mismatch
    1 : lambda x: 0,   # Insertion
    2 : lambda x: x,   # Deletion
    3 : lambda x: x,   # Skipped region
    4 : lambda x: x,   # Soft clipping
    5 : lambda x: x,   # Hard clipping
    6 : lambda x: x,   # Padding
    7 : lambda x: x,   # Match (Equal)
    8 : lambda x: x    # Mismatch (Different)
}

# Converts a 'cs' string to a DataFrame representing genome positions and operations
def cs_to_df(cs_string, pos):
    cs = split_cs_string(cs_string)
    cslist = list()
    for a, b in cs:
        low = pos
        pos += cslenfuncs[a](b)  # Update position based on operation type
        high = pos
        cslist.append([low, high, a, b])  # Add start, end, operation, and value
    csdf = pd.DataFrame(
        np.row_stack(cslist),
        columns=['low', 'high', 'ope', 'val']
    )
    csdf.loc[:, 'low'] = csdf['low'].astype(int)
    csdf.loc[:, 'high'] = csdf['high'].astype(int)
    return csdf

# Converts a CIGAR tuple list to a DataFrame
def CIGAR_to_df(cigar_tuples, pos):
    cigar_list = list()
    for symbol, length in cigar_tuples:
        low = pos
        pos += cigarlenfuncs[symbol](length)  # Update position based on operation
        high = pos
        cigar_list.append([low, high, symbol, length])  # Add start, end, operation, and value
    cigar_df = pd.DataFrame(
        np.row_stack(cigar_list),
        columns=['low', 'high', 'ope', 'val']
    )
    return cigar_df

# Processes mapped splice sites, calculates various features, and returns results as a DataFrame
def reads_to_feature(mapped_splices, start, end, coverage):
    sp_sites = mapped_splices.loc[(mapped_splices['pos_start'] > start - 10) &
      (mapped_splices['pos_end'] < end + 10) ]

    group_ids, group_num, group_std = get_sp_cluster(sp_sites)

    # Calculate mean and variance of splice site positions
    mean_sp_start = np.mean(sp_sites['pos_start'])
    mean_sp_end = np.mean(sp_sites['pos_end'])
    var_sp_start = np.std(sp_sites['pos_start'] - mean_sp_start)
    var_sp_end = np.std(sp_sites['pos_end'] - mean_sp_end)

    res = pd.DataFrame(data ={"std_start" : var_sp_start,
     "std_end": var_sp_end, "mean_start":mean_sp_start,
     "mean_end": mean_sp_end, "len_skip":mean_sp_end - mean_sp_start,
      "coverage" : coverage, "num_skip": sp_sites.shape[0],
      "skip_ratio": sp_sites.shape[0] / coverage,
       "group_num": group_num, "group_std" : group_std 
       }, index=[0], 
       )
    sp_sites['group'] = group_ids

    #print(sp_sites['skipped_bases'])

    gc_skip = np.mean(sp_sites["skipped_bases"].apply(lambda st : len(re.findall("[GC]", st)) / len(st) if len(st) > 0 else 0 ))

    res['gc_skip'] = gc_skip

    # Count occurrences of boundary positions
    count_bp_start = count_bp(sp_sites['bp_start'], True)
    count_bp_end = count_bp(sp_sites['bp_end'], False)

    # Concatenate boundary counts to the result
    res = pd.concat([res, count_bp_start, count_bp_end], axis=1)
    return res

# Groups splice sites based on proximity, returns group IDs and statistics
def get_sp_cluster(sp_sites, window_size=50):
    "Given splicing sites return groups of sp sites that are at least 100 bps apart"
    group_id = 0
    pos_temp = min(sp_sites['pos_start'])
    group_ids = [0]
    group_dic = {group_id: [pos_temp]}
    for pos_start in sp_sites['pos_start'][1:].sort_values():
        if pos_start - pos_temp > window_size:  # Create a new group if distance exceeds window size
            group_id += 1
            group_dic[group_id] = [pos_start]
            pos_temp = mean(group_dic[group_id])
        else:
            group_dic[group_id].append(pos_start)
            pos_temp = mean(group_dic[group_id])
        group_ids.append(group_id)
    
    stds = []
    # Calculate standard deviations of each group's positions
    for key in group_dic.keys():
        std_temp = np.std(group_dic[key])
        stds.append(std_temp)
    return group_ids, max(group_ids) + 1, mean(stds)

# Counts the frequency of specific boundary positions
def count_bp(bps, start):
    res = {}
    prefix = 'bp_start_' if start else 'bp_end_'
    for bp in bps.unique():
        res[prefix + bp] = 0
    for bp in bps:
        res[prefix + bp] += 1 / len(bps)  # Increment frequency count
    return pd.DataFrame(data=res, index=[0])

# Extracts exon start and end positions from a GTF file for specified chromosomes
def get_gtf_splice_pos(gtf_file, chromosomes):
    gtf = pysam.TabixFile(gtf_file)
    poslist = defaultdict(list)
    for chrom in chromosomes:
        for gtf_entry in gtf.fetch(
                chrom,
                parser=pysam.asGTF()):
            if gtf_entry.feature == 'exon':  # Only consider exons
                poslist[chrom].append((gtf_entry.start, gtf_entry.end))
    for key in poslist.keys():
        poslist[key].sort()  # Sort positions for each chromosome
    return poslist
