#!/usr/bin/env python

import glob
import pandas as pd
import os
import sys

path_to_metadata = sys.argv[1]

bam_search_pattern = sys.argv[2]

bam_file_list = glob.glob(bam_search_pattern, recursive=True)

df = pd.read_csv(path_to_metadata)

bam_dict = {}
group_dict = {}

group_num = 1
for bam_file in bam_file_list:
    sample_name = os.path.basename(bam_file).replace('_sorted_aligned_reads_with_annote.bam', '')
    try:
        genotype =df[df.FASTQFILENAME.str.startswith(sample_name)]['GENOTYPE'].values[0]
    except IndexError:
        continue
    genotype =df[df.FASTQFILENAME.str.startswith(sample_name)]['GENOTYPE'].values[0]
    strain = df[df.FASTQFILENAME.str.startswith(sample_name)]['STRAIN'].values[0]
    print(strain)
    if pd.isna(strain):
        strain = 'unknown_strain'
    geno_strain = genotype + '_' + strain

    bam_dict[bam_file] = geno_strain

    try:
        group_dict[geno_strain]
    except KeyError:
        group_dict[geno_strain] = 'group_%s' %group_num
        group_num = group_num + 1

with open('qualimap_multi.tsv', 'w') as quali_data:
    for bam_file, geno_strain in bam_dict.items():
        quali_data.write('%s\t%s\t%s\n' %(geno_strain, bam_file, group_dict[geno_strain]))
