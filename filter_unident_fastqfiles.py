#!/usr/bin/env python

import subprocess
import pandas as pd
import os
import sys

fastq_filelist_path = sys.argv[1]

fastq_filelist_df = pd.read_csv(fastq_filelist_path)

path_to_lts_sequence = sys.argv[2]

for index, row in fastq_filelist_df.iterrows():
    run_number = row['RUNNUMBER']
    run_directory = run_number + '_samples'
    fastq_basename = row['fastq_basename']
    fastq_path = os.path.join(path_to_lts_sequence, run_directory, fastq_basename)

    if not os.path.isfile(fastq_path):
        print('Cannot find: %s' %fastq_path)


    print('...counting lines in %s' %fastq_path)
    num_reads_cmd = 'zcat %s | wc -l' %fastq_path
    num_reads = subprocess.getoutput(num_reads_cmd)
    print('Done. There are %s lines' %num_reads)

    if int(num_reads) > 1*10^6:
        print('Adding %s to sample_million_reads_or_more.txt' %fastq_path)
        with open('samples_million_reads_or_more.txt', 'a') as file:
            file.write('%s\n'%fastq_path)
