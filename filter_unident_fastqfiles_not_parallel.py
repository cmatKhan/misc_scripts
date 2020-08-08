#!/usr/bin/env python

import subprocess
import pandas as pd
import os
import sys


fastq_filelist_path = sys.argv[1]
path_to_lts_sequence = sys.argv[2]

def count_fastq_reads(fastq_filelist_path, path_to_lts_sequence):

    fastq_filelist_df = pd.read_csv(fastq_filelist_path)


    for index, row in fastq_filelist_df.iterrows():
        run_number = row['RUNNUMBER']
        run_directory = run_number + '_samples'
        fastq_basename = row['fastq_basename']
        fastq_path = os.path.join(path_to_lts_sequence, run_directory, 'Brent_large', fastq_basename)

        if not os.path.isfile(fastq_path):
            print('Cannot find: %s' %fastq_path)

        num_reads_cmd = 'zcat %s | wc -l' %fastq_path
        num_reads = subprocess.getoutput(num_reads_cmd)

        if int(num_reads) > 1*10^6:
            with open('%s_samples_million_reads_or_more.txt'%run_number, 'a') as file:
                file.write('%s\n'%fastq_path)

count_fastq_reads(fastq_filelist_path, path_to_lts_sequence)
