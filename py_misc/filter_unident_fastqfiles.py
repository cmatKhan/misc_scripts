#!/usr/bin/env python

from multiprocessing import Process
import glob
import subprocess
import pandas as pd
import os
import sys

unident_fastq_sheet_dir_files = glob.glob(os.path.join(sys.argv[1],'*.csv'))
path_to_lts_sequence = sys.argv[2]


def count_fastq_reads(fastq_filelist_path, path_to_lts_sequence):

    fastq_filelist_df = pd.read_csv(fastq_filelist_path)


    for index, row in fastq_filelist_df.iterrows():
        run_number = row['RUNNUMBER']
        run_directory = run_number + '_samples'
        fastq_basename = row['fastq_basename']
        fastq_path = os.path.join(path_to_lts_sequence, run_directory, fastq_basename)

        if not os.path.isfile(fastq_path):
            print('Cannot find: %s' %fastq_path)

        num_reads_cmd = 'zcat %s | wc -l' %fastq_path
        num_reads = subprocess.getoutput(num_reads_cmd)

        if int(num_reads) > 1*10^6:
            with open('%s_samples_million_reads_or_more.txt'%run_number, 'a') as file:
                file.write('%s\n'%fastq_path)

for run_dir_sheet in unident_fastq_sheet_dir_files:
    print('starting to process %s'%run_dir_sheet)
    p = Process(target=count_fastq_reads, args=(run_dir_sheet, path_to_lts_sequence))
    p.start()
