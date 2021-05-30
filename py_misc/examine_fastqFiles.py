#!/usr/bin/env python

import pandas as pd
import os
from rnaseq_tools import utils
from rnaseq_tools.StandardDataObject import StandardData

dataframe_df = pd.read_csv('/scratch/mblab/chasem/rnaseq_pipeline/query/combined_df_20200625.csv')
report_file_path = '/scratch/mblab/chasem/rnaseq_pipeline/reports/fastqFileName_report_20200625.txt'


fastqFiles_dict = {}

sd = StandardData()

def add_to_dict(run_number, library_date, index, fastqFilename):
    """
        add entry to fastqFiles_dict
    """
    try: 
        fastqFiles_dict[run_number][library_date].append((index,fastqFilename)) 
    except KeyError: 
       fastqFiles_dict.setdefault(run_number, {library_date:[(index,fastqFilename)]})

for index,row in dataframe_df.iterrows():
    print("progress: %3.3f" %(int(index)/float(len(dataframe_df))*100))
    try:
        fastqFilename = utils.pathBaseName(str(row['fastqFileName']))
    except AttributeError:
        add_to_dict(str(row['runNumber']), str(row['libraryDate']), index, 'Error_in_fqFileName')
    else:
        fastqFilename = utils.pathBaseName(fastqFilename)+'.fastq.gz'
        run_number = str(row['runNumber'])
        if int(run_number) in sd._run_numbers_with_zeros:
            run_number = sd._run_numbers_with_zeros[int(run_number)]
        library_date = str(row['libraryDate'])
        fastq_path = os.path.join('lts_sequence','run_%s_samples'%run_number,fastqFilename)
        if not isinstance(fastqFilename, str):
            try:
                fastqFiles_dict[run_number][library_date].append((index,fastqFilename))
            except KeyError:
                fastqFiles_dict.setdefault(run_number, {library_date:[(index,fastqFilename)]})
        elif not os.path.isfile(fastq_path):
            add_to_dict(run_number, library_date, index, fastq_path)

with open(report_file_path, 'w') as write_file:
    for key, value in fastqFiles_dict.items():
        write_file.write('\n\n'+key+'\n')
        for sub_key, sub_value in value.items():
            write_file.write('\t'+'\t'+sub_key+'\n\t%s'%sub_value)

print('report in %s'%report_file_path)
