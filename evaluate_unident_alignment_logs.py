#!/usr/bin/env python

import glob
import pandas as pd
import os
import re
import sys

crypto_log_dirpath = sys.argv[1]

yeast_log_dirpath = sys.argv[2]

crypto_log_list = glob.glob(crypto_log_dirpath+'/*.log')
yeast_log_list = glob.glob(yeast_log_dirpath+'/*.log')

sample_list = [os.path.basename(x).replace('_novoalign.log','') for x in crypto_log_list]

unident_eval_df = pd.DataFrame({'sample': sample_list, 'crypto_uniq_align': None, 'yeast_uniq_align': None, 'guess': None})

unique_alignment_regex = r"(?<=Unique Alignment:\s)\s*\d*"

for index, row in unident_eval_df.iterrows():
    novoalign_name = row['sample']
    print('...checking %s'%novoalign_name)

    crypto_alignment_log_path = os.path.join(crypto_log_dirpath, novoalign_name+'_novoalign.log')
    crypto_alignment_log_file = open(crypto_alignment_log_path, 'r')
    crypto_alignment_log_text = crypto_alignment_log_file.read()

    yeast_alignment_log_path = os.path.join(yeast_log_dirpath, novoalign_name+'_novoalign.log')
    yeast_alignment_log_file = open(yeast_alignment_log_path, 'r')
    yeast_alignment_log_text = yeast_alignment_log_file.read()

    crypto_unique_alignment = int(re.findall(unique_alignment_regex, crypto_alignment_log_text)[0]) / float(50000.0)
    yeast_unique_alignment = int(re.findall(unique_alignment_regex, yeast_alignment_log_text)[0]) / float(50000.0)

    if crypto_unique_alignment > .6 and yeast_unique_alignment < .1:
        guess = 'crypto'
    elif crypto_unique_alignment < .1 and yeast_unique_alignment > .6:
        guess = 'yeast'
    else:
        guess = 'unknown'

    unident_eval_df.loc[index, 'crypto_uniq_align'] = crypto_unique_alignment
    unident_eval_df.loc[index, 'yeast_uniq_align'] = yeast_unique_alignment
    unident_eval_df.loc[index, 'guess'] = guess
    
    crypto_alignment_log_file.close()
    yeast_alignment_log_file.close()

output_path = os.path.join(sys.argv[3], 'unidentified_fastq_evaluation.csv')
print('...writing to %s' %output_path)
unident_eval_df.to_csv(output_path)
