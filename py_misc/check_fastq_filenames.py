#!/usr/bin/env python

import pandas as pd
import os
import sys

df = pd.read_csv(sys.argv[1])

for index, row in df.iterrows():
    fastq = row$fastqFileName
    run_num = str(int(float(row$runNumber)))
    path = '/lts/mblab/Crypto/rnaseq_data/lts_sequence/run_%s_samples/%s'%(run_num, fastq)
    if not os.path.exists(fastq):
        print(path)
