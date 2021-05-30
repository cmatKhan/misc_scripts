#!/usr/bin/env python

# INPUT 1: query_df path
# INPUT 2: path to the sheet to audit
# INPUT 3: path to align count dir (directory with alignment files, count files, log2 cpm, etc)

import sys
import pandas as pd
from rnaseq_tools import utils
from rnaseq_tools.CryptoQualityAssessmentObject import CryptoQualityAssessmentObject

print('...getting list of bam files')
bam_file_list = utils.extractFiles(sys.argv[3], '.bam')

print('...creating quality assessment object')
qa = CryptoQualityAssessmentObject(interactive=True)

query_path = sys.argv[1]

qual_assess_path = sys.argv[2]

print('...auditing sheet')
df = qa.auditQualAssessDataFrame(query_path, qual_assess_path, bam_file_list)

print('...writing to %s' %sys.argv[2])
df.to_csv(sys.argv[2], index=False)
