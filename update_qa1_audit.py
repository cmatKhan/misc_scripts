#!/usr/bin/env python

# INPUT 1: query_df
# INPUT 2: the sheet to audit

import sys
import pandas as pd
from rnaseq_tools.CryptoQualityAssessmentObject import CryptoQualityAssessmentObject

print('...creating quality assessment object')
qa = CryptoQualityAssessmentObject(interactive=True)

qa.query_df = pd.read_csv(sys.argv[1])

print('...auditing sheet')
df = qa.auditQualAssessDataFrame(sys.argv[2])

print('...writing to %s' %sys.argv[2])
df.to_csv(sys.argv[2], index=False)
