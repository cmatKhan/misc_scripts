#!/usr/bin/env python

import sys
from rnaseq_tools.QualityAssessmentObject import QualityAssessmentObject

print('...creating quality assessment object')
qa = QualityAssessmentObject(interactive=True)

print('...auditing sheet')
df = qa.auditQualAssess1(sys.argv[1])

print('...writing to %s' %sys.argv[1])
df.to_csv(sys.argv[1], index=False)
