#!/usr/bin/env python

"""
first argument is a path to a bam file.
second argument is path to directory for output. name of file is handled in script
Note: the genotype list is hard coded in

    outputs a dataframe with fastqFilenames from run in query_df and columns as all of the unique perturbed genotypes (minus wildtype)
           outputs into rundirectory

"""


from rnaseq_tools.QualityAssessmentObject import QualityAssessmentObject
from rnaseq_tools.OrganismDataObject import OrganismData
from rnaseq_tools import utils
import pandas as pd
import re
import os
import sys

# extract bam file
bam_file = sys.argv[1]
simple_name = utils.pathBaseName(bam_file)
# create output path
output_path = os.path.join(sys.argv[2], simple_name + '_genotype_check.csv')

qa = QualityAssessmentObject(interactive=True)

organism_object = OrganismData(organism = 'KN99', interactive=True)

genotype_list = ['CKF44_00332', 'CKF44_00896', 'CKF44_01708', 'CKF44_01948', 'CKF44_03366']

# create coverage_df with fastq_simple_name down rows and columns for each genotype
simple_name_dict = {simple_name: [-1,-1,-1,-1,-1]}
coverage_df = pd.DataFrame.from_dict(simple_name_dict, orient='index', columns=genotype_list)
coverage_df['fastqFilename'] = coverage_df.index
coverage_df.reset_index(drop=True, inplace=True)
print(coverage_df)

for index, row in coverage_df.iterrows():
    fastq_simple_name = str(row['fastqFilename'])
    print('...working on %s' %(fastq_simple_name))
    for list_index in range(len(genotype_list)):
        genotype = genotype_list[list_index]
        print('\t...working on %s which is %s/%s' %(genotype, list_index, len(genotype_list)))
        coverage_df.loc[index, genotype] = qa.calculatePercentFeatureCoverage('CDS', genotype, organism_object.annotation_file, bam_file)

if os.path.isfile(output_path):
    print('OH NO! BAM FILENAME NOT UNIQUE')

coverage_df.to_csv(output_path, index=False)
