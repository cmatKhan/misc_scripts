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
# create output path
output_path = os.path.join(sys.argv[2], utils.pathBaseName(bam_file), '_genotype_check.csv')

qa = QualityAssessmentObject(interactive=True)

organism_object = OrganismData(organism = 'KN99', interactive=True)

genotype_list = ['CKF44_00332', 'CKF44_00896', 'CKF44_01708', 'CKF44_01948', 'CKF44_03366']

# create coverage_df with fastq_simple_name down rows and columns for each genotype
coverage_df = pd.DataFrame({'fastqFilename': utils.pathBaseName(bam_file)})

# credit: https://stackoverflow.com/a/44951376/9708266
coverage_df = coverage_df.reindex(columns=[*coverage_df.columns.to_list(), *genotype_list], fill_value=-1)
coverage_df.reset_index(inplace=True, drop=True)

for index, row in coverage_df.iterrows():
    fastq_simple_name = str(row['fastqFilename'])
    print('...working on %s' %(fastq_simple_name))
    bam_file = [bam_file for bam_file in bam_file_path_list if fastq_simple_name in bam_file][0]
    for list_index in range(len(genotype_list)):
        genotype = genotype_list[list_index]
        print('\t...working on %s which is %s/%s' %(genotype, list_index, len(genotype_list)))
        coverage_df.loc[index, genotype] = qa.calculatePercentFeatureCoverage('CDS', genotype, organism_object.annotation_file, bam_file)

output_path = os.path.join(sys.argv[2], utils.pathBaseName(bam_file), '_genotype_check.csv')
if os.path.isfile(output_path):
    print('OH NO! BAM FILENAME NOT UNIQUE')

coverage_df.to_csv(output_path, index=False)
