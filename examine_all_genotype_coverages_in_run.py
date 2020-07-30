#!/usr/bin/env python

"""
    first argument, path to run align subdir (otuput of nextflow align count)
    second argument, path to query df

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

qa = QualityAssessmentObject(interactive=True)

organism_object = OrganismData(organism = 'KN99', interactive=True)

run_dir_align_filepath = sys.argv[1]

run_dir = utils.dirPath(run_dir_align_filepath)
run_number = re.search('\d+', run_dir).group(0)

query_sheet_path = sys.argv[2]

query_df = utils.readInDataframe(query_sheet_path)

# get query_df for just given runNumber
run_query_df = query_df[query_df['runNumber']==run_number]

# get unique list of genotype column in given run
genotype_list = list(run_query_df['genotype'].unique())
# drop wildtype
genotype_list.remove('CNAG_00000')
# split double KO, creates list of lists
genotype_list = [x.split('.') for x in genotype_list]
# flatten list credit: https://stackoverflow.com/a/952952/9708266
genotype_list = [x for sublist in genotype_list for x in sublist]

# extract bam file
bam_file_path_list = utils.extractFiles(run_dir_align_filepath, '.bam')

# create coverage_df with fastq_simple_name down rows and columns for each genotype
coverage_df = pd.DataFrame({'fastqFilename': run_query_df['fastqFileName']})
# covert fastq files to simple names
coverage_df.fastqFilename = coverage_df.fastqFilename.apply(lambda x: utils.pathBaseName(x))
# credit: https://stackoverflow.com/a/44951376/9708266
coverage_df = coverage_df.reindex(columns=[*coverage_df.columns.to_list(), *genotype_list], fill_value=-1)
coverage_df.reset_index(inplace=True, drop=True)

for index, row in coverage_df.iterrows():
    fastq_simple_name = str(row['fastqFilename'])
    print('...working on %s/%s: %s' %(index, len(coverage_df), fastq_simple_name))
    bam_file = [bam_file for bam_file in bam_file_path_list if fastq_simple_name in bam_file][0]
    for list_index in range(len(genotype_list)):
        genotype = genotype_list[list_index].replace('CNAG', 'CKF44')
        print('\t...working on %s which is %s/%s' %(genotype, list_index, len(genotype_list)))
        coverage_df.loc[index, genotype] = qa.calculatePercentFeatureCoverage('CDS', genotype, organism_object.annotation_file, bam_file)

all_run_genotype_path = os.path.join(run_dir, 'all_run_genotype_check.csv')

coverage_df.to_csv(all_run_genotype_path, index=False)
