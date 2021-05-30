#!/bin/bash

#SBATCH --array=1-470%47
#SBATCH --mem=20G
#SBATCH -J samtools_index
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/samtools_re_index_07_02.out

## index the files submitted in the lookup text (see create_index_lookup.sh)
## YOU MUST CHANGE THE ARRAY NUMBER TO MATCH THE LINE NUMBER OF LOOKUP FILE AND THE DATE FOR -O
## INPUT $2: the lookup file of alignment files to index

ml rnaseq_pipeline

read unindexed_align_file < <(sed -n ${SLURM_ARRAY_TASK_ID}p ${1} )

samtools index ${unindexed_align_file}
