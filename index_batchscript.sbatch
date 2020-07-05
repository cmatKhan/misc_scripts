#!/bin/bash

#SBATCH --array=1-${1}%${1}
#SBATCH --mem=20G
#SBATCH -J samtools_index
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/samtools_re_index_07_01.out

## index the files submitted in the lookup text (see create_index_lookup.sh)
## INPUT $1: the line number of lines in lookup unindexed_align_file
## INPUT $2: the lookup file of alignment files to index

ml rnaseq_pipeline

read unindexed_align_file < <(sed -n ${SLURM_ARRAY_TASK_ID}p ${2} )

samtools index ${unindexed_align_file}
