#!/bin/bash

#SBATCH --array=1-58%29
#SBATCH --mem=20G
#SBATCH -J qual_assess
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/qual_assess_revised.out

ml rnaseq_pipeline

read run_dir < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/align_count_lookup.txt )

query_sheet=/scratch/mblab/chasem/rnaseq_pipeline/query/combined_df_20200722.csv

echo "working on ${run_dir}"

quality_assess_1.py -ac ${run_dir} -o ${run_dir} -qs ${query_sheet} --interactive
