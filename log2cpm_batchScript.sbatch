#!/bin/bash

#SBATCH --array=1-57%29
#SBATCH --mem=20G
#SBATCH -J log2_cpm_count
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/log2cpm_count.out

ml rnaseq_pipeline

read raw_count_file < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/count_lookup.txt )

count_file_dir=$(dirname ${raw_count_file})

output_name=${count_file_dir}/log2_cpm.csv

log2_cpm.R -r ${raw_count_file} -g KN99 -o ${output_name}
