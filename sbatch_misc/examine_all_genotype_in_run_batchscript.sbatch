#!/bin/bash

#SBATCH --array=1-2%2
#SBATCH --mem=10G
#SBATCH -J examine_all_genotype
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/examine_all_genotype.out

ml rnaseq_pipeline

read align_dir < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/examine_genotype_filelist.txt )

echo "...working on ${align_dir}"

query_path=/scratch/mblab/chasem/rnaseq_pipeline/query/combined_df_20200724.csv


/scratch/mblab/chasem/misc_scripts/examine_all_genotype_coverages_in_run.py ${align_dir} ${query_path}
