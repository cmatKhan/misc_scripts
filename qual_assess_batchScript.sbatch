#!/bin/bash

#SBATCH --array=1-11%11
#SBATCH --mem=20G
#SBATCH -J qual_assess
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/re_assess_qual_assess_1.out

ml rnaseq_pipeline

read run_dir < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/re_assess_qual_assess_lookup.txt )

quality_assess_1.py -r ${run_dir} -o ${run_dir} -cc -qs /scratch/mblab/chasem/rnaseq_pipeline/query/nextflow_full_test.csv --interactive
