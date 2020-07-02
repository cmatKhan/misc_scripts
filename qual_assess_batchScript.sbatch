#!/bin/bash

#SBATCH --array=1-31%31
#SBATCH --mem=20G
#SBATCH -J qual_assess

ml rnaseq_pipeline

read run_dir < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/qual_assess_lookup.txt )

quality_assess_1.py -r ${run_dir} -o ${run_dir} -cc -qs /scratch/mblab/chasem/rnaseq_pipeline/query/nextflow_full_test.csv --interactive
