#!/bin/bash

#SBATCH --array=1-56%20
#SBATCH --mem=10G
#SBATCH -J qual_assess_audit
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/qual_assess_audit.out

ml rnaseq_pipeline

read qual_assess_path < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/qual_assess_audit_lookup.txt )

query_path=/scratch/mblab/chasem/rnaseq_pipeline/query/nextflow_full_test.csv

/scratch/mblab/chasem/misc_scripts/update_qa1_audit.py ${query_path} ${qual_assess_path}
