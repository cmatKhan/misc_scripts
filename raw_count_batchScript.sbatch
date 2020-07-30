#!/bin/bash

#SBATCH --array=1-58%29
#SBATCH --mem=20G
#SBATCH -J raw_count
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/raw_count.out

ml rnaseq_pipeline

read count_dir < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/count_directory_lookup.txt )

raw_counts.py -e ${count_dir} -l /scratch/mblab/chasem/rnaseq_pipeline/genome_files/KN99/KN99_gene_id_list.txt 
