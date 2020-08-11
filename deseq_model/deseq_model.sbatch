#!/bin/bash

#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=10
#SBATCH -J deseq_model
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/deseq_model.out

ml rnaseq_pipeline

raw_counts=/scratch/mblab/chasem/misc_scripts/exon_filtered_counts_for_deseq_model.csv
metadata=/scratch/mblab/chasem/misc_scripts/filtered_metadata_for_deseq_model.csv

output_path=/scratch/mblab/chasem/misc_scripts/exon_deseq_model_residuals.csv

/scratch/mblab/chasem/misc_scripts/deseq_model.R -r ${raw_counts} -m ${metadata} -o ${output_path}
