#!/bin/bash

#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=10
#SBATCH -J GENOTYPE_deseq_model
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/GENOTYPE_deseq_model.out

ml rnaseq_pipeline

raw_counts=/scratch/mblab/chasem/misc_scripts/deseq_model/results/librarydate_on_libraryprotocol_residuals/raw_residuals.csv
metadata=/scratch/mblab/chasem/misc_scripts/deseq_model/data/filtered_metadata_20200811.csv

output_dir=/scratch/mblab/chasem/misc_scripts/deseq_model/results


/scratch/mblab/chasem/misc_scripts/deseq_model/deseq_model.R -r ${raw_counts} -m ${metadata} -d ~GENOTYPE -o ${output_dir} -n genotype_on_libraryprotocol_librarydate_residuals
