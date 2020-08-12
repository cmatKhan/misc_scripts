#!/bin/bash

#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=10
#SBATCH -J LIBRARYDATE_deseq_model
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/LIBRARYDATE_deseq_model.out

ml rnaseq_pipeline

raw_counts=/scratch/mblab/chasem/misc_scripts/deseq_model/results/library_protocol_on_raw_new/raw_residuals.csv
metadata=/scratch/mblab/chasem/misc_scripts/deseq_model/data/filtered_metadata_20200811.csv

output_dir=/scratch/mblab/chasem/misc_scripts/deseq_model/results


/scratch/mblab/chasem/misc_scripts/deseq_model/deseq_model.R -r ${raw_counts} -m ${metadata} -d ~LIBRARYDATE -o ${output_dir} -n librarydate_on_libraryprotocol_residuals_new
