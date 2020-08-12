#!/bin/bash

#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=10
#SBATCH -J deseq_model
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/GENOTYPE_deseq_model.out

ml rnaseq_pipeline

raw_counts=/scratch/mblab/chasem/misc_scripts/deseq_model/data/exon_filtered_counts_20200811.csv
metadata=/scratch/mblab/chasem/misc_scripts/deseq_model/data/filtered_metadata_20200811.csv

output_dir=/scratch/mblab/chasem/misc_scripts/deseq_model/results


/scratch/mblab/chasem/misc_scripts/deseq_model/deseq_model.R -r ${raw_counts} -m ${metadata} -d ~LIBRARYPROTOCOL -o${output_dir} -n library_protocol_on_raw
