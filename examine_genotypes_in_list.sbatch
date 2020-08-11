#!/bin/bash

#SBATCH --array=1-87%20
#SBATCH --mem=10G
#SBATCH -J examine_bam_for_genotype
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/examine_bam_for_genotype.out

ml rnaseq_pipeline

read sorted_bam < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/unident_crypto_sorted_lookup.txt )

echo "...working on ${sorted_bam}"

/scratch/mblab/chasem/misc_scripts/calculate_coverages_for_bam.py ${sorted_bam} /scratch/mblab/chasem/misc_scripts/unident_output
