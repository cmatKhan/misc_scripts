#!/bin/bash

#SBATCH --array=1-11%11
#SBATCH --mem=20G
#SBATCH -J samtools_index
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/samtools_re_index_07_01.out

ml rnaseq_pipeline

read align_dir < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/re_index_lookup.txt )

for bam_file in ${align_dir}/*; do
    samtools index ${bam_file}
done
