#!/bin/bash

#SBATCH --array=1-87%20
#SBATCH --mem=10G
#SBATCH -J sort_bam
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/sort_bam.out

ml rnaseq_pipeline

read unsorted_bam < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/mblab.shared/scratch_sequence/sorted_convert_to_bam/bam_lookup.txt )

echo "...working on ${unsorted_bam}"

x=$(basename ${unsorted_bam})
sample=${x%%.bam}

samtools sort ${unsorted_bam} -o /scratch/mblab/mblab.shared/scratch_sequence/sorted_convert_to_bam/${sample}_sorted.bam
