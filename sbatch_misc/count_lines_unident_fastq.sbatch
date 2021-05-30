#!/bin/bash

#SBATCH --array=1-804%20
#SBATCH --mem=20G
#SBATCH -J unident_eval
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/unident_eval.out

ml rnaseq_pipeline

read fastq_file runnumber < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/unident_fastq_files_lookup.txt )

fastq_file_path=/lts/mblab/Crypto/rnaseq_data/lts_sequence/${runnumber}/${fastq_file}

x=$(zcat ${fastq_file_path} | head -300000 | wc -l)

if [[ $x -eq 30000 ]]; then
    echo "yessiree $x"
fi
