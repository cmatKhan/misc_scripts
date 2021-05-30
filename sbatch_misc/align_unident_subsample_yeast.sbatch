#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --array=1-382%20
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/unident_yeast_alignment.out
#SBATCH -J unident_yeast_alignment 

ml novoalign

read fastq_file < <( sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/unident_subsample_lookup.txt )

genome_index=/scratch/mblab/chasem/rnaseq_pipeline/genome_files/S288C_R64/S288C_R64.nix
output_dir_path=/scratch/mblab/mblab.shared/scratch_sequence/yeast_unident_subsample_alignment

fastq_basename=$(basename $fastq_file)
sample_name=${fastq_basename%%.fastq.gz}

# send alignment file to /dev/null. save log only

novoalign -c 8 -o SAM -d ${genome_index} -f ${fastq_file} 1> /dev/null 2> $output_dir_path/${sample_name}_novoalign.log
