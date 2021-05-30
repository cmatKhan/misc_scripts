#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --array=1-804%20
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/unident_alignment.out
#SBATCH -J unident_crypto_alignment 

ml novoalign

read fastq_file run_num< <( sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/unident_fastq_files_lookup.txt )

crypto_index=/scratch/mblab/chasem/rnaseq_pipeline/genome_files/KN99/KN99_genome_fungidb.nix
yeast_index=/scratch/mblab/chasem/rnaseq_pipeline/genome_files/S288C_R64/S288C_R64.nix
output_dir_path_crypto=/scratch/mblab/mblab.shared/scratch_sequence/unident_alignment/crypto
output_dir_path_yeast=/scratch/mblab/mblab.shared/scratch_sequence/unident_alignment/yeast
fastq_basename=$(basename $fastq_file)
sample_name=${fastq_basename%%.fastq.gz}

fastq_path=/scratch/mblab/mblab.shared/scratch_sequence/unident_fastq_all/${fastq_file}

x=$(zcat ${fastq_path} | head -1000000 | wc -l)

if [[ x -eq 1000000 ]]; then
  
  echo "aligning ${fastq_file}"
  # send alignment file to /dev/null. save log only crypto
  novoalign -c 8 -o SAM -d ${crypto_index} -f ${fastq_path} 1> /dev/null 2> ${output_dir_path_crypto}/${sample_name}_novoalign.log
  # send alignment file to /dev/null. save log only yeast
  novoalign -c 8 -o SAM -d ${yeast_index} -f ${fastq_path} 1> /dev/null 2> ${output_dir_path_yeast}/${sample_name}_novoalign.log

else

 echo ${fastq_file} has only $x reads

fi
