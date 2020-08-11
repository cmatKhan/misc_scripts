#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --array=1-131%20
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/unident_crypto.out
#SBATCH -J unident_crypto

ml novoalign

read fastq_file run_num< <( sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/unident_fastq_files_lookup.txt )

crypto_index=/scratch/mblab/chasem/rnaseq_pipeline/genome_files/KN99/KN99_genome_fungidb.nix
output_dir_path_crypto=/scratch/mblab/mblab.shared/scratch_sequence/unident_alignment/crypto/bams
fastq_basename=$(basename $fastq_file)
sample_name=${fastq_basename%%.fastq.gz}

fastq_path=/scratch/mblab/mblab.shared/scratch_sequence/unident_crypto_fastq/${fastq_file}
 
echo "aligning ${fastq_path}"
# send alignment file to /dev/null. save log only crypto
novoalign -c 8 -o SAM -d ${crypto_index} -f ${fastq_path} 1> ${output_dir_path_crypto}/${sample_name}_aligned_reads.bam 2> /dev/null
