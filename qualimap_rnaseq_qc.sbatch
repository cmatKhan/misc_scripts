#!/bin/bash

#SBATCH --array=1-1039%20
#SBATCH --mem=30G
#SBATCH -J qualimap_rnaseq_qc
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/qualimap_rnaseq_qc.out

read bam_file annote_file strandedness < <( sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/qualimap_rnaseq_qc_lookup.txt )
run_dir=$(dirname $(dirname $bam_file))
genome=/scratch/mblab/chasem/rnaseq_pipeline/genome_files/KN99/KN99_genome_fungidb.fasta

if [[ ! -f "$bam_file" ]]; then
    bam_filepath_with_leading_zero=${bam_file/run_/run_0}
    bam_file=${bam_filepath_with_leading_zero}
fi

if [[ ! -f ${bam_file} ]]; then
    echo "${bam_file} DNE"

else

    ml rnaseq_pipeline

    filtered_bam=/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/bam_${SLURM_ARRAY_TASK_ID}.bam
    echo "filtering bam file ${bam_file}"
    samtools view ${bam_file} | grep -v "CNAG_[[:digit:]]*"| samtools view -bS -T ${genome} > ${filtered_bam}
    echo "qualimapping ${bam_file}"
    qualimap rnaseq --java-mem-size=16G -bam ${filtered_bam} -gtf ${annote_file} -p ${strandedness} -outdir ${run_dir}
    echo "done with ${filtered_bam}. Removing it. check ${run_dir} for output"
    rm ${filtered_bam}
fi
