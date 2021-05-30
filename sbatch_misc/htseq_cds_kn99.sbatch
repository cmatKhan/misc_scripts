#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH --array=1-1039%20
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/htseq_cds_kn99.out
#SBATCH -J htseq_cds_kn99 

ml htseq 

read simple_name bam_filepath strandedness annotation_file < <( sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/kn99_induction_cds_counts_lookup.txt )

if [[ ! -f "$bam_filepath" ]]; then
    bam_filepath_with_leading_zero=${bam_filepath/run_/run_0}
    bam_filepath=${bam_filepath_with_leading_zero}
fi

if [[ ! -f "$bam_filepath" ]]; then
    echo "${bam_filepath} DNE!"
fi
output_fullpath_count_file=/scratch/mblab/chasem/rnaseq_pipeline/experiments/cds_counts/${simple_name}_read_counts.tsv
output_fullpath_log=/scratch/mblab/chasem/rnaseq_pipeline/experiments/cds_counts/${simple_name}_htseq_log.txt

echo "count reads on CDS in ${simple_name}"

htseq-count -f bam -s ${strandedness} -t CDS -i gene_id ${bam_filepath} ${annotation_file} 1> ${output_fullpath_count_file} 2> ${output_fullpath_log} 
