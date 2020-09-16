#!/bin/bash

#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=10
#SBATCH -J parallel_tester
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/parallel_tester.out

ml rnaseq_pipeline

/scratch/mblab/chasem/misc_scripts/deseq_model/sbatch_scripts/parallel_tester.R
