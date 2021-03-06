#!/bin/bash

# subsample bam file in powers of two between 11 and 26 with 5 replicates per
# power for subsample sizes less than the full size of the library. The last size subsampled
# will be the power of two < full size. Full sized library will be saved as
# ${basename_no_ext}_full_${protein_coding_library_size}_protein_coding_subset.bam

# next, a file will be created with a list of the sample names, sbatch scripts will be created
# for both counting and quantifying exon coverage, and sbatch scripts will be submitted
# results will be one directory up in ${name}_count and ${name}_exon_coverage

# INPUT $1: (full, not filtered) bam file from novoalign
# INPUT $2: genome .fasta
# INPUT $3: name of bam (maybe a run_number). Files will be made in $PWD eg) 674_samples, 674_count, 674_exon_coverage, if $3 is 674
# INPUT $4: strandedness of library
# INPUT $5: annotation file

# NOTE: the htseq-count options are -i exon -t gene. You will need to change this manually in the batchfiles if not appropriate

if (( $# < 5 )); then
  echo 'Less than five arguments were passed. Try head -19 /path/to/script for instructions'
  exit 1
fi

# name inputs
input_bam=$1
genome=$2
name=$3
strandedness=$4
annotation_file=$5

# create sample_dir and job_scripts in $PWD
sample_dir=${name}_samples
mkdir -p ${sample_dir}
# make job_scripts directory
mkdir -p job_scripts
# copy this file into job_scripts as well as cmd to execute a given run of this script
rsync -aHv $0 ./job_scripts
echo "$0 $1 $2 $3 $4 $5" > job_scripts/execution_cmd.sh

# create the variables necessary to name the output bamfiles
no_path_input_bam=${input_bam##*/}
basename_no_ext=${no_path_input_bam%_sorted_aligned_reads_with_annote.bam}
filtered_bam_file=${sample_dir}/${basename_no_ext}_protein_coding_counted.bam

echo "filtering bam file for unique, unambiguous, protein coding reads"
# filter bam and deposit in PWD
samtools view -F 4 -q 10 ${input_bam} | grep CKF44 | grep -v alignment_not_unique | grep -v ambiguous | grep -v no_feature | grep -v CNAG | samtools view -T ${genome} -bS > ${filtered_bam_file}

# get full size of filtered protein coding read library
protein_coding_library_size=$(samtools view ${filtered_bam_file} | wc -l)
echo "size of filtered bam is: ${protein_coding_library_size}"


# re-name the full protein coding library
full_subset_moniker=full_${protein_coding_library_size}
# remember to include the ${sample_dir} in the path
protein_coding_library=${sample_dir}/${basename_no_ext}_${full_subset_moniker}_protein_coding_subset.bam
mv ${filtered_bam_file} ${protein_coding_library}
echo "filtered bam in ${protein_coding_library}"

create samples. largest protein coding count in crypto is 5.6*10^7
for power in {11..26};do
  for seed in {1..5};do

    # set power given power of 2 to variable n
    power_of_two=$(echo "2^${power}"|bc)
    echo "power of two is ${power_of_two}"

   # get the fraction of the input bam file that the given power_of_two corresponds to
    frac=$( echo "scale=10 ; ${power_of_two} / ${protein_coding_library_size}" | bc )
    echo "fraction of protein coding library size is ${frac}"

    # drop leading zero from $frac and replace with $seed
    if (( $(echo "$frac < 1" | bc -l) )); then
      nozero=$(echo $frac | sed 's/^0*//')
      seed_frac=${seed}${nozero}
      # use sametools with the given seed (see outer loop) and the fraction of the library to subset thereads
      echo "subsetting out ${power_of_two} reads, replicate ${seed}, with seed_frac ${seed_frac}"
      samtools view -bs ${seed_frac} ${protein_coding_library} > ${sample_dir}/${basename_no_ext}_${seed}_${power}_protein_coding_subset.bam
  fi
  done
done

# create list of files in job_scripts directory
subset_bam_file_list=job_scripts/${name}_subset_file_list.txt
find $PWD/${sample_dir} -name "*.bam" -type f > ${subset_bam_file_list}
num_samples=$( cat ${subset_bam_file_list} | wc -l)

# create ${name}_count directory
mkdir -p ${name}_count
# create ${name}_exon coverage directory
mkdir -p ${name}_exon_coverage directory
# create slurm log directory
mkdir -p sbatch_log

# create python script for sbatch exon coverage script
printf "#!/usr/bin/env python\n" > job_scripts/exonic_coverage.py

printf "import sys\n" >> job_scripts/exonic_coverage.py
printf "from rnaseq_tools.CryptoQualityAssessmentObject import CryptoQualityAssessmentObject\n" >> job_scripts/exonic_coverage.py

printf "qa = CryptoQualityAssessmentObject(interactive=True)\n" >> job_scripts/exonic_coverage.py

printf "qa.calculateExonicCoverage(sys.argv[1], sys.argv[2])" >> job_scripts/exonic_coverage.py

# make the python script executable
chmod +x job_scripts/exonic_coverage.py

# create sbatch script for exon coverage
printf "#!/bin/bash\n\n" > job_scripts/exon_coverage_${name}.sbatch

printf "#SBATCH --mem=20G\n#SBATCH --array=1-${num_samples}%%20\n#SBATCH -D ./\n#SBATCH -o sbatch_log/exon_coverage_slurm.out\n\n" >> job_scripts/exon_coverage_${name}.sbatch

printf "ml rnaseq_pipeline\n\n" >> job_scripts/exon_coverage_${name}.sbatch

printf "read bam_file < <( sed -n \${SLURM_ARRAY_TASK_ID}p ${subset_bam_file_list} )\n\n" >> job_scripts/exon_coverage_${name}.sbatch

printf "job_scripts/exonic_coverage.py \${bam_file} ${name}_exon_coverage" >> job_scripts/exon_coverage_${name}.sbatch

# create sbatch script for counts
printf "#!/bin/bash\n\n" > job_scripts/count_${name}.sbatch

printf "#SBATCH --mem=20G\n#SBATCH --array=1-${num_samples}%%20\n#SBATCH -D ./\n#SBATCH -o sbatch_log/count_slurm.out\n\n" >> job_scripts/count_${name}.sbatch

printf "ml htseq\n\n" >> job_scripts/count_${name}.sbatch

printf "read bam_file < <( sed -n \${SLURM_ARRAY_TASK_ID}p ${subset_bam_file_list} )\n" >> job_scripts/count_${name}.sbatch

printf "no_path_input_bam=\${bam_file##*/}\n" >> job_scripts/count_${name}.sbatch

printf "sample_name=\${no_path_input_bam%%_protein_coding_subset.bam}\n\n" >> job_scripts/count_${name}.sbatch

printf "htseq-count -f bam -s ${strandedness} -t exon -i gene \${bam_file} ${annotation_file} 1> ${name}_count/\${sample_name}_read_count.tsv 2> ${name}_count/\${sample_name}_htseq.log" >> job_scripts/count_${name}.sbatch

# submit sbatch scripts
echo "submitting job_scripts/count_${name}.sbatch"
sbatch job_scripts/count_${name}.sbatch
echo "submitting job_scripts/exon_coverage_${name}.sbatch"
sbatch job_scripts/exon_coverage_${name}.sbatch
