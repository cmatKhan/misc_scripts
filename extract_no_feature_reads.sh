
# first, remove tab from front of .sam output by htseq
sed '/\t//' <the htseq sam> > outfile

# next, paste the column onto the end of each line in the bam grep for no feature and write out

samtools view <bam> | paste - b > outbam | grep no_feature > outfile
