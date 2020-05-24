#!/usr/bin/env python3

"""
    remove the sequences identified by their header in the filter_file.txt
 
        eg) if there are chromosomes 1 to 3 in fasta in format >chr1
                                                               ATAGCA...
                                                               >chr2
                                                               ATAGCCC...
                                                               >chr2
                                                               AAAATCC...

            and the filter_file.txt has the following: chr1
                                                       chr3

            the output will be >chr2
                               ATAGCCC...

    usage: remove_scaffold_from_genome_fasta.py -g genome.fasta -f filter_file.txt

    written_by: Kamil S Jaron (https://bioinformatics.stackexchange.com/a/3940)
    modified_by: chase.mateusiak@gmail.com
    dependencies: Biopython, argparse
"""
from Bio import SeqIO
import argparse
import os
import sys

def main(argv):
    """
        main method
    """
    print('...parsing input')
    args = parse_args(argv)
    try:
        if not os.path.isfile(args.genome):
            raise FileNotFoundError('GenomeFileNotFound')
    except FileNotFoundError:
        print('The path to the genome file is not valid. Check and try again')
    else:
        genome_fasta_path = args.genome
    try:
        if not os.path.isfile(args.filter_file):
            raise FileNotFoundError('FilterFileNotFound')
    except FileNotFoundError:
        print('The path to the filter file is not valid. Check and try again.')
    else:
        filter_file_path = args.filter_file
    finally:
        if not (os.path.isfile(args.genome) and os.path.isfile(args.filter_file)):
            sys.exit()

    # read filter_file_path and create a set
    print('...reading filter_file')
    with open(filter_file_path, 'r') as file:
        scaffold_to_remove_set = set(line.strip() for line in file)

    # create generator from fasta file
    genome_generator = SeqIO.parse(genome_fasta_path, 'fasta')
    # loop over genome, if the header is in the filter file
    print('...filtering input %s for scaffolds in %s' %(genome_fasta_path, filter_file_path))
    for seq_record in genome_generator:
        try:
            # seq_record.name is in scaffold_to_remove_set, pop it out of the set and continue in loop
            scaffold_to_remove_set.remove(seq_record.name)
        except KeyError:
            print(seq_record.format('fasta').strip())
    if len(scaffold_to_remove_set) != 0:
        print('%s not found in %s' %(scaffold_to_remove_set, genome_fasta_path))


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="filter scaffold from genome fasta")
    parser.add_argument("-g", "--genome", required=True,
                        help="[REQUIRED] genome in fasta format")
    parser.add_argument("-f", "--filter_file", required=True,
                        help="[REQUIRED] file with the headers to remove (see doc string at top of file, use head, to see example/description)")

    args = parser.parse_args(argv[1:])
    return args


if __name__ == "__main__":
    main(sys.argv)
