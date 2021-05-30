#!/usr/bin/env bash

#  input:
#      -h --help should display this docstring
#      -g --gff an annotation file
#      -f --fasta genome .fa
#      -o --output output filename



last_line_docstring=16


main(){
  # main method, called at bottom of script after all functions read in

  # parse cmd line input
  parseArgs "$@"
  # verify cmd line input
  checkInput

  agat_sp_statistics.pl -g ${gff} -f ${fasta} -o ${output}
}

checkInput(){
  # check input, raise errors
  # TODO: should this go to 2 or 1?
  if [[ ! -e $gff ]]; then
      echo "RunFastQCInputError: fastq ${gff} file does not exist"
      exit 1
  fi
}

parseArgs(){
  #    parse cmd line input, set global variables
  #    usage: Assuming this is called from main, and the entire cmd line argument array was passed to main, parseArgs "$@"
  #    input: cmd line input passed in main method via $@

  while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
    -h | --help )
      head -${last_line_docstring} $0
      exit
      ;;
    -g | --gff )
      shift; gff=$1
      ;;
    -g | --gff )
      shift; gff=$1
      ;;
    -f | --fasta )
      shift; fasta=$1
      ;;
    -o | --output )
      shift; gff=$1
      ;;
  esac; shift; done
  if [[ "$1" == '--' ]]; then shift; fi

}

main "$@"

