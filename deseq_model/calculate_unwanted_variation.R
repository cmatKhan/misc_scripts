#!/usr/bin/env Rscript

# calculate k unwanted variation columns using RUVseq, output as columnar data
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(RUVSeq))

main = function(parsed_cmd_line_args){
  
  # parse cmd line arguments
  print('parsing cmd line input')
  count_log_norm_path = parsed_cmd_line_args$count_log_norm_path
  residual_log_norm_path = parsed_cmd_line_args$residual_log_norm_path
  protein_coding_gene_path = parsed_cmd_line_args$protein_coding_gene_path
  num_unwanted_variation_columns = as.double(parsed_cmd_line_args$num_unwanted_variation_columns)
  output_dirpath = parsed_cmd_line_args$output_directory
  output_name = parsed_cmd_line_args$output_name
  
  # read in data
  print('reading in data')
  norm_count_with_pseudocount_matrix = as.matrix(read_csv(count_log_norm_path))
  residual_log_norm_matrix = as.matrix(read_csv(residual_log_norm_path))
  protein_coding_gene_vector = read_csv(protein_coding_gene_path, col_names = FALSE)$X1
  # attach gene column to norm and residual matrix as rownames
  rownames(norm_count_with_pseudocount_matrix) = protein_coding_gene_vector
  rownames(residual_log_norm_matrix) = protein_coding_gene_vector
  
  # calculate unwanted variation columns
  print('calculating unwanted variation')
  ruvr_object = RUVr(norm_count_with_pseudocount_matrix, protein_coding_gene_vector, num_unwanted_variation_columns, residual_log_norm_matrix, isLog = TRUE)

  # write out
  output_name_path = paste(output_dirpath, paste0(output_name, '.csv'), sep='/')
  print(paste0('writing unwanted variation sheet to: ', output_name_path))
  write_csv(as_tibble(ruvr_object$W), output_name_path)
  
} # end main()

parseArguments <- function() {
  # parse and return cmd line input
  
  option_list <- list(
    make_option(c('-c', '--count_log_norm_path'),
                help='normalized, log2 counts (genes x samples)'),
    make_option(c('-r', '--residual_log_norm_path'), 
                help='residuals in normalized, logged scale'),
    make_option(c('-g', '--protein_coding_gene_path'), 
                help='path to protein coding only genes (plus markers)'),
    make_option(c('-k', '--num_unwanted_variation_columns'), 
                help='number of unwanted variation covariates to calculate (k in the RUVseq vignette/argument description)'),
    make_option(c('-o', '--output_directory'), 
                help='path to directory to output results'),
    make_option(c('-n', '--output_name'),
                help='name of the ruvr output sheet, eg genotype_model_10_unwanted_cofactors (no extension)'))
  
  args <- parse_args(OptionParser(option_list=option_list))
  return(args)
} # end parseAarguments

main(parseArguments()) # call main method

# input_list = list()
# 
# input_list['count_log_norm_path'] = '/mnt/htcf_scratch/chasem/misc_scripts/deseq_model/results/libraryprotocol_librarydate_genotype_model/log2_norm_counts.csv'
# input_list['residual_log_norm_path'] = '/mnt/htcf_scratch/chasem/misc_scripts/deseq_model/results/libraryprotocol_librarydate_genotype_model/log2_norm_space_residuals.csv'
# input_list['protein_coding_gene_path'] = '~/code/cmatkhan/misc_scripts/deseq_model/data/KN99_protein_coding_gene_list_with_markers.csv'
# input_list['num_unwanted_variation_columns'] = 10
# input_list['output_directory'] = '~/code/cmatkhan/misc_scripts/deseq_model/results'
# input_list['output_name'] = 'unwated_variation_test'
# 
# main(input_list)