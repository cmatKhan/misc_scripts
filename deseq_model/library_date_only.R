#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(DESeq2))
suppressMessages(library(tidyverse))
suppressMessages(library(BiocParallel))
register(MulticoreParam(10))

main = function(args){
  # main method of script
  
  print('...Parsing cmd line arguments')
  parsed_cmd_line_args = args
  raw_counts_df_path = parsed_cmd_line_args$raw_counts
  metadata_df_path = parsed_cmd_line_args$metadata
  coefficients_output = parsed_cmd_line_args$coefficients_output_full_path
  residuals_output = parsed_cmd_line_args$residuals_output_full_path
  
  print('...reading in raw counts')
  raw_counts_df = read_csv(raw_counts_df_path)
  print('...reading in metdata')
  metadata_df = read_csv(metadata_df_path)
  
  print('...construct deseq model')
  design_formula = '~LIBRARYDATE'
  deseq_model = generateDeseqModel(raw_counts_df, metadata_df, formula(design_formula))
  
  print('...extracting the model coefficients')
  coef_df = coef(deseq_model)
  write_csv(coef_df, coefficients_output)
  
  print('...adding a pseudocount +1 and taking log2 of normalized counts')
  log2_norm_counts = log2(counts(deseq_model, normalized=TRUE) + 1)
  
  # create a copy of log2_norm_counts to store the computed residuals
  residual_df = as_tibble(cbind(log2_norm_counts))
  
  for (i in seq(1,ncol(residual_df))){
    
    # extract sample name from column
    sample = colnames(residual_df)[i]
    
    # calculate the intercept + common variable (in this case, libraryprotocol)
    model_prediction = as.integer(coef_df[, 'Intercept']) + as.integer(coef_df[, 'LIBRARYDATE'])
    
    residual_df[,i] = abs(log2_norm_counts[,i] - model_prediction)
  }
  
  print(paste0('writing to output: ',output_path ))
  write_csv(residual_df, residuals_output)
  
} # end main()



generateDeseqModel = function(raw_count_df, metadata_df, design_formula){
  
  # generate deseq dataset (summarized experiment object) from count, metadata and design formula
  dds = DESeqDataSetFromMatrix(countData = raw_count_df, colData = metadata_df, design = design_formula)
  
  # construct the deseq model
  deseq_model = DESeq(dds, parallel=TRUE)
  
  return(deseq_model)
}

parseArguments <- function() {
  # parse and return cmd line input
  
  option_list <- list(
    make_option(c('-r', '--raw_counts'),
                help='raw count matrix (genes x samples).'),
    make_option(c('-m', '--metadata'),
                help='metadata with all samples corresponding to the columns of the count data x metadata. must include the columns in the design formula'),
    make_option(c('-a', '--coefficients_output_full_path'),
                help='full output path eg /scratch/mblab/chasem/rnaseq_pipeline/experiments/myexperiment_residuals.csv'),
    make_option(c('-b', '--residuals_output_full_path'),
                help='full output path eg /scratch/mblab/chasem/rnaseq_pipeline/experiments/myexperiment_residuals.csv'))
  
  args <- parse_args(OptionParser(option_list=option_list))
  return(args)
} # end parseAarguments

main(parseArguments()) # call main method
