#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(DESeq2))
suppressMessages(library(tidyverse))
#suppressMessages(library(factoextra))
suppressMessages(library(BiocParallel))
register(MulticoreParam(10))

main = function(args){
  # main method of script
  
  print('...Parsing cmd line arguments')
  parsed_cmd_line_args = args
  raw_counts_df_path = parsed_cmd_line_args$raw_counts
  metadata_df_path = parsed_cmd_line_args$metadata
  output_dir = parsed_cmd_line_args$output_directory
  output_name = parsed_cmd_line_args$name
  # error if no tilde in first position TODO!!
  design_formula = formula(parsed_cmd_line_args$design_formula)
  
  # create output directory
  output_path = paste(output_dir,output_name, sep='/')
  print(paste0('output directory created at: ', output_path))
  dir.create(output_path)
  
  print('...reading in raw counts')
  raw_counts_df = read_csv(raw_counts_df_path)
  
  print('...reading in metdata')
  metadata_df = read_csv(metadata_df_path)
  metadata_df$LIBRARYDATE = as.Date(metadata_df$LIBRARYDATE, format="%m.%d.%y")

  print('...factoring design formula columns')
  metadata_df = factorFormulaColumnsInMetadata(design_formula, metadata_df)
  
  print('...constructing design matrix from formula')
  model_matrix = model.matrix(design_formula, metadata_df)
  writeOutDataframe(output_path, 'model_matrix', as_tibble(model_matrix))
  
  print('...construct deseq model')
  deseq_model = generateDeseqModel(raw_counts_df, metadata_df, design_formula)
  writeOutDataframe(output_path, 'size_factors', as_tibble(sizeFactors(deseq_model))) # added 20200812
  
  print('...extracting coefficient matrix')
  coefficient_df = coef(deseq_model)
  rownames(coefficient_df) = rownames(raw_counts_df)
  coefficient_df = as_tibble(coefficient_df)
  writeOutDataframe(output_path, 'coefficients', coefficient_df)
  
  print('...calculating model predictions')
  model_predictions = calculateModelPredictions(model_matrix, coefficient_df, rownames(raw_counts_df), metadata_df$FASTQFILENAME)
  writeOutDataframe(output_path, 'model_predictions', as_tibble(model_predictions))  # added 20200812
  
  print('...adding a pseudocount +1 and taking log2 of normalized counts')
  norm_counts_plus_pseudo = counts(deseq_model, normalized=TRUE) + 1
  log2_norm_counts = log2(norm_counts_plus_pseudo)
  writeOutDataframe(output_path, 'log2_norm_counts', as_tibble(log2_norm_counts)) # added 20200812
  
  # calculate residuals
  residual_norm_space_df = norm_counts_plus_pseudo - apply(model_predictions, 2, function(x) 2**x) # unlog the model predictions
  residual_norm_space_df = as_tibble(residual_norm_space_df)
  residual_norm_space_df[is.na(residual_norm_space_df)] = 0 # set NA to 0
  writeOutDataframe(output_path, 'normalized_residuals', residual_norm_space_df)
  
  # return residuals to raw count scale
  unlogged_unnormalized_residual_df = round(as_tibble(unNormalize(residual_norm_space_df, sizeFactors(deseq_model)))) # NOTE: DESEQ only accepts ints. maybe able to submit directly to nbinomwald
  writeOutDataframe(output_path, 'raw_residuals', unlogged_unnormalized_residual_df)
  
  # calculate r_squared
  print('...calculating correlation coefficient')
  tss = calculateTotalSumOfSquares(norm_counts_plus_pseudo)
  rss = calculateResidualSumOfSquares(residual_norm_space_df)
  r_squared = 1 - (rss/tss)
  r_squared_df = tibble(name=output_name, r_2=r_squared)
  writeOutDataframe(output_path, 'r_squared', r_squared_df)

  # calculate principal components on residuals after log2 transformation (add psueocount)
  residuals_prcomp_object = prcomp(residual_norm_space_df)
  residuals_pc_df = as_tibble(residuals_prcomp_object$rotation)
  residuals_pc_df$FASTQFILENAME = rownames(residuals_prcomp_object$rotation)
  residuals_pc_df = dplyr::inner_join(residuals_pc_df, metadata_df, on=FASTQFILENAME)
  writeOutDataframe(output_path, 'log_normalized_residual_pc', residuals_pc_df)
  
  # plot
  createPlots(residuals_prcomp_object, residuals_pc_df, design_formula, output_path, output_name)
  print(paste0('Done with ', as.character(design_formula)[2]))
  
} # end main()

factorFormulaColumnsInMetadata = function(design_formula, df){
  
  # extract columns from the design formula as a list
  formula_str = as.character(design_formula)
  column_list = str_split(formula_str, '\\+')
  column_list = lapply(column_list, trimws)[-1]
  df[unlist(column_list)] = lapply(df[unlist(column_list)], factor)
  
  return(df)

}

generateDeseqModel = function(raw_count_df, metadata_df, design_formula){
  
  # generate deseq dataset (summarized experiment object) from count, metadata and design formula
  dds = DESeqDataSetFromMatrix(countData = raw_count_df, colData = metadata_df, design = design_formula)
  
  # construct the deseq model
  deseq_model = DESeq(dds, parallel=TRUE)
  
  return(deseq_model)
  
} # end generateDeseqModel()

calculateModelPredictions = function(model_matrix, coefficient_matrix, gene_list, sample_list){
  
  # model matrix is m=sample  n=model predictors, coefficient matrix is m=gene n=predictors
  y_hat = model_matrix %*% t(coefficient_matrix)
  # return to m=gene n=sample
  y_hat = t(y_hat)
  # name rows and columns
  rownames(y_hat) = gene_list
  colnames(y_hat) = sample_list
  
  return(y_hat)
  
} # end calculateModelPredictions()

calculateTotalSumOfSquares = function(log2_norm_count_df){
  
  deviation_df = log2_norm_count_df - rowMeans(log2_norm_count_df) # need to check operation here 20200812. problem in r2 calc
  
  squared_deviation_df = deviation_df**2
  squared_deviation_df[is.na(squared_deviation_df)] = 0
  
  return(sum(colSums(squared_deviation_df)))
  
} # end calculateTotalSumOfSquares()

calculateResidualSumOfSquares = function(residuals_df){
  
  squared_residuals = residuals_df**2
  squared_residuals[is.na(squared_residuals)] = 0
  
  return(sum(colSums(squared_residuals)))
  
} # end calculateResidualSumOfSquares()

unNormalize = function(norm_residuals_df, size_factors){

  # un-normalize the residuals
  unnormalized_unlogged_residuals = as_tibble(norm_residuals_df)
  # multiply normalized counts by size factor
  for (j in seq(1,length(size_factors))){
    unnormalized_unlogged_residuals[,j] = unnormalized_unlogged_residuals[,j] * size_factors[j]
  }

  return(unnormalized_unlogged_residuals)
  
} # end unNormalize()

writeOutDataframe = function(output_path, chart_name, df){
  
  # output path
  csv_output_path = paste(output_path, paste0(chart_name, '.csv'), sep='/')
  # tell user whats what  
  print(paste0('writing sheet: ', csv_output_path))
  # write
  write_csv(df, csv_output_path)
  
} # end writeOutDataframe()

createPlots = function(prcomp_object, residual_pc_df, design_formula, output_path, output_name){
  # TODO: GENERALIZE TO CREATE LIST OF PLOTS -- PASS LIST OF COLUMN VARIABLES IN TO PLOT
  # note: design_formula passed as a formula object
  
  graph_title = paste0(as.character(design_formula)[2], '_', output_name)
  
  g_librarydate = ggplot(residual_pc_df, aes(PC1,PC2))+geom_point(aes(color=LIBRARYDATE))+ggtitle(graph_title)
  g_libraryprotocol = ggplot(residual_pc_df, aes(PC1,PC2))+geom_point(aes(color=LIBRARYPROTOCOL))+ggtitle(graph_title)
  g_genotype = ggplot(residual_pc_df, aes(PC1,PC2))+geom_point(aes(color=GENOTYPE))+ggtitle(graph_title)+theme(legend.position = "none")
  
  library_date_output = paste(output_path, 'pca_by_library_date.pdf', sep='/')
  ggsave(filename = library_date_output, plot = g_librarydate, device='pdf', height=8, width=12)
  library_date_output = paste(output_path, 'pca_by_library_protocol.pdf', sep='/')
  ggsave(filename = library_date_output, plot = g_libraryprotocol, device='pdf', height=8, width=12)
  library_date_output = paste(output_path, 'pca_by_genotype.pdf', sep='/')
  ggsave(filename = library_date_output, plot = g_genotype, device='pdf', height=8, width=12)
  
  # scree_plot_output = paste(output_path, 'scree_plot.pdf', sep='/')
  # pdf(scree_plot_output)
  # plot(fviz_eig(prcomp_object), main = graph_title)
  # dev.off() # this library unavailable on cluster in rnaseq_pipeline env (collision w/ some other package)
  
} # end createPlots()

parseArguments <- function() {
  # parse and return cmd line input
  
  option_list <- list(
    make_option(c('-r', '--raw_counts'),
                help='raw count matrix (genes x samples)'),
    make_option(c('-m', '--metadata'), 
                help='metadata with all samples corresponding to the columns of the count data x metadata. must include the columns in the design formula'),
    make_option(c('-d', '--design_formula'), 
                help='eg ~LIBRARYDATE+GENOTYPE currently does not accept variables with continuous data. base level will be the first level in the factored column(s). Currently not set up for interaction terms'),
    make_option(c('-o', '--output_directory'), 
                help='path to directory to output results'),
    make_option(c('-n', '--name'),
                help='name of results subdirectory outputed in the path above eg if comparing library date libraryDate_model might be the name'))
  
  args <- parse_args(OptionParser(option_list=option_list))
  return(args)
} # end parseAarguments

main(parseArguments()) # call main method