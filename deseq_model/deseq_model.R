#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(DESeq2))
suppressMessages(library(tidyverse))
#suppressMessages(library(factoextra))
suppressMessages(library(BiocParallel))
register(MulticoreParam(10))

main = function(parsed_cmd_line_args){
  # main method of script
  
  print('...Parsing cmd line arguments')
  raw_counts_df_path = parsed_cmd_line_args$raw_counts
  metadata_df_path = parsed_cmd_line_args$metadata
  ruvr_unwanted_covariation_path = parsed_cmd_line_args$ruvr_unwanted_covaration_path
  num_unwanted_covariates = as.double(parsed_cmd_line_args$num_unwanted_covariates)
  protein_coding_gene_path = parsed_cmd_line_args$protein_coding_gene_path
  genotype_results_flag = parsed_cmd_line_args$genotype_results_flag
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
  model_matrix = createModelMatrix(design_formula, metadata_df)
  if (!(is.null(ruvr_unwanted_covariation_path) | is.null(num_unwanted_covariates))){
    
    unwanted_covariate_matrix = as.matrix(read_csv(ruvr_unwanted_covariation_path))
    model_matrix = cbind(model_matrix, unwanted_covariate_matrix[,1:num_unwanted_covariates])
   
  }
  writeOutDataframe(output_path, 'model_matrix', as_tibble(model_matrix))
  
  print('...construct deseq model')
  dds = createDeseqDataObject(raw_counts_df, metadata_df, model_matrix)
  deseq_model = generateDeseqModel(dds)
  writeOutDataframe(output_path, 'size_factors', as_tibble(sizeFactors(deseq_model))) # added 20200812
  
  print('...extracting coefficient matrix')
  coefficient_df = coef(deseq_model)
  rownames(coefficient_df) = rownames(raw_counts_df)
  coefficient_df = as_tibble(coefficient_df)
  writeOutDataframe(output_path, 'coefficients', coefficient_df)
  
  print('...calculating model predictions')
  model_predictions = calculateModelPredictions(model_matrix, coefficient_df, rownames(raw_counts_df), metadata_df$FASTQFILENAME)
  writeOutDataframe(output_path, 'model_predictions', as_tibble(model_predictions))  # added 20200812
  
  if (!is.null(genotype_results_flag)){
    if (genotype_results_flag == TRUE ){
      print('...writing out genotype results')
      writeGenotypeResults(deseq_model, output_path)
    }
  }
  
  print('...adding a pseudocount +1 and taking log2 of normalized counts')
  norm_counts_plus_pseudo = counts(deseq_model, normalized=TRUE) + 1
  log2_norm_counts = log2(norm_counts_plus_pseudo)
  writeOutDataframe(output_path, 'log2_norm_counts', as_tibble(log2_norm_counts)) # added 20200812
  
  # calculate residuals
  residual_log_norm_space_df = log2_norm_counts - model_predictions
  residual_log_norm_space_df = as_tibble(residual_log_norm_space_df)
  residual_log_norm_space_df[is.na(residual_log_norm_space_df)] = 0 # set NA to 0
  writeOutDataframe(output_path, 'log2_norm_space_residuals', residual_log_norm_space_df)
  
  residual_norm_space_df = norm_counts_plus_pseudo - apply(model_predictions, 2, function(x) 2**x) # unlog the model predictions
  residual_norm_space_df = as_tibble(residual_norm_space_df)
  residual_norm_space_df[is.na(residual_norm_space_df)] = 0 # set NA to 0
  writeOutDataframe(output_path, 'normalized_space_residuals', residual_norm_space_df)
  
  # return residuals to raw count scale
  unlogged_unnormalized_residual_df = round(as_tibble(unNormalize(residual_norm_space_df, sizeFactors(deseq_model)))) # NOTE: DESEQ only accepts ints. maybe able to submit directly to nbinomwald
  writeOutDataframe(output_path, 'raw_space_residuals', unlogged_unnormalized_residual_df)
  
  # calculate r_squared
  print('...calculating correlation coefficient')
  tss = calculateTotalSumOfSquares(norm_counts_plus_pseudo)
  rss = calculateResidualSumOfSquares(residual_norm_space_df)
  r_squared = 1 - (rss/tss)
  r_squared_df = tibble(name=output_name, r_2=r_squared)
  writeOutDataframe(output_path, 'r_squared', r_squared_df)

  # calculate principal components on residuals
  normalized_space_residual_pc = calculatePrincipalComponents(residual_norm_space_df, metadata_df)
  writeOutDataframe(output_path, 'normalized_space_residual_pc', normalized_space_residual_pc)
  
  log2_norm_space_residual_pc = calculatePrincipalComponents(residual_log_norm_space_df, metadata_df)
  writeOutDataframe(output_path, 'log2_normalized_space_residual_pc', log2_norm_space_residual_pc)
  
  # plot
  createResidualPcaPlots(log2_norm_space_residual_pc, design_formula, output_path, output_name)
  createDeseqPcaPlots(deseq_model, output_path)
  print(paste0('Done with ', as.character(design_formula)[2]))
  
} # end main()

convertDesignFormulaToColumnList = function(design_formula){
  
  # extract columns from the design formula as a list
  formula_str = as.character(design_formula)
  column_list = str_split(formula_str, '\\+')
  column_list = lapply(column_list, trimws)[-1]
  
  return(column_list)

} # end convertDesignFormulaToColumnList()

factorFormulaColumnsInMetadata = function(design_formula, df){
  
  # extract columns from the design formula as a list
  column_list = convertDesignFormulaToColumnList(design_formula)
  df[unlist(column_list)] = lapply(df[unlist(column_list)], factor)
  
  if ('LIBRARYPROTOCOL' %in% unlist(column_list)){
    df$LIBRARYPROTOCOL = relevel(df$LIBRARYPROTOCOL, ref='E7420L')
  }
  if ('LIBRARYDATE' %in% unlist(column_list)){
    df$LIBRARYDATE = relevel(df$LIBRARYDATE, ref=max(levels(df$LIBRARYDATE)))
  }
  
  return(df)

} # end factorFormulaColumnsInMetadata

createModelMatrix = function(design_formula, metadata_df){
  
  column_list = convertDesignFormulaToColumnList(design_formula)
  model_matrix = model.matrix(design_formula, metadata_df)
  
  libraryprotocol_librarydate_flag = libraryProtocolDateTest(column_list)
  
  if (libraryprotocol_librarydate_flag){
    col_to_drop = paste0('LIBRARYDATE', min(levels(metadata_df$LIBRARYDATE)))
    col_to_drop_index = match(col_to_drop, colnames(model_matrix))
    model_matrix = model_matrix[,-col_to_drop_index]
  }
  
  return(model_matrix)

} # end createModelMatrix()

libraryProtocolDateTest = function(column_list){
  
  library_date_flag = FALSE
  if ('LIBRARYDATE' %in% unlist(column_list)) {
    library_date_flag = TRUE
  }
  
  library_protocol_date_flag = FALSE
  if (library_date_flag == TRUE){
    if ('LIBRARYPROTOCOL' %in% unlist(column_list)){
      library_protocol_date_flag = TRUE
    }
  }
  
  return(library_protocol_date_flag)
  
} # end libraryProtocolDateTest()

createDeseqDataObject = function(raw_count_df, metadata_df, model_matrix){
  
  # generate deseq dataset (summarized experiment object) from count, metadata and design formula
  dds = DESeqDataSetFromMatrix(countData = raw_count_df, colData = metadata_df, design = model_matrix)
  
  return(dds)
  
} # end createDeseqDataObject

generateDeseqModel = function(dds){
  
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

calculatePrincipalComponents = function(residuals_df, metadata_df){
  
  # calculate principal components on residuals after log2 transformation (add psueocount)
  residuals_prcomp_object = prcomp(residuals_df)
  # convert to tibble
  residuals_pc_df = as_tibble(residuals_prcomp_object$rotation)
  # add rownames
  residuals_pc_df$FASTQFILENAME = rownames(residuals_prcomp_object$rotation)
  # join with metadata
  residuals_pc_df = dplyr::inner_join(residuals_pc_df, metadata_df, on=FASTQFILENAME)
  
  return(residuals_pc_df)
  
} # end calculatePrincipalComponents

edgerResidualCalculation = function(raw_counts, size_factors, model_matrix){
  
  # create DGEList object
  edge_r_object = y <- DGEList(counts=raw_counts, lib.siz = colSums(raw_counts), norm.factors = size_factors)
  # estimate common, trended and tagwise dispersions (see page 21 of edgeRUsersGuide())
  edge_r_object = estimateDisp(edge_r_object, model_matrix)
  # fit model
  edge_r_model = glmFit(edge_r_object, model_matrix)
  # extract deviance
  deviance_residuals = residuals(edge_r_model, type="deviance")
  
  return(deviance_residuals)
  
}

createResidualPcaPlots = function(log_space_residual_pc_df, design_formula, output_path, output_name){
  # TODO: GENERALIZE TO CREATE LIST OF PLOTS -- PASS LIST OF COLUMN VARIABLES IN TO PLOT
  # note: design_formula passed as a formula object
  
  graph_title = paste0(as.character(design_formula)[2], '_', output_name)
  
  g_librarydate = ggplot(log_space_residual_pc_df, aes(PC1,PC2))+geom_point(aes(color=LIBRARYDATE))+ggtitle(graph_title)
  g_libraryprotocol = ggplot(log_space_residual_pc_df, aes(PC1,PC2))+geom_point(aes(color=LIBRARYPROTOCOL))+ggtitle(graph_title)
  g_genotype = ggplot(log_space_residual_pc_df, aes(PC1,PC2))+geom_point(aes(color=GENOTYPE, size=GENOTYPE=='CNAG_00000'), alpha=.5)+ggtitle(graph_title)+theme(legend.position = "none")
  
  library_date_output = paste(output_path, 'residual_pca_by_library_date.pdf', sep='/')
  print('writing library_date PCA')
  ggsave(filename = library_date_output, plot = g_librarydate, device='pdf', height=8, width=12)
  library_protocol_output = paste(output_path, 'residual_pca_by_library_protocol.pdf', sep='/')
  print('writing library_protocol PCA')
  ggsave(filename = library_protocol_output, plot = g_libraryprotocol, device='pdf', height=8, width=12)
  genotype_output = paste(output_path, 'residual_pca_by_genotype.pdf', sep='/')
  print('writing genotype PCA')
  ggsave(filename = genotype_output, plot = g_genotype, device='pdf', height=8, width=12)
  
} # end createPcaPlots()

createDeseqPcaPlots = function(deseq_model, output_path){
  
  deseq_model_vst = vst(deseq_model)
  
  pca_genotype_output_path = paste(output_path, paste0('pca_genotype_vst.pdf'), sep='/')
  pdf(pca_genotype_output_path)
    pca_g = plotPCA(deseq_model_vst, intgroup='GENOTYPE') + coord_fixed() + 
      geom_point(aes(size=GENOTYPE=='CNAG_00000'), alpha = .05) + theme(legend.position = "none") + ggtitle('genotype (with vst)')
    plot(pca_g)
  dev.off()
  
  pca_librarydate_output_path = paste(output_path, paste0('pca_library_date_vst.pdf'), sep='/')
  pdf(pca_librarydate_output_path)
    pca_g = plotPCA(deseq_model_vst, intgroup='LIBRARYDATE') + coord_fixed() + ggtitle('library date (with vst)')
    plot(pca_g)
  dev.off()
  
  pca_libraryprotocol_output_path = paste(output_path, paste0('pca_library_protocol_vst.pdf'), sep='/')
  pdf(pca_libraryprotocol_output_path)
    pca_g = plotPCA(deseq_model_vst, intgroup='LIBRARYPROTOCOL') + coord_fixed() + ggtitle('library protocol (with vst)')
    plot(pca_g)
  dev.off()
  
} # createDeseqPcaPlots()

createScreePlot = function(prcomp_object, design_formula, output_path, output_name){
  
  graph_title = paste0(as.character(design_formula)[2], '_', output_name)
  
  scree_plot_output = paste(output_path, 'scree_plot.pdf', sep='/')
  pdf(scree_plot_output)
  plot(fviz_eig(prcomp_object), main = graph_title)
  dev.off() # this library unavailable on cluster in rnaseq_pipeline env (collision w/ some other package)
  
} # end createScreePlot()

writeGenotypeResults = function(deseq_model, output_path){
  # write out results for for each genotype against all wildtype, and by library protocol
  # if date flag is true, also include contrast against each genotype by date
  
  de_results_directory = paste(output_path, 'de_results', sep='/')
  dir.create(de_results_directory)
  
  # for genotype in resultsNames
  genotype_list = resultsNames(deseq_model)[startsWith(resultsNames(deseq_model), 'GENOTYPE')]
  for (genotype in genotype_list){
    
    # create results table
    results_table_output_path = paste(de_results_directory, paste0('de_results_', genotype, '.csv'), sep='/')
    results_table = results(deseq_model, name = genotype)
    write_csv(as_tibble(results_table), results_table_output_path)
    
    # plot ma
    ma_plot_output_path = paste(de_results_directory, paste0('ma_plot_', genotype, '.pdf'), sep='/')
    pdf(ma_plot_output_path)
        DESeq2::plotMA(results_table, ylim = c(-6,6), main=genotype)
    dev.off()
    
    pval_hist_output_path = paste(de_results_directory, paste0('pval_hist_', genotype, '.pdf'), sep='/')
    pdf(pval_hist_output_path)
        g = ggplot(as(results_table, "data.frame"), aes(x = pvalue)) +
          geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
        plot(g)
    dev.off()
    
  }
  
} # end writeGenotypeResults()

# writeeResults = function(model_matrix, deseq_model, output_path){
#   # example method, extrating genotype and creating contast vector
#   
#   # for genotype in resultsNames
#   genotype_list = resultsNames(deseq_model)[startsWith(resultsNames(deseq_model), 'GENOTYPE')]
#   for (genotype in genotype_list){
#     genotype_matrix_model = model_matrix[model_matrix[, genotype] == 1,]
#     coefficient_model_vector = as.numeric(apply(genotype_matrix_model, 2,
#                                                 function(genotype_model_matrix_column) 
#                                                   ifelse(1 %in% genotype_model_matrix_column, 1, 0)))
#     results_table = as.tibble(results(deseq_model, contrast = coefficient_model_vector, tidy=TRUE))
#     output_full_path = paste0(output_path, '/', 'results', '/', 'de_results_lib_prep_date_', genotype, '.csv')
#     write_csv(results_table, output_full_path)
#   }
#   
# }

# notes = function(){
#   
#   x = La.svd(as.matrix(residual_log_norm_space_df), nu=6969, nv=196)
#   diag(sing_d) = x$d
#   binder = matrix(rep(0, (6969-196) * 196), nrow=6969-196, ncol=196)
#   rbind(sing_d, binder)
#   sing_d_full = rbind(sing_d, binder)
#   w = x$u %*% sing_d_full
#   
# }

parseArguments <- function() {
  # parse and return cmd line input
  
  option_list <- list(
    make_option(c('-r', '--raw_counts'),
                help='raw count matrix (genes x samples)'),
    make_option(c('-m', '--metadata'), 
                help='metadata with all samples corresponding to the columns of the count data x metadata. must include the columns in the design formula'),
    make_option(c('-u', '--ruvr_unwanted_covaration_path'),
                help='path to sheet containing sample x k columns of unwanted variation'),
    make_option(c('-k', '--num_unwanted_covariates'),
                help='number of unwanted covariates to include in the model'),
    make_option(c('-g', '--genotype_results_flag'), action='store_true',
                help='set -g (no input) to write out all genotype results to subdiretory of results directory'),
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

#for testing
# input_list = list()
# input_list['raw_counts'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/data/test_2_counts.csv'
# input_list['metadata'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/data/test_2_metadata.csv'
# input_list['design_formula'] = '~LIBRARYPROTOCOL+LIBRARYDATE+GENOTYPE'
# input_list['ruvr_unwanted_covaration_path'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/results/fullrank_test//unwanted_variation.csv'
# input_list['num_unwanted_covariates'] = 3
# input_list['genotype_results_flag'] = TRUE
# input_list['output_directory'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/results'
# input_list['name'] = 'fullrank_test_3'
# 
# main(input_list)