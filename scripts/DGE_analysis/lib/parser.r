library(optparse)
make_parser <- function(parser) {
  parser <- add_option(parser,
                       c("-i","--count_dir"),
                       help = paste(c("directory of count-matrices.",
                                      "matrix name should be on form",
                                      '"count_data-ID_Replicate.tsv"'),
                                    collapse = " "))
  
  
  parser <- add_option(parser,
                       c("-o", "--output_dir"),
                       default = getwd(),
                       type = "character",
                       help = paste(c("directory to save output into.",
                                      "if none specified cwd will be used."),
                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-f","--feature_dir"),
                       help = paste(c("directory of feature-files.",
                                      "file name should be on form",
                                      '"count_data-ID_Replicate.tsv"'),
                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-s", "--select_for"),
                       default = NULL,
                       type = "character",
                       help = paste(c("file containing the sample ids.",
                                      "for those samples to be studied in",
                                      "the analysis"),
                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-g", "--gene_file"),
                       default =  NULL,
                       help = paste(c("list of genes to be included in DGE analysis.",
                                      "each row shall contain name of gene.",
                                      'if no list provided all genes will be used.'),
                                    collapse = " "))
  
  
  parser <- add_option(parser,
                       c("-d", "--design_file"),
                       help = paste(c("design file to be used in DGE-analysis.",
                                      "explictily states what features from the",
                                      "feature files that should be used."),
                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-a", "--remove_ambigious"),
                       default = FALSE,
                       action = "store_true",
                       help = paste(c("remove unambigious reads",
                                      "if not included all genes names",
                                      "will be kept."),
                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-fn", "--feature_name"),
                       default = "tumor",
                       help = paste(c("name of feature that will be studied.",
                                      "should be equivalent to column name in",
                                      "feature file."),
                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-md", "--max_dist"),
                       default = 2,
                       help = paste(c("maximum distance for a spot",
                                      "for two spots to be considered",
                                      "neighbours."),
                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-k", "--k_members"),
                       default = 20,
                       help = paste(c("number of spots to",
                                      "use for pooling of each sample"),
                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-n", "--n_samples"),
                       default = 20,
                       help = paste(c("number of samples to generate",
                                      "from the given section."),
                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-db", "--debug"),
                       default = FALSE,
                       action = "store_true",
                       help = paste(c("set logger level to debug",
                                      "otherwise info is used."),
                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-tn", "--top_n"),
                       default = 100,
                       type = "integer",
                       help = paste(c("save top n genes and p_adj values",
                                      "to file for easy analysis."),
                                    collapse = " "))
  
  #TODO : Implement option to run on multiple cores
  parser <- add_option(parser,
                       c("-w", "--workers"),
                       help = paste(c("number of cores/workers to be used.",
                                      "if non given half of maximum will be used.",
                                      "must be an integer"),
                                    collapse = " "),
                       type = "integer",
                       default = floor(detectCores()/2)
  )
  return(parser);
}
