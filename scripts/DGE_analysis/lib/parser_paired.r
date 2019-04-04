library(optparse)


make_parser <- function(parser,tag) {
  
  
  parser <- add_option(parser,
                       c("-c","--count_file"),
                       default = c(),
                       help = paste(c("Path to count file(s).",
                                      "Three options available",
                                      "1. Single directory : all files in",
                                      "directory will be loaded.",
                                      "2. Single File Path : File will be loaded",
                                      "3. Multiple Paths : All files will be loaded"
                                      ),
                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-f","--feature_file"),
                       default = c(),
                       help = paste(c("Path to feature file(s)",
                                      "Note how order of files must match",
                                      "with that of the count files",
                                      "Same options as for count_file",
                                      "argument exists"
                                      ),
                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-o", "--output_dir"),
                       default = getwd(),
                       type = "character",
                       help = paste(c("directory to save output into.",
                                      "if none specified cwd will be used."),
                                    collapse = " "))
  
parser <- add_option(parser,
                     c("-m", "--method"),
                     default = 'deseq2',
                     type = "character",
                     help = paste(c("Choose from available",
                                    "DGE-methods. Currently supported",
                                    'are "edgeR" and "deseq2" '),
                                    collapse = " ")
                      )

#  parser <- add_option(parser,
#                       c("-s", "--stamp"),
#                       default = '',
#                       type = "character",
#                       help = paste(c("stamp to add to output",
#                                      "can be any identifier. Appended to",
#                                      "standard-output-filename."),
#                                    collapse = " "))
  
  parser <- add_option(parser,
                       c("-ft","--filter_tumors"),
                       default = FALSE,
                       action = 'store_true',
                       help = paste(c("Only keep tumor spots.",
                                      "Will reduce memory requirements",
                                      "Recommended to use if comaptible with",
                                      "analysis."),collapse = " ")
  )
  
  parser <- add_option(parser,
                       c("-d","--design"),
                       default = -1,
                       help =paste(c("choose from available design formulas",
                                     "use flag without argument to see legend",
                                     "of defined formulas. Custom design formula",
                                     "can be provided via the --custom_design flag"),
                                   collapse = " ")
  )
  
  parser <- add_option(parser,
                       c('-cr','--contrast'),
                       default = NULL,
                       help = paste(c("contrast to to use in DESeq or ",
                                      'edgeR analysis.2',
                                      'must be on form "factor" "nominator"', 
                                      '"denominator"',
                                      "if a pre-defined design formula",
                                      "is used this will be ignored."
                       ),collapse = " ") 
  )
  
  parser <- add_option(parser,
                       c("-cd","--custom_design"),
                       default = NULL,
                       help =paste(c("choose from available design formulas",
                                     "use flag without argument to see legens",
                                     "of defined formulas. Custom design formula"),
                                   collapse = " ")
  )
  
    parser <- add_option(parser,
                       c("-z","--zero_coef"),
                       default = NA,
                       type = "integer",
                       help =paste(c("index of terms to set to",
                                     "zero in reduced formula, used",
                                     "in the LRT test. 1 based",
                                     "indexing is used."),
                                   collapse = " ")
  )
  
  parser <- add_option(parser,
                       c("-ssn","--subsample_number"),
                       default = 0,
                       type = "integer",
                       help =paste(c("integer giving the number of samples which", 
                                     "should be taken from each section provided",
                                     "the ratio between variables of the feature",
                                     "will be kept equal to the full set."),
                                   collapse = " ")
  )
  
  parser <- add_option(parser,
                       c("-ssf","--subsample_feature"),
                       default = "tumor",
                       type = "character",
                       help =paste(c("name of feature which subsampling",
                                     "should be performed w.r.t."),
                                   collapse = " ")
  )
  
  parser <- add_option(parser,
                       c("-a", "--remove_ambigious"),
                       default = FALSE,
                       action = "store_true",
                       help = paste(c("remove unambigious reads",
                                      "if not included all genes names",
                                      "will be kept."), collapse = " ")
  )
  
  parser <- add_option(parser,
                       c("-l", "--logname"),
                       default = paste(c(tag,".log"),collapse=""),
                       help = paste(c("full name of where to write",
                                      "log file. If none specified default",
                                      "with date and timestamp will be used"),collapse = " ")
  )
  
  parser <- add_option(parser,
                       c("-r", "--rsave"),
                       action = "store_true",
                       default = FALSE,
                       help = paste(c("Save R objects of",
                                      "DGE analysis results. Will be saved in",
                                      "output directory"), collapse = " ")
  )
    
  parser <- add_option(parser,
                       c('-fc','--fancy'),
                       default = FALSE,
                       action = "store_true",
                       help = paste(c("Save additional file where NA padj",
                                      "have been filtered out and filtering",
                                      "of genes based on padj value is perfomed.",
                                      "Also include Gene Symbols"),
                                      collapse = " ") 
  )
  
  parser <- add_option(parser,
                       c('-p','--pval'),
                       default = 0.01,
                       help = paste(c("Include if Gene Symbols",
                                      "for top n-genes should be",
                                      "included in ouput"),collapse = " ") 
  )
  
  parser <- add_option(parser,
                       c('-mps','--min_per_spot'),
                       type = "integer",
                       default = 20,
                       help = paste(c("Threshold for number",
                                      "of transcripts required",
                                      "to be present within a spot",
                                      "in order to keep this in analysis.",
                                      "Set to zero to keep all spots"
                                      ),collapse = " ") 
  )
  
  parser <- add_option(parser,
                       c('-mpg','--min_per_gene'),
                       default = 0.0,
                       help = paste(c("Threshold for number",
                                      "of transcripts required",
                                      "to have been mapped to a gene",
                                      "in order to keep this in analysis.",
                                      "Set to zero to keep all spots"
                       ),collapse = " ") 
  )
  
  parser <- add_option(parser,
                       c('-v','--volcano'),
                       default = FALSE,
                       action = "store_true",
                       help = paste(c("generate a volcano",
                                      "plot based on the results",
                                      "from the analysis."
                       ),collapse = " ") 
  )
  
  parser <- add_option(parser,
                       c('-zd','--zone_distance'),
                       default = NA,
                       type = "double",
                       help = paste(c("distance to nearest",
                                      "tumor spot in zone",
                                      "generation,"
                       ),collapse = " ") 
  )
  
  parser <- add_option(parser,
                       c('-q','--zone_method'),
                       default = 'three_levels',
                       help = paste(c("method to use when",
                                      "generating zones"
                       ),collapse = " ") 
  )
  
  parser <- add_option(parser,
                       c('-cnd','--condition_on'),
                       default = NA,
                       help = paste(c("conditon analysis on",
                                      "certain spots. Needs to",
                                      "Two arguments required",
                                      "[1] name of column containing covariates",
                                      "[2] covariate to condition on"
                       ),collapse = " ") 
  )
  
  
  
  return(parser)

}