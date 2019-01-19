library(optparse)

log_name <- function() {
  timestamp <- gsub(pattern = ':| ', replacement = "-", x = as.character(Sys.time())) 
  return(paste(c('DGE_analysis.',timestamp,'.',runif(1),".log"),collapse = ""))
}

make_parser <- function(parser) {
  
  
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
  
  #parser <- add_option(parser,
  #                     c('-s','--samples'),
  #                     default = c(),
  #                     help = paste(c("samples to be used in analysis",
  #                                    "if only patient ids and not replicates",
  #                                    "are provided all replicates of that patient",
  #                                    "will be used"),collapse = " ")
  #)
  
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
                       c("-cd","--custom_design"),
                       default = NULL,
                       help =paste(c("choose from available design formulas",
                                     "use flag without argument to see legens",
                                     "of defined formulas. Custom design formula"),
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
                       default = log_name(),
                       help = paste(c("full name of where to write",
                                      "log file. If none specified default",
                                      "with date and timestamp will be used"),collapse = " ")
  )
  
  parser <- add_option(parser,
                       c("-pr", "--paranoid"),
                       action = "store_true",
                       default = FALSE,
                       help = paste(c("use paranoid mode. R objects of",
                                      "DESeq2 analysis will be saved in",
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
  
  return(parser)

}