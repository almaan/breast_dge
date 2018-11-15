suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("scran"))

suppressPackageStartupMessages(library("docstring"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("futile.logger"))



extract_names <- function(string){
  #' Extracts sample id and replicate provided
  #' the count_data-ID_replicate.tsv filename
  
  sample_name <- strsplit(strsplit(string,'-')[[1]][2],'[.]')
  return(sample_name[[1]][1])
}

merge_matrices <- function(M1,M2, na_fill = 0.0){
  #' merge two matrices inplace
  #' will expand matrix M1 as to include the union of columns from M1 and M2
  #' All unaccounted genes (not present in a matrix) will be
  #' assigned a value specified by the user
  #' 
  #' Input: 
  #'    - M1 : matrix with that rows of M2 will be vertical appended to
  #'    - M2 : matrix that is appended to M1 (vertically)
  #'    - na_fill : value to fill unspecified counts with. 0.0 is default
  #'    
  #' Output: 
  #'    - M1 : concatenated matrix
  
  r1_names <- rownames(M1)
  r2_names <- rownames(M2)
  M1 <- bind_rows(M1,M2)
  rownames(M1) <- c(r1_names,r2_names)
  M1[is.na(M1)] = na_fill
  return(M1)
}


prepare_counts <-function(main_pth,single_sample){
  #'
  full_pth <- paste(main_pth,single_sample,sep = "/")
  sample_name <- extract_names(single_sample)
  pasCts <- as.data.frame(read.csv(full_pth,sep = '\t',header = TRUE, row.names = 1))
  rownames(pasCts) <- sapply(rownames(pasCts), function(x) paste(sample_name,x, sep = '_'))
  return(pasCts)
}

prepare_coldata <- function(feature_file_dir,single_sample) {
  
  sample_name <- extract_names(single_sample)
  patient_id <- strsplit(sample_name,'-')[1]
  replicate <- strsplit(sample_name,'-')[2]
  
  full_pth <- paste(paste(feature_file_dir, paste('feature_data',sample_name,sep='-'), sep = '/'),'tsv', sep = '.')
  
  cdata <- as.data.frame(read.csv(full_pth, sep = '\t', header = TRUE))
  cdata['id'] <- patient_id
  cdata['replicate'] <- replicate
  
  stem_name <- apply(cdata[c("x_coord","y_coord")], 1, paste, collapse = "x")
  
  rownames(cdata) <- sapply(stem_name, function(x) paste(sample_name,x,sep='_'))

  return(cdata)
}


main <- function() {
  
  count_file_names <- sort(list.files(COUNT_DATA_DIR, pattern = '*count*'))
  
  for (num in 1:length(count_file_names)) {
    if (num == 1) {
      
      cnt <- prepare_counts(COUNT_DATA_DIR, count_file_names[num])
      coldata <- prepare_coldata(FEATURE_FILE_DIR, count_file_names[num])
      
    } else {
      cnt <- merge_matrices(cnt,prepare_counts(COUNT_DATA_DIR,count_file_names[num]))
      coldata <- merge_matrices(coldata, prepare_coldata(FEATURE_FILE_DIR, count_file_names[num]))
    }
  } 
  
  size_factors <- SingleCellExperiment(list(counts = t(cnt)))
  size_factors <- computeSumFactors(size_factors)
  size_factors = size_factors@int_colData@listData$size_factor

  dds <- DESeqDataSetFromMatrix(countData = t(cnt),
                                colData = coldata,
                                design = design)
  
  
  sizeFactors(dds) <- size_factors
  
  #for easy testing
  keep <- rowSums(counts(dds)) >= 40
  
  dds <- dds[keep,]
  dds <- DESeq(dds, 
               parallel = TRUE,
               BPPARAM = N_WORKERS)
  
  res <-results(dds)
  res <- results(dds,
                 contrast=c("f1","tumor","non-tumor"))
  
  resLFC <- lfcShri

  print(res)
     
}


parser <- OptionParser()
parser <- add_option(parser,
                     c("-i","--count_dir"),
                     help = paste(c("directory of count-matrices.",
                                    "matrix name should be on form",
                                    '"count_data-ID_Replicate.tsv"'),
                                  collapse = " "),
                    )

parser <- add_option(parser,
                     c("-f","--feature_dir"),
                     help = paste(c("directory of feature-files.",
                                    "file name should be on form",
                                    '"count_data-ID_Replicate.tsv"'),
                                  collapse = " "),
                      )

#TODO: Implement gene list parsing

# parser <- add_option(parser,
#                      c("-g", "--gene_list"),
#                      help = paste(c("list of genes to be included in DGE analysis.",
#                                     "each row shall contain name of gene.",
#                                     'if no list provided all genes will be used.'),
#                                   collapse = " "),
#                      )


parser <- add_option(parser,
                     c("-c", "--design"),
                     help = paste(c("design file to be used in DGE-analysis.",
                                    "explictily states what features from the",
                                    "feature files that should be used."),
                                  collapse = " "),
                      )

parser <- add_option(parser,
                     c("-w", "--workers"),
                     help = paste(c("number of cores/workers to be used.",
                                    "if non given half of maximum will be used.",
                                    "must be an integer"),
                                  collapse = " "),
                     default = floor(detectCores()/2)
)

args <- parse_args(parser)

source(args$design)

if (!exists(x = "design")) {
  print("Please provide a proper design-file")
  print("Exiting")
  quit(0,status = 0)
}



#matirces for DESeq2 must be gene x sample
#print version from github git descrie 


#COUNT_DATA_DIR <- args$count_dir
#FEATURE_FILE_DIR = args$feature_dir

#TODO : Consider wether default option for design should be includeded.

# 
# stp <- TRUE
# if (stp == TRUE) {
#   quit(0, status = 0)
# }

if (!interactive()) {
  
  COUNT_DATA_DIR <<- args$count_dir
  FEATURE_FILE_DIR <<- args$feature_dir
  N_WORKERS <<- args$workers
  main()
  #  SELECTED_GENES <- args$
} else {
  print("Run via terminal for proper activation.")
  
}
