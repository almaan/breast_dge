library("dplyr")
library("DESeq2")
library("scran")

library("docstring")
library("optparse")


extract_names <- function(string){
  #' Extracts sample id and replicate provided
  #' the count_data-ID_replicate.tsv filename
  
  sample_name <- strsplit(strsplit(string,'-')[[1]][2],'[.]')
  return(sample_name[[1]][1])
}

merge_matrices <- function(M1,M2){
  #' merge two matrices inplace
  r1_names <- rownames(M1)
  r2_names <- rownames(M2)
  M1 <- bind_rows(M1,M2)
  rownames(M1) <- c(r1_names,r2_names)
  M1[is.na(M1)] = 0.0
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
  full_pth <- paste(paste(feature_file_dir, paste('feature_data',sample_name,sep='-'), sep = '/'),'tsv', sep = '.')
  print(full_pth)
  pasFctr <- as.data.frame(read.csv(full_pth, sep = '\t', header = TRUE))
  part_1 <- apply(pasFctr[c("x_coord","y_coord")], 1, paste, collapse = "x")
  rownames(pasFctr) <- sapply(part_1, function(x) paste(sample_name,x,sep='_'))
  return(pasFctr)
}


suppressPackageStartupMessages(library("optparse"))


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
  
  print(sample_n(cnt,10))
  print(sample_n(coldata,10))
   
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
#                      required = FALSE,
#                      help = paste(c("list of genes to be included in DGE analysis.",
#                                     "each row shall contain name of gene.",
#                                     'if no list provided all genes will be used.'),
#                                   collapse = " "),
#                      )

#TODO: Implement design file parsing

# parser <- add_option(parser,
#                      c("-c", "--design"), 
#                      required = TRUE,
#                      help = paste(c("design file to be used in DGE-analysis.",
#                                     "explictily states what features that from the",
#                                     "feature files that should be used.",
#                                     "if non specified all features will be used"),
#                                   collapse = " "),
#                       )

args <- parse_args(parser)



#matirces for DESeq2 must be gene x sample
#print version from github git descrie 


#COUNT_DATA_DIR <- args$count_dir
#FEATURE_FILE_DIR = args$feature_dir


# 
# stp <- TRUE
# if (stp == TRUE) {
#   quit(0, status = 0)
# }

if (!interactive()) {
  
  COUNT_DATA_DIR <<- args$count_dir
  FEATURE_FILE_DIR <<- args$feature_dir
  main()
  #  SELECTED_GENES <- args$
} else {

  COUNT_DATA_DIR <- "/home/alma/ST-2018/CNNp/data/pre-tile-extraction"
  FEATURE_FILE_DIR = '/home/alma/ST-2018/CNNp/DGE/data/test_data'
  
  sample_names <- sort(list.files(COUNT_DATA_DIR, pattern = '*count*'))

  # For test only
  mat <- prepare_counts(COUNT_DATA_DIR,sample_names[1])
  cdata <- prepare_coldata(FEATURE_FILE_DIR, sample_names[1])
  
  
  sce <- SingleCellExperiment(list(counts = t(mat)))
  sce <- computeSumFactors(sce)
  sce = sce@int_colData@listData$size_factor
  
  # test
  dds <- DESeqDataSetFromMatrix(countData = t(mat),
                                colData = cdata,
                                design = ~ f1)
  sizeFactors(dds) <- sce
  keep <- rowSums(counts(dds)) >= 40
  dds <- dds[keep,]
  dds <- DESeq(dds, parallel = TRUE)
  res <-results(dds)
  res <- results(dds, name="f1_tumor_vs_non.tumor")
  res <- results(dds, contrast=c("f1","tumor","non-tumor"))
  resultsNames(dds)
  resLFC <- lfcShrink(dds, coef="f1_tumor_vs_non.tumor", type="apeglm")  
}
