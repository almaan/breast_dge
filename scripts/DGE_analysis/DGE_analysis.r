library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(scran)
library(stringr)
library(optparse)
library(dplyr)

load_matrix <- function(main_pth, file_name, nrows = -1){
  mat <- read.csv(paste(c(main_pth,file_name),collapse="/"), sep = "\t", nrows = nrows, header = TRUE, row.names = 1)
  return(mat)
}
load_features <- function(main_pth,file_name) {
  feat <- read.csv(paste(c(main_pth,file_name),collapse="/"), sep = "\t",
                   header = TRUE, row.names = 1, stringsAsFactors = TRUE)
  return(feat)
}

get_id <- function(path) {
  return(str_extract(path,pattern = "\\d{4,5}_\\w{2}"))
}

clean_matrix <- function(mat, min_spot = 100, min_gene = 0.05){
  keep_genes <- colMeans(mat !=0) >= min_gene
  keep_spots <- rowSums(mat) >= min_spot
  mat <- mat[keep_spots,keep_genes]
  return(mat)
}


generate_matrices <- function(count_pth,feature_pth,select_for ){
  pattern_cnt <- "count.*tsv"
  pattern_fea <-"*tsv"
  count_files <- sort(list.files(count_pth,pattern = pattern_cnt))
  
  if (!is.null(select_for)) {
    selected_samples <- readLines(file(select_for,"r"))
    selected_samples <- selected_samples[sapply(selected_samples, function(x) x != "")]
    unlink(select_for)
    count_files <- count_files[sapply(count_files, function(x) get_id(x) %in% selected_samples)]
  }
  
  feature_files <- list.files(feature_pth, pattern = pattern_fea)
  count_matrix <- data.frame()
  feature_matrix <- data.frame()
  
  
  for(k in c(1:length(count_files))) {
    cmat <- load_matrix(count_pth,count_files[k])
    tag <- get_id(count_files[k])
    
    fmat_name <- grep(paste(c(".*",tag,".*"),collapse = ""),feature_files)
    fmat <- load_features(feature_pth,feature_files[fmat_name])
    inter <- intersect(rownames(fmat),rownames(cmat))
    fmat <-fmat[inter,]
    fmat[c('patient','replicate')] <- unlist(strsplit(tag,'_'))
    cmat <- cmat[inter,]
    count_matrix <- bind_rows(count_matrix,cmat)
    feature_matrix <- bind_rows(feature_matrix,fmat)
  }
  
  rownames(count_matrix) <- c(1:dim(count_matrix)[1])
  rownames(feature_matrix) <- c(1:dim(feature_matrix)[1])
  
  feature_matrix<-feature_matrix[rownames(count_matrix),]
  
    
  count_matrix[is.na(count_matrix)] <- 0
  count_matrix <- clean_matrix(count_matrix)
  count_matrix <-t(count_matrix)
  
  return(list(count_matrix = count_matrix, feature_matrix = feature_matrix))
}


DESeq_pipline <- function(count_matrix, feature_matrix, design_formula) {
 
  scran_size_factors <- computeSumFactors(count_matrix, positive = TRUE)
  se <- SummarizedExperiment(assays = count_matrix, colData = feature_matrix)
  remove(count_matrix)
  remove(feature_matrix)
  
  dds <- DESeqDataSet(se, design = design_formula)
  sizeFactors(dds) <- scran_size_factors
  dds$tumor <- relevel(dds$tumor, "tumor")
  dds<-DESeq(dds, 
             #parallel = TRUE,
             #BPPARAM = bpparam("SnowParam")
             )
  return(dds)
}

main <- function(count_input_dir,
                 feature_input_dir,
                 select_for,
                 output_dir,
                 design_file){
  
  design_formula <- as.formula(readLines(file(design_file,"r")))
  close(file(design_file))
  
  matrices <- generate_matrices(count_pth = count_input_dir, feature_pth = feature_input_dir,select_for)
  dds <- DESeq_pipline(matrices$count_matrix, matrices$feature_matrix, design_formula = design_formula)
  res <- results(dds)
  write.csv(as.data.frame(res), file = paste(c(output_dir,"DGE_analysis_result.tsv"), collapse = "/"), sep = "\t")
 
     
}

snowparam <- SnowParam(workers = 4, type = "SOCK")
register(snowparam, default = TRUE)
registered()

parser <- OptionParser()
parser <- add_option(parser,
                     c("-i","--count_dir"),
                     help = paste(c("directory of count-matrices.",
                                    "matrix name should be on form",
                                    '"count_data-ID_Replicate.tsv"'),
                                  collapse = " "),
)


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

#TODO - allow for selection of specific genes
#parser <- add_option(parser,
#                     c("-g", "--gene_file"),
#                     default =  NULL,
#                     help = paste(c("list of genes to be included in DGE analysis.",
#                                    "each row shall contain name of gene.",
#                                    'if no list provided all genes will be used.'),
#                                  collapse = " "))


parser <- add_option(parser,
                     c("-d", "--design_file"),
                     help = paste(c("design file to be used in DGE-analysis.",
                                    "explictily states what features from the",
                                    "feature files that should be used."),
                                  collapse = " "))

parser <- add_option(parser,
                     c("-w", "--workers"),
                     help = paste(c("number of cores/workers to be used.",
                                    "if non given half of maximum will be used.",
                                    "must be an integer"),
                                  collapse = " "),
                     type = "integer",
                     default = floor(detectCores()/2)
)

args <- parse_args(parser)

main(count_input_dir = args$count_dir,
     feature_input_dir = args$feature_dir,
     design_file = args$design_file,
     output_dir =args$design_file,
     select_for = args$select_for)


