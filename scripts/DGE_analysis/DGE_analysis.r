sh <- suppressPackageStartupMessages

sh(library(SummarizedExperiment))
sh(library(DESeq2))
sh(library(ggplot2))
sh(library(scran))
sh(library(stringr))
sh(library(optparse))
sh(library(dplyr))
sh(library(futile.logger))
source(paste(c(getwd(),"lib/parser.r"),collapse = "/"))
source(paste(c(getwd(),"lib/poolf.r"), collapse = "/"))

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

clean_matrix <- function(mat, min_spot = 100, min_gene = 0.05, remove_ambigious){
  keep_genes <- colMeans(mat !=0) >= min_gene
  keep_spots <- rowSums(mat) >= min_spot
  mat <- mat[keep_spots,keep_genes]
  if (remove_ambigious){
    n_before <- dim(mat)[2]
    mat <-mat[,!grepl('ambiguous',colnames(mat))]
    n_after <- dim(mat)[2]
    flog.info(sprintf("Removed %d genes out of %d total", (n_before-n_after), n_before))
  }
  return(mat)
}


generate_matrices <- function(count_pth,
                              feature_pth,
                              select_for,
                              feature_name,
                              k_members,
                              n_samples,
                              max_dist,
                              gene_file,
                              remove_ambigious){
  
  pattern_cnt <- "count.*tsv"
  pattern_fea <-"*tsv"
  count_files <- sort(list.files(count_pth,pattern = pattern_cnt))
  
  if (max_dist < 0) {
    max_dist <- Inf
  }
  
  if (!is.null(gene_file)){
    gene_names <- read.table(gene_file, sep = "\n", header = FALSE)
  }
  
  if (!is.null(select_for)) {
    selected_samples <- readLines(file(select_for,"r"))
    selected_samples <- selected_samples[sapply(selected_samples, function(x) x != "")]
    count_files <- count_files[sapply(count_files, function(x) any(sapply(selected_samples, function (y) grepl(y,x))))]
    
  }
  
  flog.info(sprintf("A total of %d sections with unique ID's will be used", length(count_files)))
  
  feature_files <- list.files(feature_pth, pattern = pattern_fea)
  count_matrix <- data.frame()
  feature_matrix <- data.frame()
  
  
  for(k in c(1:length(count_files))) {
    cmat <- load_matrix(count_pth,count_files[k])
    tag <- get_id(count_files[k])
    flog.info(sprintf("Sampling form sample %s", tag))
    
    if (!is.null(gene_file)){
      select_gene <- intersect(gene_names, colnames(cmat)) 
      cmat <- cmat[,gene_names]
    }
    
    fmat_name <- grep(paste(c(".*",tag,".*"),collapse = ""),feature_files)
    fmat <- load_features(feature_pth,feature_files[fmat_name])
    inter <- intersect(rownames(fmat),rownames(cmat))
    fmat <-fmat[inter,]
    cmat <- cmat[inter,]
    
    matl <- make_pseudo(cmat,fmat,select = feature_name, lim = max_dist, 
                        k_members = k_members,
                        n_samples = n_samples)
    
    cmat <- matl$pseudo_cnt
    fmat <- matl$pseudo_feat
    
        
    fmat['patient'] <- unlist(strsplit(tag,'_'))[1]
    #fmat['replicate'] <- unlist(strsplit(tag,'_'))[2]
    fmat['replicate'] <- tag #unsure about best way to do this
    
    count_matrix <- bind_rows(count_matrix,cmat)
    feature_matrix <- bind_rows(feature_matrix,fmat)
  }
  rownames(count_matrix) <- c(1:dim(count_matrix)[1])
  rownames(feature_matrix) <- c(1:dim(feature_matrix)[1])
    
  count_matrix[is.na(count_matrix)] <- 0
  count_matrix <- clean_matrix(count_matrix, remove_ambigious = remove_ambigious)
  
  feature_matrix<-feature_matrix[rownames(count_matrix),]
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

#TODO Make a results formatting file printing top genes. Top into one output folder and zip.

main <- function(count_input_dir,
                 feature_input_dir,
                 select_for,
                 output_dir,
                 design_file,
                 feature_name,
                 k_members,
                 n_samples,
                 max_dist,
                 gene_file,
                 remove_ambigious){
  
  design_formula <- readLines(file(design_file,"r"))
  design_formula <- as.formula(design_formula[!grepl('#',design_formula)])
  flog.info(paste(c("Using design : ", design_formula ),collapse = " "))
  close(file(design_file))
  
  matrices <- generate_matrices(count_pth = count_input_dir,
                                feature_pth = feature_input_dir,
                                select_for = select_for,
                                feature_name = feature_name,
                                k_members = k_members,
                                n_samples = n_samples,
                                max_dist = max_dist,
                                gene_file = gene_file,
                                remove_ambigious = remove_ambigious)
  
  dds <- DESeq_pipline(matrices$count_matrix, matrices$feature_matrix, design_formula = design_formula)
  res <- results(dds)
  write.csv(as.data.frame(res), file = paste(c(output_dir), collapse = "/"))
 
     
}

#snowparam <- SnowParam(workers = 4, type = "SOCK")
#register(snowparam, default = TRUE)
#registered()

parser <- OptionParser()
parser <- make_parser(parser)
args <- parse_args(parser)

flog.threshold(DEBUG)
flog.info("Starting DGE analysis")
flog.info(sprintf("Will be using %d neighbours for sampling generating %d samples per section", 
                  args$k_members, args$n_samples))

main(count_input_dir = args$count_dir,
     feature_input_dir = args$feature_dir,
     design_file = args$design_file,
     output_dir =args$output_dir,
     select_for = args$select_for,
     feature_name = args$feature_name,
     gene_file = args$gene_file,
     k_members = args$k_members,
     n_samples = args$n_samples,
     max_dist = args$max_dist,
     remove_ambigious = args$remove_ambiguous)


