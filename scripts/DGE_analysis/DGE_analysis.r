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

get_git_version <- function(){
  sys.args <- commandArgs(trailingOnly = FALSE)
  execpth <- dirname(gsub('--file=',"",sys.args[grepl('--file=',sys.args, perl = TRUE)],ignore.case = FALSE))
  cmd <- sprintf("cd %s; git describe --tags; cd %s", execpth, getwd())
  version <- system(command = cmd, intern = TRUE)
  return(version)
}

save_tmp <- function(count_matrix, feature_matrix, output_dir){
  flog.info('save temporary feature matrix')
  write.table(feature_matrix,file = paste(c(dirname(output_dir),'tmp.feature.file.tsv'),collapse="/"), sep = "\t")
  flog.info("save temporary count matrix")
  write.table(count_matrix,file = paste(c(dirname(output_dir),'tmp.count.file.tsv'),collapse="/"), sep = "\t")
}

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

clean_matrix <- function(mat, min_sample = 0.05, min_gene = 100, remove_ambigious){
  #if matrix is formatted n_samples x n_genes
  keep_samples <- rowMeans(mat !=0) >= min_sample
  keep_genes <- colSums(mat) >= min_gene
  mat <- mat[keep_samples,keep_genes]
  if (remove_ambigious){
    n_before <- dim(mat)[2]
    mat <-mat[,!grepl('.*ambig.*',colnames(mat))]
    n_after <- dim(mat)[2]
    flog.info(sprintf("Removed %d ambigious genes out of %d total", (n_before-n_after), n_before))
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
                              remove_ambigious,
                              outdir,
                              output_dir){
  
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
    
    #TODO: Fix this ugly shit but works
    patient <- unlist(strsplit(tag,'_'))[1]
    fmat['patient'] <- patient
    if (k > 1) {
      if (patient %in% feature_matrix[['patient']]) {
         repl <- repl + 1
      } else {
        repl <- 1
      }
    } else {
      repl <- 1
    }
    print(repl)
    fmat['replicate'] <- paste(c("R",repl),collapse="")
   # fmat['replicate'] <- tag #unsure about best way to do this
    
    count_matrix <- bind_rows(count_matrix,cmat)
    feature_matrix <- bind_rows(feature_matrix,fmat)
    
  #save results temporarily only used in larger runs
   if (k %% 20 == 0) {
      save_tmp(count_matrix,feature_files,output_dir)
      }
    }
  
  rownames(count_matrix) <- c(1:dim(count_matrix)[1])
  rownames(feature_matrix) <- c(1:dim(feature_matrix)[1])
  
  count_matrix[is.na(count_matrix)] <- 0
  count_matrix <- clean_matrix(count_matrix, remove_ambigious = remove_ambigious)
  feature_matrix<-feature_matrix[rownames(count_matrix),]
  feature_matrix['replicate'] <- as.factor(feature_matrix[['replicate']])
  feature_matrix['patient'] <- as.factor(feature_matrix[['patient']])
  feature_matrix[feature_name] <- as.factor(feature_matrix[[feature_name]]) 
  
  count_matrix <-t(count_matrix)
  save_tmp(count_matrix, feature_matrix,output_dir)
  
  return(list(count_matrix = count_matrix, feature_matrix = feature_matrix))
}


DESeq_pipline <- function(count_matrix, feature_matrix, design_formula) {
  #' Perform DESeq2 DGE analysis with prepared matrices
  
  #use scran as to avoid error with zero counts
  scran_size_factors <- computeSumFactors(count_matrix, positive = TRUE)
  #construct SE object for compatibility
  se <- SummarizedExperiment(assays = list(counts = count_matrix), colData = feature_matrix)
  #release memory
  remove(count_matrix)
  remove(feature_matrix)
  
  #construct 
  mm <- model.matrix(design_formula, colData(se))
  mm <- mm[,colSums(mm) != 0]
  print(colSums(mm))
  dds <- DESeqDataSet(se,design = mm)
  sizeFactors(dds) <- scran_size_factors
  dds$tumor <- relevel(dds$tumor, "non")
  print("came here")
  dds<-DESeq(dds, full = mm)
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
                                remove_ambigious = remove_ambigious,
                                output_dir = output_dir)
  
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

if (args$debug){
  flog.threshold(DEBUG)
}else {
  flog.threshold(INFO)
}

flog.info(sprintf("Using github version %s", get_git_version()))
flog.info("Program iniated with argumens")
flog.info(args)
flog.info("__Starting DGE analysis__")
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


