sh <- suppressPackageStartupMessages

sh(library(SummarizedExperiment))
sh(library(DESeq2))
sh(library(ggplot2))
sh(library(scran))
sh(library(stringr))
sh(library(optparse))
sh(library(dplyr))
sh(library(futile.logger))

#load functions
source(paste(c(getwd(),"lib/parser.r"),collapse = "/"))
source(paste(c(getwd(),"lib/poolf.r"), collapse = "/"))
source(paste(c(getwd(),"lib/utils.r"), collapse ="/"))

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
  
  #correct max-distance in neighbour search to infinity if negative 
  if (max_dist < 0) {
    max_dist <- Inf
  }
  
  #check if list of genes to be used is specified
  if (!is.null(gene_file)){
    gene_names <- read.table(gene_file, sep = "\n", header = FALSE)
  }
  #patterns to use for grepping of count and feature matrices
  pattern_cnt <- "count.*tsv"
  pattern_fea <-"*tsv"
  #get all filenames of count matrices within input folder
  count_files <- sort(list.files(count_pth,pattern = pattern_cnt))
  feature_files <- list.files(feature_pth, pattern = pattern_fea)
  
  #select only subset of files specified to be analyzed in input folder
  #if non specified all will be analyzed
  if (!is.null(select_for) & file.exists(select_for)) {
    selected_samples <- readLines(file(select_for,"r"))
    selected_samples <- selected_samples[sapply(selected_samples, function(x) x != "")]
    count_files <- count_files[sapply(count_files, function(x) any(sapply(selected_samples, function (y) grepl(y,x))))]
    
  }
  
  flog.info(sprintf("A total of %d sections with unique ID's will be used", length(count_files)))
  
  #create objects
  count_matrix <- data.frame()
  feature_matrix <- data.frame()
  
  #loop over all count matrices to be analyzed
  tfull <- 0.0
  for(k in c(1:length(count_files))) {
    tstart <- Sys.time()
    cmat <- load_matrix(count_pth,count_files[k]) #load single count matrix
    tag <- get_id(count_files[k]) #get patient id of count matrix, used for feature file pairing
    eta <- ifelse(k>1,as.character(round(rtime,2)),'--')
    flog.info(sprintf("Pooling from sample %s : %d / %d | ETA : %smin", tag, k, length(count_files), eta))
    
    #select only specifed genes if such are given
    if (!is.null(gene_file)){
      select_gene <- intersect(gene_names, colnames(cmat)) 
      cmat <- cmat[,gene_names]
    }
    
    #load feature matrix and match rows with count matrix
    fmat_name <- grep(paste(c(".*",tag,".*"),collapse = ""),feature_files) #find matching feature file
    fmat <- load_features(feature_pth,feature_files[fmat_name]) #load matching feature file
    inter <- intersect(rownames(fmat),rownames(cmat)) #intersection of sample names between count and feature file
    fmat <-fmat[inter,] #only use samples present with both features and counts available
    cmat <- cmat[inter,] #only use samples present with both features and counts available
    
    #generate "pseudo samples" by pooling from original matrices
    matl <- make_pseudo(cmat,fmat,select = feature_name,
                        lim = max_dist, 
                        k_members = k_members,
                        n_samples = n_samples)
    
    #extract results for easy appending of properties
    cmat <- matl$pseudo_cnt
    fmat <- matl$pseudo_feat
    
    patient <- unlist(strsplit(tag,'_'))[1]
    #set patient and replicate for pseudo samples
    fmat['replicate'] <- k
    fmat['patient'] <- patient

    #for nesting of patients and replicates in DESeq2    
    if (ifelse(k > 1,patient %in% feature_matrix[['patient']],FALSE)) {
      #if not first iteration and seen patient
      repl <- repl + 1
    } else {
      #if first iteration or new patient
      repl <- 1
    }
    fmat['replicate.nested'] <- repl
    
    #add pseudo-samples to full count matrix and feature matrix
    count_matrix <- bind_rows(count_matrix,cmat)
    feature_matrix <- bind_rows(feature_matrix,fmat)
  
  #save results temporarily, only used in larger runs
  #TODO: once code works, add unlink in the end to remove tmp files
  if (k %% 20 == 0) {
      save_tmp(count_matrix,feature_files,output_dir)
  }
  #estimate remaining time
  tfull <- tfull + as.numeric((Sys.time() - tstart))/60
  rtime <- tfull/k * (length(count_files)-k)
  }
  
  #set rownames for easy manipulation
  rownames(count_matrix) <- c(1:dim(count_matrix)[1])
  rownames(feature_matrix) <- c(1:dim(feature_matrix)[1])
  
  #assign missing gene counts (not present) as 0 rather than na, from row-bind
  count_matrix[is.na(count_matrix)] <- 0
  #clean matrix as to remove low count samples and genes with very low abundancy
  #also remove genes that are ambigiously assigned if requested
  count_matrix <- clean_matrix(count_matrix, remove_ambigious = remove_ambigious)
  feature_matrix<-feature_matrix[rownames(count_matrix),]

  #treat features as factors for categorical analysis
  #feature_matrix['replicate'] <- as.factor(feature_matrix[['replicate']])
  feature_matrix['replicate'] <- as.factor(feature_matrix[['replicate']])
  feature_matrix['replicate.nested'] <- as.factor(feature_matrix[['replicate.nested']])
  feature_matrix['patient'] <- as.factor(feature_matrix[['patient']])
  
  feature_matrix[feature_name] <- as.factor(feature_matrix[[feature_name]])
  
  #transpose as to get proper orientation, generates n_genes x n_samples
  count_matrix <-t(count_matrix)
  #save matrix tuple temporarily
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
  #remove(count_matrix)
  #remove(feature_matrix)
  
  #construct model matrix using specified design formulae
  mm <- model.matrix(design_formula, colData(se))
  #identify elements where no interaction (all zero) are present
  all.zero <- apply(mm, 2, function(x) all(x==0))
  #get indices for zero columns
  zidx <- which(all.zero)
  if (length(zidx) > 0) {
   # if any zero columns, remove these for linear independency
    #recommended procedure by Love
    mm <- mm[,-zidx]
  }
  
  dds <- DESeqDataSet(se,design = design_formula)
  flog.debug("Successfully prepared DESeq Dataset")
  sizeFactors(dds) <- scran_size_factors
  flog.debug("Successfully computes scran size factors")
  dds$tumor <- relevel(dds$tumor, "non")
  flog.debug('releveled tumor as to contraste against "non"')
  dds<-DESeq(dds, full = mm)
  #dds<-DESeq(dds)
  return(dds)
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

flog.appender(appender.tee(gsub('\\.csv','\\.log',args$output_dir)),name = "ROOT")


flog.info(banner())
flog.info(timestamp(quiet = TRUE))
flog.info(sprintf("Using github version %s", get_git_version()))

#TODO: format better argument print
#flog.info("Program iniated with argumens")
#flog.info(args)

flog.info(sprintf("Will be using %d neighbours for sampling generating %d samples per section | max_ dist is %d", 
                  args$k_members, args$n_samples, args$max_dist))

design_formula <- readLines(file(args$design_file,"r"))
design_formula <- as.formula(design_formula[!grepl('#',design_formula)])
flog.info(paste(c("Using design : ", design_formula ),collapse = " "))
close(file(args$design_file))

matrices <- generate_matrices(count_pth = args$count_dir,
                              feature_pth = args$feature_dir,
                              select_for = args$select_for,
                              feature_name = args$feature_name,
                              k_members = args$k_members,
                              n_samples = args$n_samples,
                              max_dist = args$max_dist,
                              gene_file = args$gene_file,
                              remove_ambigious = args$remove_ambigious,
                              output_dir = args$output_dir)

dds <- DESeq_pipline(matrices$count_matrix,
                     matrices$feature_matrix,
                     design_formula = design_formula)

res <- data.frame(results(dds))
topnres <- subset(res[order(res$padj)[1:min(args$top_n,dim(res)[1])],],select = 'padj')

write.csv(res, file = paste(c(args$output_dir), collapse = "/"))
write.csv(topnres, file = gsub('\\.csv',sprintf('\\.top%d\\.csv',args$top_n),args$output_dir))
flog.info("DGE Analysis Successfully Completed")