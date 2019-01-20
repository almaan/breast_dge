#!/usr/bin/Rscript
sh <- suppressPackageStartupMessages

sh(library(SummarizedExperiment))
sh(library(DESeq2))
sh(library(ggplot2))
sh(library(scran))
sh(library(stringr))
sh(library(optparse))
sh(library(dplyr))
sh(library(futile.logger))

# get script path for library load
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))

# load modules
#source(paste(c(scriptPath,"lib/utils.r"), collapse ="/"))
source(paste(c(scriptPath,"lib/designs.r"), collapse ="/"))
source(paste(c(scriptPath,"lib/parser_paired.r"), collapse ="/"))
source(paste(c(scriptPath,"lib/modOptparse.r"),collapse ="/"))

# Function Space ---------------------------------


get_session_tag <- function() {
  timestamp <- gsub(pattern = ':| ', replacement = "-", x = as.character(Sys.time())) 
  return(paste(c('DGE_analysis.',timestamp,'.',round(runif(1)*10e4,0)),collapse = ""))
}

sym_diff <- function(a,b) { 
  # returns elements of two sets a and b which are 
  # present in one set but not both
  return(setdiff(union(a,b), intersect(a,b)))
  }

exit <- function() {
  # exit function
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}

banner <- function(){
  # completely unecessary banner.
  txt <- paste(c("",
                 " __..___.     .__ .__ .___  .__.      .",          
                 "(__   |   ___ |  \\[ __[__   [__]._  _.|  . __* __",
                 ".__)  |       |__/[_./[___  |  |[ )(_]|\\_|_) |_)",
                 "                                       ._|      "), collapse ="\n")
  return(txt)
}

get_clean_indices <- function(mat, remove_ambigious,
                         min_per_spot, min_per_gene){
  
  #' Returns indices of rows and columns to keep in
  #' provided count matrix after cleaning specified
  #' by user. Matrix should be formatted n_samples x n_genes
  #' 
  #' param : mat - count_matrix
  #' param : remove_ambigious - remove all ambigiously mapped genes
  #' param : min_per_spot : threshold for total spot transcript count
  #' param : min_per_gene : threshold for avearage gene transcript count
  #' 
  #' returns: list with items row_idx and col_idx   
  
  oldrownames <- rownames(mat)
  keep_spots <- rowSums(mat) >= min_per_spot 
  keep_genes <- colMeans(mat) >= min_per_gene
  
  if (remove_ambigious){
    
    n_before <- length(keep_genes)
    keep_genes[grepl('.*ambig.*',colnames(mat))] <- FALSE 
    n_after <- sum(keep_genes)
   
     if(n_after != n_before) {
      flog.info(sprintf("Removed %d ambigious genes out of %d total", (n_before-n_after), n_before))
      } else {
      flog.info("No ambigious reads detected. All genes are kept.")
    }
  }
  
  return(list(row_idx = keep_spots, col_idx = keep_genes))
}


DESeq_pipline <- function(count_matrix, feature_matrix, design_formula) {
  #' Perform DESeq2 DGE analysis with prepared matrices
  #' use scran as to avoid error with zero counts
  
  scran_size_factors <- computeSumFactors(count_matrix, positive = TRUE)
  flog.info("Successfully computed scran size factors")
  # construct SE object for compatibility
  se <- SummarizedExperiment(assays = list(counts = count_matrix), colData = feature_matrix)
  flog.debug("Summarized Experiment object generated")
  
  # construct model matrix using specified design formulae
  mm <- model.matrix(design_formula, colData(se))
  # identify elements where no interaction (all zero) are present
  all.zero <- apply(mm, 2, function(x) all(x==0))
  # get indices for zero columns
  zidx <- which(all.zero)
  
  if (length(zidx) > 0) {
    # if any zero columns, remove these for linear independency
    # recommended procedure by Love
    mm <- mm[,-zidx]
  }
  
  # generate a DESeqDataSet object from Summarized Experiment Object
  dds <- DESeqDataSet(se,design = design_formula)
  flog.debug("Successfully prepared DESeq Dataset")
  sizeFactors(dds) <- scran_size_factors
  
  # relevel tumor factors as to always contrast
  # tumor to non-tumor
  if ('tumor' %in% colnames(feature_matrix)) {
    dds$tumor <- relevel(dds$tumor, "non")
  flog.info('releveled tumor as to contraste against "non"')
  }
  # perform DESeq-analysis
  dds<-DESeq(dds, full = mm)
  return(dds)
}


pseudo_replicate <- function(fm) {
  #' Returns vector of new replivate names
  #' DESeq2 requires that columns are linearly independent
  #' in design matrix. Each replicate identifier will be 
  #' assignes a numeric value from 1:n_replicates. 
  #' n_replicates is number of replicates of one patient.
  
  ids <- unique(fm[['id']]) #all unique patient id's
  new <- matrix("",dim(fm)[1]) #copy of current replicates

    # iterate over all patients in dataframe
  for (sid in ids) { #  iterate over patient
    reps <- unique(fm[fm['id']== sid,'replicate'])
    psid = 1
    for (r in reps) { #  iterate over replicate

      idx <- (fm['id'] == sid) & (fm['replicate'] == r)
      new[idx] = paste(c(sid,as.character(psid)),collapse=".")
      psid = psid + 1
    
    }
  }
  # return as factors
  return(as.factor(new))
}

generate_matrices <- function(path_feat,
                              path_cnt,
                              min_per_spot,
                              min_per_gene,
                              filter_tumors = FALSE,
                              remove_ambigious = FALSE) {
  
  #' Returns a list of a concatenated feature matrix and count matrix
  #' Data from single replicates and patients are joined into larger data frame
  
  
  clist <- list()
  feature_matrix <- data.frame(stringsAsFactors = TRUE)
  
  gene_names <- c()
  num_spots <- 0
  
  for (num in 1:length(path_cnt)) {
    
    flog.info(paste(c("Reading feature file :", path_feat[num]),collapse = ' '))
    # load feature file
    fmat <- read.csv(path_feat[num], sep = "\t",
                     header = TRUE, row.names = 1, stringsAsFactors = TRUE)
    
    # extract only relevant columns
    #TODO: make these non static
    fmat <- fmat[,c('tumor','id','replicate')]
    
    # if only tumor-annotated spots are to be analyzed
    if (filter_tumors) {
      fmat <- fmat[fmat['tumor'] == 'tumor',] #remove non-tumor spots
    }
    
    flog.info(paste(c("Reading count file :", path_cnt[num]),collapse = ' '))
    # load count-matrix
    cmat <- read.csv(path_cnt[num], sep = "\t",
                     header = TRUE, row.names = 1)
    
    
    if(length(rownames(fmat)) != length(rownames(cmat))) {
      flog.warn(paste(c("feature matrix and column matrix ",
                           "conatins different number of entries.",
                           "Analysis will continue but unexpected results",
                           "may arise as a consequence of this."
      ),collapse = " "))
    }
    
    # find spots present in count and feature file. Set to same order
    inter <- intersect(rownames(fmat),rownames(cmat))
    
    cmat <-cmat[inter,]
    fmat <- fmat[inter,]
    
    gene_names <- union(gene_names, unlist(colnames(cmat)))
    num_spots <- num_spots + length(inter)
    
    clist <- append(clist, list(cmat))
    feature_matrix <- rbind(feature_matrix,as.data.frame(fmat))
    
    
  } 
  
  count_matrix <- data.frame(matrix(0, num_spots, length(gene_names)))
  colnames(count_matrix) <- gene_names
  pos <- 1
  
  for (i in 1:length(path_cnt)) {
    nsp <- dim(clist[[1]])[1] - 1
    count_matrix[pos:(pos+ nsp),colnames(clist[[1]])] <- clist[[1]]
    clist[1] <- NULL
    pos <- pos + nsp 
  
  }
  
  # Rename replicates to obtain linear independent columns in downstream DESeq2 analysis
  feature_matrix['pseudo.replicate'] <- pseudo_replicate(feature_matrix)
  # clean count matrix
  count_matrix[is.na(count_matrix)] <- 0
  
  indices <- get_clean_indices(count_matrix,
                             min_per_spot = min_per_spot,
                             min_per_gene = min_per_gene,  
                             remove_ambigious = remove_ambigious)
  
  count_matrix <- count_matrix[indices$row_idx,indices$col_idx]
  # sync feature and count matrix after cleanup
  feature_matrix <- feature_matrix[indices$row_idx,]
  rownames(feature_matrix) <- rownames(count_matrix)

  return(list(count_matrix = count_matrix, feature_matrix = feature_matrix))
}

# Parsing ---------------------------------

session_tag <- get_session_tag()

parser <- OptionParser()
parser <- make_parser(parser,tag = session_tag) #  add options, delegated to external function
# parse arguments; modified to support multiargs
args <- splitMultipleArgs(parse_args(parser, args = allowMultipleArgs())) 
# set logger to write both to stdout and file
flog.appender(appender.tee(paste(c(args$output_dir,args$logname),
                                 collapse="/")),name = "ROOT")
flog.info(banner()) #aestethics

# get filenames
if (length(args$count_file)!=length(args$feature_file)) {
  flog.warning(paste(c("The number of count files does not match",
                        "the number of feature files. Exiting"), collapse =" "))
  exit()
}

# get full paths to count and feature files
# TODO: make more elegant
if (ifelse(length(args$count_file) == 1,
           (dir.exists(args$count_file) & 
            dir.exists(args$feature_file)),FALSE)) {
  
  count_pth <- list.files(args$count_file)
  feature_pth <- list.files(args$feature_file)
} else {
  count_pth <- args$count_file
  feature_pth <- args$feature_file
}


flog.info(paste(c("count matrix files(s) used : ",count_pth ),collapse = "\n"))
flog.info(paste(c("feature matrix files(s) used : ",count_pth ),collapse = "\n"))

# set design to use
if (args$design < 1 & length(args$custom_design) < 1) {
      # if flag used but no value given (default -1)
      # print pre-constructed formulas and exit
      pre_design_formulas(args$design)
      exit()
      
  } else if (length(args$custom_design) >= 1) {
    # if a custom designed is provided 
    design_formula <- as.formula(args$custom_design)
  } else {
    design_formula <- as.formula(pre_design_formulas(args$design))
  }

flog.info(paste(c("design formula :", design_formula),collapse =" "))

# Main Body ---------------------------------

# generate concatenated (full) feature matrix and count matrix
matrices <- generate_matrices(path_feat = args$feature_file, 
                              path_cnt = args$count_file,
                              min_per_spot = args$min_per_spot,
                              min_per_gene = args$min_per_gene, 
                              remove_ambigious = args$remove_ambigious, 
                              filter_tumors = args$filter_tumors 
                              )

flog.info("Successfully generated joint feature and count matrices")
flog.info("Initiate DESeq2")

# initate DESeq2 with generated matrices and provided design formula
dds <- DESeq_pipline(t(matrices$count_matrix), matrices$feature_matrix,design_formula)
results_dds <- results(dds)

remove(matrices) #  release memory
flog.info(paste(c("Successfully Completed DESeq2 analysis. Saving Output to files")))

# Save Results ---------------------------------

res <- data.frame(results_dds)

try(write.csv(res, file = paste(c(args$output_dir,paste(c(session_tag,'tsv'),collapse='.')),
                                                              collapse = "/"), sep = "\t"))


if (args$paranoid) {
  try(save(results_dds, file = paste(c(args$output_dir,paste(c('res',session_tag,'r'),
                                                      collapse='.')),collapse = "/")))
  try(save(dds, file = paste(c(args$output_dir,paste(c('dds',session_tag,'r'),
                                              collapse='.')), collapse = "/")))
}

# if fancy output is to be produced as well
if (args$fancy) {
  sh(library(org.Hs.eg.db)) 
  sh(library(AnnotationDbi))
  # select for genes with padj equal or less to specifed value
  topnres <- res[res$padj <= args$pval & !(is.na(res$padj)),]
  # include gene symbols in output
  topnres$symbol <-mapIds(org.Hs.eg.db, keys=rownames(topnres),
                            column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  
  try(write.csv(topnres, file = paste(c(args$output_dir,paste(c(session_tag,'fancy.tsv'),collapse='.')),
                                      collapse = "/"), sep = "\t"))

}

# TODO : Implement volcano plot function print as visualization of results