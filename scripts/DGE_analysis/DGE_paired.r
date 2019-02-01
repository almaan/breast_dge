#!/usr/bin/Rscript
sh <- suppressPackageStartupMessages

sh(library(SummarizedExperiment))
sh(library(DESeq2))
sh(library(ggplot2))
sh(library(scran))
sh(library(stringr))
sh(library(optparse))
sh(library(plyr))
sh(library(futile.logger))
sh(library(BiocParallel))
sh((library(edgeR)))

# TODO: add uption to use Multicoreparam
register(MulticoreParam(8))

# get script path for library load
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))

# load modules
#source(paste(c(scriptPath,"lib/utils.r"), collapse ="/"))
source(paste(c(scriptPath,"lib/designs.r"), collapse ="/"))
source(paste(c(scriptPath,"lib/parser_paired.r"), collapse ="/"))
source(paste(c(scriptPath,"lib/modOptparse.r"),collapse ="/"))
source(paste(c(scriptPath,"lib/volcano.r"),collapse ="/"))
source(paste(c(scriptPath,"lib/zone_generation.r"),collapse ="/"))

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

subsample_w_preserved_ratio <- function(labels, n_samples) {
  #' Allows for subsamplic of labels where the ratio between
  #' the different labels are kept to the best extent, meaning
  #' that having at least one occurance of all present labels 
  #' will be favored over having a ratio more similar to the 
  #' original.
  #' 
  #' param : labels - a Nx1 vector of labels for N samples
  #' param : n_samples - number of samples to draw 
  #' 
  #' returns : idx_keep - numerical indices of those samples 
  #' to keep. Index refers back to orignal label order.
   
  # get number of indices to draw from each label
  n_elements <- ceiling(table(labels)/sum(table(labels)) * n_samples) #  ceiling used to avoid rounding to zero
  var.names <- names(n_elements)
  
  # Adjust if total number of labels exceed specified
  while(sum(n_elements) > n_samples) {
    n_elements[which.max(n_elements)] <- n_elements[which.max(n_elements)] - 1  
  }
  
  # sample indices to keep from respective label
  idx_keep <- c()
  for (var in var.names) {
    idx_keep <- c(idx_keep, sample(which(labels == var),n_elements[[var]])) #  use which to keep original index reference
  }
  
  return(idx_keep)
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
  #' param : min_per_gene : threshold number of non-zero occurances
  #' 
  #' returns: list with items row_idx and col_idx   
  
  oldrownames <- rownames(mat)
  keep_spots <- rowSums(mat) >= min_per_spot 
  keep_genes <- colSums(mat > 0)/dim(mat)[1] >= min_per_gene
  
  if (remove_ambigious){
    
    n_before <- sum(keep_genes)
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

edgeR_pipeline <- function(count_matrix,
                           feature_matrix,
                           design_formula
                           ) {
  # Should accept count_matrix on format n_genes x n_samples
  dgList <- DGEList(counts= count_matrix, genes = rownames(count_matrix))
  dgList <- calcNormFactors(dgList, method = "TMMwzp")
  designMat <- model.matrix(design_formula, data = feature_matrix, robust = TRUE)
  flog.info("Column Names Design Matrix: ")
  flog.info(colnames(designMat))
  dgList <- estimateDisp(dgList, designMat)
  coef <- dim(designMat)[2]
  
  fit <- glmFit(dgList, designMat)
  lrt <- glmLRT(fit, coef = coef)
  flog.info(paste(c("Used comparision :",as.character(lrt$comparison)),collapse=" "))

  lrt <- topTags(lrt, n = dim(lrt$table)[1],p.value = 0.01)
    
  return(list(fit = fit, res = lrt))
  
  }

DESeq_pipline <- function(count_matrix,
                          feature_matrix,
                          design_formula,
                          contrast,
                          coef) {
  #' Perform DESeq2 DGE analysis with prepared matrices
  
  # construct SE object for compatibility
  se <- SummarizedExperiment(assays = list(counts = count_matrix), colData = feature_matrix)
  flog.debug("Summarized Experiment object generated")
  
  
  # generate a DESeqDataSet object from Summarized Experiment Object
  dds <- DESeqDataSet(se, design = design_formula)
  flog.debug("Successfully prepared DESeq Dataset")
  
  
  decomp_formula <- unlist(strsplit(as.character(design_formula)[2],'\\+'))
  
  if (length(decomp_formula) > 1) {
    reduced_formula <- as.formula(paste(c("~",paste(decomp_formula[-coef],collapse="+")),collapse = ''))
  } else {
    reduced_formula <- as.formula("~1")
  }
  
  # perform DESeq-analysis
  dds<-DESeq(dds, 
             sfType = 'poscounts', # to handle zero counts 
             minmu = 1e-6, #  recommended for SC analysis
             #TODO: look into if LRT is better option
             test = 'LRT', #  recommended for SC analysis, Likelihood Ratio Test
             reduced = reduced_formula,
             useT = TRUE,   #  recommended for SC analysis, use t-distribution as null dist in Wald stat.
             minReplicatesForReplace = Inf, # do not replace outliers
             parallel = TRUE
             )
  
  
  
  if (!is.null(contrast)) {
    results_dds <- results(dds, contrast = contrast, minmu = 10e-6)
    # TODO: evaluate if lfcShrinkage is necessary
    #try(results_dds <-  lfcShrink(dds,res = results_dds,contrast = contrast))  
  } else {
    reults_dds <- results(dds)
  }
  
  flog.info("Head of results object :")
  flog.info(head(results_dds))  
  return(list(fit = dds, res = results_dds))
}


generate_matrices <- function(path_feat,
                              path_cnt,
                              min_per_spot,
                              min_per_gene,
                              ss_number,
                              ss_feature,
                              keep_cols,
                              zone_distance,
                              filter_tumors = FALSE,
                              remove_ambigious = FALSE) {
  
  #' Returns a list of a concatenated feature matrix and count matrix
  #' Data from single replicates and patients are joined into larger data frame
  
  clist <- list()
  feature_matrix <- data.frame(stringsAsFactors = TRUE)
  
  patient_names <- c("")
  gene_names <- c()
  
  num_spots <- 0
  num_patient <- 0  
    
  for (num in 1:length(path_cnt)) {
    
    flog.info(paste(c("Reading feature file :", path_feat[num]),collapse = ' '))
    # load feature file
    fmat <- read.csv(path_feat[num], sep = "\t",
                     header = TRUE, row.names = 1, stringsAsFactors = TRUE)
    
    
    # if only tumor-annotated spots are to be analyzed
    if (filter_tumors) {
      fmat <- fmat[fmat['tumor'] == 'tumor',] #remove non-tumor spots
    }
    
    if (ss_number > 0 & ss_number <= dim(fmat)[1]) {
      fmat <- fmat[subsample_w_preserved_ratio(labels = fmat[[ss_feature]], n_samples = ss_number),]
    }
    
    if (!(is.na(zone_distance))) {
      zones <- make_zones(crd = cbind(fmat$xcoord,fmat$ycoord),
                          labels = fmat[["tumor"]],
                          zone_method = "three_levels",
                          ulim = zone_distance,
                          foci_label = "tumor")
      
      fmat$zones <- zones 
      print(unique(fmat$zone))
    }
    
    flog.info(paste(c("Reading count file :", path_cnt[num]),collapse = ' '))
    # load count-matrix
    cmat <- read.csv(path_cnt[num], sep = "\t",
                     header = TRUE, row.names = 1)
    
    
    if(length(rownames(fmat)) != length(rownames(cmat)) & ss_number == 0) {
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
    
    #add new gene names
    gene_names <- union(gene_names, unlist(colnames(cmat))) # union to avoid duplicate
    #add total number of spots
    num_spots <- num_spots + length(inter)
    
    #add count matrix to list
    clist <- append(clist, list(cmat))
    feature_matrix <- rbind(feature_matrix,as.data.frame(fmat))
    
    
  } 
  
  # give unique replicate names
  feature_matrix['replicate'] <- paste(feature_matrix$patient,feature_matrix$replicate,sep=".")
  
  # Assemble full count matrix --------
  # create blank empty matrix of dimension n_spots x n_genes
  count_matrix <- data.frame(matrix(0, num_spots, length(gene_names)))
  # set proper column names
  colnames(count_matrix) <- gene_names
  
  # fill matrix from stored count matrices
  row.start <- 1
  for (i in 1:length(path_cnt)) {
    add.n.rows <- dim(clist[[1]])[1] - 1 #  subtraction to not overcount
    count_matrix[row.start:(row.start+ add.n.rows),colnames(clist[[1]])] <- clist[[1]] # fill interval in count matrix
    clist[1] <- NULL #  dequeue matrix, to free up memory
    row.start <- row.start + add.n.rows #  move to blank row in count_matrix  
  }
  
  # remove potential NA and 
  if(any(is.na(count_matrix))) {
    flog.warn("Count Matrix contains NA elements. These are set to zero")
    count_matrix[is.na(count_matrix)] <- 0
  }
  # get indices for "cleaned" matrix
  indices <- get_clean_indices(count_matrix,
                             min_per_spot = min_per_spot,
                             min_per_gene = min_per_gene,  
                             remove_ambigious = remove_ambigious)
  
  count_matrix <- count_matrix[indices$row_idx,indices$col_idx]
  
  # sync feature and count matrix
  feature_matrix <- feature_matrix[indices$row_idx,]
  # put turn all columns into factors
  feature_matrix[] <- lapply(feature_matrix,factor)
  
  if (!is.na(zone_distance)) {
    zones <- as.numeric(levels(feature_matrix$zones))[feature_matrix$zones]
    uni_zones <- unique(zones)
    print(uni_zones)
    new_levels <- c(uni_zones[-which(uni_zones == 1)], 1)
    feature_matrix$zones <- factor(zones, levels = new_levels)
    
    if (any(feature_matrix$zones == 2)) {
      feature_matrix$zones <- relevel(feature_matrix$zones,ref = "2")
      flog.info("Using Peripheral region as base-level")
    } else {
      feature_matrix$zones <- relevel(feature_matrix$zones,ref = "0")
      flog.info("Using Tumor region as base-level")
    }
  }
  
  # relevel tumor factors as to contrast tumor to non-tumor
  if (any('tumor' == colnames(feature_matrix))) {
    try(feature_matrix$tumor <- relevel(feature_matrix$tumor, "non"))
    flog.info('releveled tumor as to contrast against "non"')
  }
  
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

flog.info(paste(c("Command line arguments :", paste(args,collapse=" ")), collapse = "\n"))

# get filenames
if (length(args$count_file)!=length(args$feature_file)) {
  flog.warn(paste(c("The number of count files does not match",
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


# information for logging
flog.info(paste(c("count matrix files(s) used : ",count_pth ),collapse = "\n"))
flog.info(paste(c("feature matrix files(s) used : ",feature_pth ),collapse = "\n"))
flog.info(paste(c("Removing genes with lower than", args$min_per_gene,"average occurance"),collapse=" "))
flog.info(paste(c("Removing spots with less than", args$min_per_spot,"total transcripts"),collapse=" "))

if (args$subsample_number > 0) {
  flog.info(sprintf("Will subsample data taking %d spots from each section w.r.t to feature %s",
                    args$subsample_number, args$subsample_feature))  
}


# set design to use
if (args$design < 1 & length(args$custom_design) < 1) {
      # if flag used but no value given (default -1)
      # print pre-constructed formulas and exit
      pre_design_formulas(args$design)
      exit()
      
  } else if (length(args$custom_design) >= 1) {
    # if a custom designed is provided 
    design_formula <- as.formula(args$custom_design)
    contrast <- args$contrast
  } else {
    pre_expression <- pre_design_formulas(args$design)
    design_formula <- as.formula(pre_expression$design)
    contrast <- pre_expression$contrast
  }

if(is.na(args$zero_coef)) {
  coef <- as.numeric(length(strsplit(as.character(design_formula)[2],'\\+')))
} else {
  coef <- as.numeric(args$zero_coef)
}

flog.info(paste(c("design formula :", design_formula),collapse =" "))
if (!is.null(contrast)) {
  flog.info(paste(c("contrast triple :", paste(contrast,collapse =" ")),collapse =" "))
}
# Main Body ---------------------------------

keep_cols <- all.vars(design_formula)

# generate concatenated (full) feature matrix and count matrix
matrices <- generate_matrices(path_feat = args$feature_file, 
                              path_cnt = args$count_file,
                              min_per_spot = args$min_per_spot,
                              min_per_gene = args$min_per_gene,
                              keep_cols = keep_cols,
                              ss_number = args$subsample_number,
                              ss_feature = args$subsample_feature,
                              zone_distance = args$zone_distance,
                              remove_ambigious = args$remove_ambigious,
                              filter_tumors = args$filter_tumors 
                              )


flog.info("Successfully generated joint feature and count matrices")
flog.info(sprintf("Will be using %d spots and %d genes in analysis",
                  dim(matrices$count_matrix)[1],dim(matrices$count_matrix)[2]))
flog.info("Initiate Differential Gene Expression analysis")

# initate DESeq2 with generated matrices and provided design formula
if (args$method == 'deseq2') {
  flog.info("Using method DESeq2")
  deseqobj <- DESeq_pipline(t(matrices$count_matrix), 
                               matrices$feature_matrix,
                               design_formula,
                               contrast,
                               coef)
  results_dge <- deseqobj$res
  fit_dge <- deseqobj$fit
  pval_col <- "padj"
  lfc <- "log2FoldChange"
  
} else if (args$method == 'edgeR') {
  flog.info("Using method edgeR")
  edgeRobj   <- edgeR_pipeline(t(matrices$count_matrix),
                                 matrices$feature_matrix,
                                 design_formula
                                 )
  results_dge <- edgeRobj$res
  fit_dge <- edgeRobj$fit
  pval_col <- "FDR"
  lfc <- "logFC"
} else {
  flog.error("Method not supported. Exiting")
  exit()
}


flog.info(paste(c("Successfully Completed DGE analysis. Now Saving results to files")))

# Save Results ---------------------------------

res <- data.frame(results_dge)

if (args$rsave) {
  try(save(fit_dge, file = paste(c(args$output_dir,paste(c('r.res',session_tag,'r'),
                                                     collapse='.')),collapse = "/")))
  flog.info("Saved R object")
}

# if fancy output is to be produced as well
if (args$fancy) {
  sh(library(org.Hs.eg.db)) 
  sh(library(AnnotationDbi))
  # select for genes with padj equal or less to specifed value
  topnres <- res[res[pval_col] <= args$pval & !(is.na(res[pval_col])),]
  
  # only do fancy if stat. signigicant genes have been detected 
  if (dim(topnres)[1] > 1) {
    # include gene symbols in output
    topnres$symbol <-mapIds(org.Hs.eg.db, keys=rownames(topnres),
                            column="SYMBOL", keytype="ENSEMBL", multiVals="first")
    
    try(write.csv(topnres, file = paste(c(args$output_dir,paste(c(session_tag,'fancy.tsv'),collapse='.')),
                                                                                          collapse = "/")))
    flog.info("Save Fancy results")  
    }
  
}

if (args$volcano) {
  filtered_idx <- which(!is.na(res[pval_col]))
  volcano.plot <- volcano_plot(res[filtered_idx,lfc],res[filtered_idx,pval_col])
  try(ggsave(filename = paste(c(args$output_dir,paste(c(session_tag,'volcano.png'),collapse='.')),collapse = "/"), volcano.plot))
  flog.info("Saved volcano plot based on results")
}

