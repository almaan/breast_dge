
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
source(paste(c(getwd(),"lib/utils.r"), collapse ="/"))
source(paste(c(getwd(),"lib/designs.r"), collapse ="/"))


##----------Functions-----------##
sym_diff <- function(a,b) setdiff(union(a,b), intersect(a,b))

exit <- function() {
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}

banner <- function(){
  txt <- paste(c("",
                 " __..___.     .__ .__ .___  .__.      .",          
                 "(__   |   ___ |  \\[ __[__   [__]._  _.|  . __* __",
                 ".__)  |       |__/[_./[___  |  |[ )(_]|\\_|_) |_)",
                 "                                       ._|      "), collapse ="\n")
  return(txt)
}


clean_matrix <- function(mat, remove_ambigious = FALSE,
                         min_sample = 100 , min_gene = 0.1){
  #matrix should be formatted n_samples x n_genes
  oldrownames <- rownames(mat)
  keep_samples <- rowSums(mat) >= min_sample 
  keep_genes <- colMeans(mat) >= min_gene
  mat <- mat[keep_samples,keep_genes]
  
  if (remove_ambigious){
    n_before <- dim(mat)[2]
    mat <-mat[,!grepl('.*ambig.*',colnames(mat))]
    n_after <- dim(mat)[2]
    flog.info(sprintf("Removed %d ambigious genes out of %d total", (n_before-n_after), n_before))
  }
  flog.info(paste(c("Removed samples with index :", sym_diff(rownames(mat),oldrownames)),collapse = " "))
  
  return(mat)
}



DESeq_pipline <- function(count_matrix, feature_matrix, design_formula) {
  #Perform DESeq2 DGE analysis with prepared matrices
  #use scran as to avoid error with zero counts
  scran_size_factors <- computeSumFactors(count_matrix, positive = TRUE)
  flog.debug("size factors extracted")
  #construct SE object for compatibility
  se <- SummarizedExperiment(assays = list(counts = count_matrix), colData = feature_matrix)
  flog.debug("Summarized Experiment object generated")
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
  #TODO: implement relevel if tumor is in columns
  #dds$tumor <- relevel(dds$tumor, "non")
  #flog.debug('releveled tumor as to contraste against "non"')
  dds<-DESeq(dds, full = mm)
  #dds<-DESeq(dds)
  return(dds)
}

get_filenames <- function(path, sample_list) {
  files <- list.files(path,pattern = '*.tsv')
  samples <- c()
  for (sample in sample_list) {
    samples <- c(samples,files[sapply(files, function(x) grepl(sample,x))])
  }
  return(samples)
}

pseudo_replicate <- function(fm) {
  ids <- unique(fm[['id']])
  new <- fm[['replicate']]
  for (sid in ids) {
    reps <- unique(fm[fm['id']== sid,'replicate'])
    psid = 1
    for (r in reps) {
      new[(fm['id'] == sid) & (fm['replicate'] == r)] = psid
      psid = psid + 1
    }
  }
  return(as.factor(new))
}

generate_matrices <- function(path_feat, path_cnt,samples_feat, samples_cnt,
                              filter_tumors = FALSE, remove_ambigious = FALSE) {
  
  if (length(samples_feat) != length(samples_cnt)) {
    warning("Inconsistent number of files provided")
  }
  
  count_matrix <- data.frame()
  feature_matrix <- data.frame()
  #assuming samples are ordered by row
  #assuming genes are ordered by column
  for (num in 1:length(samples_feat)) {
    flog.info(paste(c("Working with file :", samples_feat[num]),collapse = ' '))
    fmat <- read.csv(paste(c(path_feat,samples_feat[num]),collapse="/"), sep = "\t",
                     header = TRUE, row.names = 1, stringsAsFactors = TRUE)
    
    fmat <- fmat[,c('tumor','id','replicate')]
    
    #TODO:something problematic here!
    #if only tumors are to be analyzed, reduces mem.usage
    if (filter_tumors) {
      #fmat <- fmat[fmat['tumor'] == 'tumor',]
    }
    
    cmat <- read.csv(paste(c(path_cnt,samples_cnt[num]),collapse="/"), sep = "\t",
                     header = TRUE, row.names = 1)
    
    inter <- intersect(rownames(fmat),rownames(cmat))
    cmat <-cmat[inter,]
    fmat <- fmat[inter,]
    
    count_matrix <- bind_rows(count_matrix,cmat)
    feature_matrix <- bind_rows(feature_matrix,fmat)
    
  } 
  
  feature_matrix['pseudo.replicate'] <- pseudo_replicate(feature_matrix)
  count_matrix[is.na(count_matrix)] <- 0
  count_matrix <- clean_matrix(count_matrix, remove_ambigious = remove_ambigious)
  feature_matrix <- feature_matrix[intersect(rownames(feature_matrix),rownames(count_matrix)),]
 

  return(list(count_matrix = count_matrix, feature_matrix = feature_matrix))
}

##-------------PARSER-----------##

parser <- OptionParser()


#currently only supports single file directory
# TODO: implement multiple directory
#idea here https://stackoverflow.com/questions/13790610/passing-multiple-arguments-via-command-line-in-r

print(banner())

parser <- add_option(parser,
                     c("-c","--count_dir"),
                     default = c(),
                     help = paste(c("directory of count-matrices.",
                                    "matrix name should be on form",
                                    '"count_data-ID_Replicate.tsv"'),
                                  collapse = " "))

parser <- add_option(parser,
                     c("-f","--feature_dir"),
                     default = c(),
                     help = paste(c("directory of feature-files.",
                                    "file name should be on form",
                                    '"count_data-ID_Replicate.tsv"'),
                                  collapse = " "))

parser <- add_option(parser,
                     c("-o", "--output_dir"),
                     default = getwd(),
                     type = "character",
                     help = paste(c("directory to save output into.",
                                    "if none specified cwd will be used."),
                                  collapse = " "))

parser <- add_option(parser,
                     c('-s','--samples'),
                     default = c(),
                     help = paste(c("samples to be used in analysis",
                                    "if only patient ids and not replicates",
                                    "are provided all replicates of that patient",
                                    "will be used"),collapse = " ")
)

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

args <- parse_args(parser)
                     
if (FALSE){                     
    if (args$design < 1 & length(args$custom_design) < 1) {
        pre_design_formulas(args$design)
        exit()
    }else if (length(args$custom_design) >= 1) {
    design_formula <- as.formula(args$custom_design)
    } else {
    design_formula <- pre_design_formulas(args$design)
    }
}

if (grepl(',',args$samples)) {
  sample_ids <- gsub(' ','',args$samples)
  sample_ids <-unlist(strsplit(sample_ids,','))
} else {
  sample_ids <- args$samples                    
}
##-----------Main--------------#


#sample_id1 <- '24105'
#sample_id2 <- '23377'
#sample_ids <- c(sample_id1, sample_id2)
#path_feat <- '/home/alma/ST-2018/CNNp/DGE/res/results_extraction6'
#path_cnt <-'/home/alma/ST-2018/CNNp/data/count_data'

#vital that files are in same order upon load
samples_feat <- get_filenames(args$feature_dir,sample_ids) 
samples_cnt <- get_filenames(args$count_dir,sample_ids)

print(args$samples)
print(samples_feat)
flog.info("Successfully loaded samples. 
          Assuming naming is consistent between counts and features.")

matrices <- generate_matrices(path_feat = args$feature_dir, path_cnt = args$count_dir,
                              samples_feat = samples_feat, samples_cnt = samples_cnt,
                              remove_ambigious = args$remove_ambigious,
                              filter_tumors = args$filter_tumors)

feature_matrix = matrices$feature_matrix
count_matrix = matrices$count_matrix
flog.info("Successfully generated joint feature and count matrices")
#design_formula <- as.formula('~ id:pseudo.replicate + pseudo.replicate + id ')
design_formula <- as.formula('~ pseudo.replicate + pseudo.replicate:tumor + tumor')
flog.info("Initiate DESeq2 Pipline")
start_time <- Sys.time()
dds <- DESeq_pipline(t(count_matrix), feature_matrix,design_formula)
end_time <- Sys.time()

flog.info(paste(c("Successfully Completed DESeq2 analysis in >>",as.numeric(end_time - start_time),"seconds"),
                collapse = " "))

res <- data.frame(results(dds))
topnres <- subset(res[order(res$padj)[1:min(100,dim(res)[1])],],select = 'padj')
 
write.csv(res, file = paste(c(args$output_dir,paste(c(args$samples,'tsv'),collapse='.')), collapse = "/"))
write.csv(topnres, file = paste(c(args$output_dir,paste(c(args$samples,'top100.tsv'),collapse='.')), collapse = "/"))