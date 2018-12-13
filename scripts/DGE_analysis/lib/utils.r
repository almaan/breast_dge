library(futile.logger)
get_git_version <- function(){
  #' prints the git version of the used source file
  sys.args <- commandArgs(trailingOnly = FALSE) #get input arguments inclusive path
  has_filename <- sys.args[grepl('--file=',sys.args, perl = TRUE)] #get the argument containing filename
  
  filepath <- gsub('--file=',"", #replace anchor --file with non to only get path
                   has_filename,
                   ignore.case = FALSE)
  
  execpth <- dirname(filepath) #get directory of sourcefile
  
  cmd <- sprintf("cd %s; git describe --tags; cd %s", execpth, getwd()) #command to execute
  version <- system(command = cmd, intern = TRUE) #intern = TRUE to catch return from system call
  
  return(version)
}

save_tmp <- function(count_matrix, feature_matrix, output_dir){
  #' saving temporary file in-case of crash as well as debug purpose
  flog.info('save temporary feature matrix')
  write.table(feature_matrix,file = paste(c(dirname(output_dir),'tmp.feature.file.tsv'),collapse="/"), sep = "\t")
  flog.info("save temporary count matrix")
  write.table(count_matrix,file = paste(c(dirname(output_dir),'tmp.count.file.tsv'),collapse="/"), sep = "\t")
}

load_matrix <- function(main_pth, file_name, nrows = -1){
  #' Load count matrix expects format samples x genes
  mat <- read.csv(paste(c(main_pth,file_name),collapse="/"), sep = "\t", nrows = nrows, header = TRUE, row.names = 1)
  return(mat)
}
load_features <- function(main_pth,file_name) {
  #' load feature file complementary to count matrix
  feat <- read.csv(paste(c(main_pth,file_name),collapse="/"), sep = "\t",
                   header = TRUE, row.names = 1, stringsAsFactors = TRUE)
  return(feat)
}

get_id <- function(path) {
  #'get sample id. Expects a 4-5 character long patient id with the
  #'two character replicate id appended to this separated by an underscore
  return(str_extract(path,pattern = "\\d{4,5}_\\w{2}"))
}

clean_matrix <- function(mat, remove_ambigious = TRUE, min_sample = 0.1 , min_gene = 1){
  #matrix should be formatted n_samples x n_genes
  keep_samples <- rowMeans(mat > 0) >= min_sample 
  keep_genes <- colMeans(mat > 0) >= min_gene
  mat <- mat[keep_samples,keep_genes]
  if (remove_ambigious){
    n_before <- dim(mat)[2]
    mat <-mat[,!grepl('.*ambig.*',colnames(mat))]
    n_after <- dim(mat)[2]
    flog.info(sprintf("Removed %d ambigious genes out of %d total", (n_before-n_after), n_before))
  }
  return(mat)
}


banner <- function(x){
  #' Completely unnecessary banner
  txt <- paste(c("",
                 " __..___.     .__ .__ .___  .__.      .",          
                 "(__   |   ___ |  \\[ __[__   [__]._  _.|  . __* __",
                 ".__)  |       |__/[_./[___  |  |[ )(_]|\\_|_) |_)",
                 "                                       ._|      "), collapse ="\n")
  return(txt)
}
