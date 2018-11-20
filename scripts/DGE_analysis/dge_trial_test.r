suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("scran"))

suppressPackageStartupMessages(library("docstring"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("futile.logger"))

system("git describe --tags")

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


prepare_counts <-function(main_pth,single_sample, lower_bound = 0){
  #'Load and format the count-matrix of a .tsv file
  #'Input file must have genes as columns and spots as rows
  #'Will remove all spots with a total count less than lower_bound
  #'If a specific gene-list is provided by the user, only these
  #'genes will be read, else all union of genes found in all samples will
  #'be read.
  #'
  #'Input : 
  #'   - main_pth : directory where .tsv count matrix files are located
  #'   - single_sample : file-name of count matrix. Format spots x genes
  #'   - lower_bound : lower bound of total count of a spot. Default 0
  #'   
  #'Output : 
  #'   - tcnt : count matrix dataframe with spots as rows, genes as columns 
  
  full_pth <- paste(main_pth,single_sample,sep = "/")
  sample_name <- extract_names(single_sample)
  
  flog.info(paste(c('Including', sample_name, "into count matrix"), collapse = " "))
  
  if (!typeof(SELECTED_GENES) == "logical"){
    
    header <- unlist(strsplit(readLines(file(full_pth), n = 1),'\t'))
    close(file(full_pth))
    
    gene_names <- which(header %in% SELECTED_GENES)
    colclass <- replicate(length(header), NULL)
    colclass[gene_names] <- NA
    colclass[1] <- NA
    
    remove(header)
  } 
  
  tcnt <- as.data.frame(read.csv(full_pth,
                                   sep = '\t',
                                   header = TRUE,
                                   colClasses = colclass,
                                   row.names = 1,
                                   ))
  
  rownames(tcnt) <- sapply(rownames(tcnt), function(x) paste(sample_name,x, sep = '_'))
  
  prior_n_spots <- dim(tcnt)[1]
  tcnt <- tcnt[rowSums(tcnt) > lower_bound,]
  post_n_spots <- dim(tcnt)[1]
  
  flog.info(paste(c("A total of", 
                    prior_n_spots - post_n_spots,
                    "/",
                    prior_n_spots,
                    "spots were removed",
                    "due to low counts"),
                  collapse = " ")
            )
  
  return(tcnt)
}

prepare_coldata <- function(feature_file_dir,single_sample) {
  #' Prepare coldata dataframe for DESeq2
  #' 
  #'  Input : 
  #'    - feature_file_dir : directory of feature_file_dir specified for sample
  #'    - single_sample : file-name of sample to be loaded
  #'     
  #'  Output : 
  #'     - cdata : colData data-frame to be used in DESeq2
  
  sample_name <- extract_names(single_sample)
  flog.info(paste(c('Including', sample_name, "into feature matrix"), collapse = " "))
  
  patient_id <- unlist(strsplit(sample_name,'_'))[1]
  replicate <- unlist(strsplit(sample_name,'_'))[2]
  
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
      
      cnt <- prepare_counts(COUNT_DATA_DIR,
                            count_file_names[num])
      
      coldata <- prepare_coldata(FEATURE_FILE_DIR,
                                 count_file_names[num])
      
    } else {
      cnt <- merge_matrices(cnt,
                            prepare_counts(COUNT_DATA_DIR,
                            count_file_names[num]))
      
      coldata <- merge_matrices(coldata, 
                                prepare_coldata(FEATURE_FILE_DIR,
                                count_file_names[num]))
    }
  } 
  
  
  flog.info("Initate Size Factor Estimation using scran")
  
  scran_size_factors <- computeSumFactors(as.matrix(t(cnt)), positive = TRUE)
  
  flog.info("Size Factors Computed Using scran")
  
  # TODO : Look into how to remedy this ugly shit
  
  keep_spots <- (scran_size_factors > 0.0)
  scran_size_factors <- scran_size_factors[keep_spots]
  cnt <- cnt[keep_spots,]
  
  coldata <- coldata[rownames(cnt),]
  
  if (sum(rownames(coldata) != rownames(cnt)) == 0) {
    flog.info("Features are coherent with Count Matrix")
  } else {
    flog.error("Features are not coherent with Count Matrix")
    quit("no",0)
  }
  
  
  dds <- DESeqDataSetFromMatrix(countData = t(cnt),
                                colData = coldata,
                                design = design)
  
  flog.info("DESeq2 Matrix generated")
  
  sizeFactors(dds) <- scran_size_factors
  
  flog.info("Transfered size factors from scran to DESeq2")  
  
  dds <- DESeq(dds, 
              parallel = TRUE,
              BPPARAM = bpparam("SnowParam"),
              )
  
  flog.info("DESeq estimation successfull")
  
  res <-results(dds)
  res <- results(dds,contrast=c("f1","tumor","non-tumor"))
  

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
                    c("-o", "--output_dir"),
                    default = getwd(),
                    type = "character",
                    help = paste(c("directory to save output into.",
                                   "if none specified cwd will be used."),
                                 collapse = " "),
                    
                    )

parser <- add_option(parser,
                     c("-f","--feature_dir"),
                     help = paste(c("directory of feature-files.",
                                    "file name should be on form",
                                    '"count_data-ID_Replicate.tsv"'),
                                  collapse = " "),
                      )


parser <- add_option(parser,
                     c("-g", "--gene_file"),
                     default =  NULL,
                     help = paste(c("list of genes to be included in DGE analysis.",
                                    "each row shall contain name of gene.",
                                    'if no list provided all genes will be used.'),
                                  collapse = " "),
                     )


parser <- add_option(parser,
                     c("-d", "--design"),
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
                     type = "integer",
                     default = floor(detectCores()/2)
)

args <- parse_args(parser)

source(args$design)

if (!exists(x = "design")) {
  print("Please provide a proper design-file")
  print("Exiting")
  quit(0,status = 0)
}


#TODO : Consider whether default option for design should be includeded.

if (!interactive()) {
  
  COUNT_DATA_DIR <<- args$count_dir
  FEATURE_FILE_DIR <<- args$feature_dir
  OUTPUT_DIR <<- args$output_dir
  
  snowparam <- SnowParam(workers = args$workers, type = "SOCK")
  register(snowparam, default = TRUE)
  registered()
  
  
  flog.info(paste(c(args$workers, "cores are being used for computation"), collapse = " "))
    
  flog.threshold(DEBUG)
  #TODO : Consider if STDOUT or log-file should be used
  #flog.appender(appender.file(paste(getwd(),'test_evaluate.logger',sep='/')))
  
    
  if (! typeof(args$gene_file) == "logical") {
    
    SELECTED_GENES <<- read.csv(args$gene_file, header = FALSE)[,1]
    flog.info(paste(c("Using the genes specified in file", args$gene_file), collapse = " "))
  } else {
    SELECTED_GENES <<- FALSE
    flog.info("Using union of all genes.")
  }
  
  main()

} else {
  
  print("Run via terminal for proper activation.")
  
}