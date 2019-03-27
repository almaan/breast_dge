#!/usr/bin/Rscript

#' Script used to analyze micoenvironment and distal tissue
#' Takes an r-object containing the edgeR output and performs a LRT-test
#' will contrast microenvironment against the tumor tissue by default
#' But can be changed as to compare against distal tissue as well
#' 
#' Assumptions
#' Zones are defined as
#' 0 - distal tissue
#' 1 - tumor micro environment
#' 2 - tumor

sh <- suppressPackageStartupMessages


sh(library(org.Hs.eg.db)) 
sh(library(AnnotationDbi))
sh(library(edgeR))
sh(library(optparse))


parser <- OptionParser()

parser <- add_option(parser, c('-f','--files'),
                     default = NULL,
                     help = paste(c('Directly Specify which',
                                  'filename(s) of R-objects',
                                  'to be used.'),collapse=' '
                                  ))

parser <- add_option(parser, c('-l','--list'),
                     default = NULL,
                     help = paste(c('Specify list of',
                                  'filenames of R-objects',
                                  'to be used. Convinient',
                                  'when multiple analyses have',
                                  'to be done'
                                  ),collapse = ' ')
)

parser <- add_option(parser, c('-t','--test'),
                     default = 'tumor',
                     help = paste(c('specify what region',
                                  'the TME should be',
                                  'compared to. Options are',
                                  'either "tumor" or "distal".',
                                  'Default is set to "tumor'
                                  ),collapse = ' ')
)

args <- parse_args(parser)
print(args$files)
print(args$list)

if (!is.null(args$files)) {
  print('Using Files Directly from argument')
  filenames <-args$files
} else if (!(is.null(args$list))) {
  print(sprintf('Reading r-object files from %s',args$list))
  filenames <- as.vector(unlist(read.delim(args$list,
                                           sep = '\n', 
                                           header = F,
                                           stringsAsFactors = F)))
  }



#loop over filenames as file
for (file in filenames) {
  print(basename(file))
  print(sprintf('Working with file %s', file))
  load(file)
  cnames <- colnames(fit_dge)
  print(cnames)
  tme <- which(grepl('zones1',cnames))
  tumor <-which(grepl('zones2',cnames))
  if (args$test == 'tumor') {
    print('Contrasting TME against Tumor Tissue')
    contrast <- rep(0,length(cnames))
    contrast[tme] <- 1
    contrast[tumor] <- -1
    lrt <- glmLRT(fit_dge, contrast = contrast)
    oname_type <- 'tme_vs_tumor'
  } else if (args$test == 'distal') {
    print(sprintf('Contrasting TME against Distal Tissue'))
    lrt <- glmLRT(fit_dge, coef = tme)
    oname_type <- 'tme_vs_distal'
  }
  
  topnres <- topTags(lrt, n = dim(lrt$table)[1], p.value = 0.01)
  symbols <-mapIds(org.Hs.eg.db,
                   keys=rownames(topnres$table),
                   column="SYMBOL",
                   keytype="ENSEMBL",
                   multiVals="first")
  
  fancy <- topnres$table
  fancy$symbol <- symbols
  print(sprintf('Identified %d signigicantly DE genes', dim(fancy)[1]))
  
  oname_base <- gsub('r.res.|.r','',basename(file))
  oname_filename <- paste(c(oname_base,oname_type,'fancy.tsv'),collapse = '.')
  oname_complete_path <- paste(c(dirname(file),oname_filename),collapse = '/')
  
  write.csv(x = fancy, file = oname_complete_path)
  print(sprintf('Wrote results to >> %s',oname_complete_path))
}
