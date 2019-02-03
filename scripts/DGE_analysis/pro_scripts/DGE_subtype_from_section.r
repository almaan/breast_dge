#!/usr/bin/Rscript
library(edgeR)
library(optparse)
library(futile.logger)

parser <- OptionParser()

parser <- add_option(parser, 
                     c("-d","--dge_file"),
                     type = "character",
                     help = c("name of dge result file."))

parser <- add_option(parser, 
                     c("-p","--p_cutoff"),
                     type = "double",
                     default = 0.01,
                     help = c("p-value cutoff"))


parser <- add_option(parser, c("-o", "--output"),
                            type = "character",
                            default = NA,
                            help = c("specify name of output file. If none, default is used"))

args <- parse_args(parser)

flog.appender(appender.tee(paste(c(dirname(args$dge_file),"subtype_from_section.log"),
                                 collapse="/")),name = "ROOT")


flog.info(sprintf("using DGE-file >>  %s", args$dge_file))

load("/home/alma/ST-2018/CNNp/DGE/data/Robjects/subtype_patient_list.r") #  load subtype object
load(args$dge_file)

flog.info(sprintf("The following subtype association is used %s : %s", names(subtype_list), subtype_list))

cnames <- colnames(fit_dge)

flog.info(sprintf("The follwoing colnames are used >> %s ", paste(cnames, collapse =", ")))

subtypes <- names(subtype_list)


for (subtype in subtypes) {
  ctr <- replicate(length(cnames),0)
  idx_sub <- grepl(paste(subtype_list[[subtype]],collapse ="|"), cnames)
  idx_other <- !idx_sub
  ctr[idx_sub] <- 1/sum(idx_sub)
  ctr[idx_other] <- -1/sum(idx_other)
  
  if (is.na(args$output)) {
    oname <- paste(c(dirname(args$dge_file),paste(c(subtype,gsub('r.res.|\\.r','',basename(args$dge_file)),"fancy.tsv"),collapse = ".")), collapse = "/")
  } else {
    oname <- paste(c(dirname(args$output),paste(c(subtype,basename(args$output)),collapse = ".")), collapse = "/")  
  }
  
  flog.info(sprintf("Will save subtype %s results to file %s",subtype, oname))
  flog.info(sprintf("Will be using contrast >> %s", paste(ctr,collapse = ",")))
  lrt <- glmLRT(fit_dge, contrast = ctr)
  flog.info("glmLRT fitted completed")
  lrt$table$FDR <- p.adjust(lrt$table$PValue, method = "BH")
  lrt$table <- lrt$table[lrt$table$FDR < args$p_cutoff,]
  flog.info("adjusted Pvalues added to results using method . Filtered for adj.Pvalues <= %f", args$p_cutoff)
  write.table(lrt$table, file=oname,sep ="\t")  
}

