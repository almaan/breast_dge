#!/usr/bin/Rscript

library(edgeR)
library(optparse)
library(org.Hs.eg.db)
library(AnnotationDbi)

parser <- OptionParser()

parser <- add_option(parser, c("-r","--res_file"),
                     type ="character",
                     help = c(""))

args <- parse_args(parser)

ipth <- args$res_file
odir <- dirname(ipth)
odir <- paste(c(odir,"comp"), collapse = "/")
bname <- gsub("r\\.res\\.|\\.r","",basename(ipth))

load(ipth)


cnames <-colnames(fit_dge)
len <- length(cnames)
base <- rep(-1.0/(len-1),len)

for (i in 1:len) {
  print(paste(c("Working with >> ",cnames[i], " | sample ",i, "/",len),collapse = ""))
  contrast <- base
  contrast[i] <- 1.0
  start <- Sys.time()
  lrt <- glmLRT(fit_dge, contrast = contrast)
  end = Sys.time()
  print("Fitting Complete")
  print(end-start)
  
  fancy <- topTags(lrt, n = dim(lrt)[1], adjust.method = "BH", p.value= 0.01)
  try(fancy$symbol <- mapIds(org.Hs.eg.db, keys=rownames(fancy),column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
  
  fname <- paste(c(bname,".",cnames[i],"_vs_all.fancy.tsv"), collapse = "")
  oname <- paste(c(odir,fname),collapse = "/")
  print("Successfully saved file!")
  write.csv(fancy,oname)
}
  