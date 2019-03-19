#!/usr/bin/Rscript

sh <- suppressPackageStartupMessages
sh(library(edgeR))
sh(library(optparse))
sh(library(org.Hs.eg.db))

# Functions ---------------------

get_contrast_combinations <- function(x) {
  xs <- sort(x)
  ncomb <- choose(length(xs),2)
  pos <- c()
  neg <- c()
  for (p in 1:(ncomb-1)){
    for (n in (p+1):ncomb) {
      pos <- c(pos,xs[p])
      neg <- c(neg,xs[n])
    }
  }
  return(data.frame(pos = pos, neg = neg))
  
}

# Parser ---------------------

parser <- OptionParser()

parser <- add_option(parser, c("-r","--res_file"),
                     type ="character",
                     help = c(""))

parser <- add_option(parser, c("-t","--tag"),
                     type = "character",
                     default = "tumor_id",
                     help = c(""))

parser <- add_option(parser, c("-p","--pvalue"),
                     type = "double",
                     help = c(""))


parser <- add_option(parser, c("-a","--get_all"),
                     type ="logical",
                     default = F,
                     action = "store_true",
                     help = c(""))

args <- parse_args(parser)

# Variable space -------------

tag <- args$tag
pval = args$pvalue
pth <- args$res_file
bname <- gsub("r\\.res\\.|\\.r","",basename(pth))
odir <- dirname(pth)
load(pth)

# Main -----------------------

cnames <- colnames(fit_dge)
pre_contrast <- rep(0,length(cnames))
tumor_pos <- grep(tag,cnames)

if (args$get_all){
  comb <- expand.grid(tumor_pos,tumor_pos)
  colnames(comb) <- c("pos","neg")
  comb <- comb[!(comb$pos == comb$neg),]
} else {
  comb <- get_contrast_combinations(tumor_pos)
}

for (x in 1:dim(comb)[1]) {
    p <- comb$pos[x]
    n <- comb$neg[x]
    crt <- rep(0,length(cnames))
    crt[c(p,n)] <- c(1,-1)
    lrt <- glmLRT(fit_dge, contrast = crt)
    fancy <- topTags(lrt, n = dim(lrt)[1], adjust.method = "BH", p.value= pval)
    try(fancy$symbol <- mapIds(org.Hs.eg.db, keys=rownames(fancy),column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
    fname <- paste(c(bname,".",cnames[p],"_vs_",cnames[n],".fancy.tsv"), collapse = "")
    oname <- paste(c(odir,fname),collapse = "/")
    write.csv(fancy,oname)
}

for (x in tumor_pos) {
  lrt <- glmLRT(fit_dge, coef = x)
  fancy <- topTags(lrt, n = dim(lrt)[1], adjust.method = "BH", p.value= pval)
  try(fancy$symbol <- mapIds(org.Hs.eg.db, keys=rownames(fancy),column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
  fname <- paste(c(bname,".",cnames[x],"_vs_non_tumor",".fancy.tsv"), collapse = "")
  oname <- paste(c(odir,fname),collapse = "/")
  write.csv(fancy,oname)
}
