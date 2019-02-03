#!/usr/bin/Rscript

library(ggplot2)
library(gridExtra)

load("/home/alma/ST-2018/CNNp/DGE/data/Robjects/matrices.r")
load("/home/alma/ST-2018/CNNp/DGE/data/Robjects/subtype_patient_list.r")


# Functions ------------------

select_genes <- function(pvals, cutoff = 0.01, n_genes = 100) {
  idx <- pvals <= cutoff
  adj.num <- min(sum(idx),n_genes)
  idx <- order(pvals)[1:adj.num]
  return(idx)
}


# Settings --------------------
N_genes <- 1000
p_cutoff <- 0.01
logfc <- "logFC"
pval <- "FDR"

# Tumor within all  -------------------
tumor_within_all_pth <- "/home/alma/ST-2018/CNNp/DGE/res/DGEresults/section_and_tumor_within_all/DGE_analysis.2019-02-01-14-25-38.40231.fancy.tsv"
tumor_within_all <- read.csv(tumor_within_all_pth, sep = ",", head = TRUE, row.names = 1)

# Tumor within subtypes -----------------
subtype_path <- "/home/alma/ST-2018/CNNp/DGE/res/DGEresults/tumor_and_section_within_subtype"
subtype_files <- paste(subtype_path, list.files(subtype_path,pattern = "*fancy*"), sep = "/")
base_names <- sapply(list.files(subtype_path, pattern = "*fancy*"), function(x) gsub("fancy.tsv","",x))
log_names <- paste(subtype_path,paste(base_names,"log",sep = ""),sep="/")
subtypes <- names(subtype_list)

subtype_gene_sets <- list() 

for (filenum in 1:length(subtype_files)) {
  tmp <- read.csv(subtype_files[filenum], sep = ",", head = TRUE, row.names = 1)
  sbt <- names(which(sapply(subtypes, function(x) any(grepl(x,readLines(log_names[filenum]))))))
  subtype_gene_sets[[sbt]] <- rownames(tmp)
}

subtypes <- names(subtype_gene_sets)

# Compare ----------------

inter_tumor <- list()
inter_sub <- c()
n_sub <- 1
n_bus <- 1
mat <- matrix(0,nrow = length(subtypes), ncol = length(subtypes))

for (sub in subtypes) {
  inter_tumor[[sub]] <- length(intersect(unlist(subtype_gene_sets[sub]),unlist(rownames(tumor_within_all))))
  n_bus <- 1
  for (bus in subtypes) {
    mat[n_sub,n_bus] <- length(intersect(unlist(subtype_gene_sets[sub]),unlist(subtype_gene_sets[bus])))
    inter_sub <- c(inter_sub,ifelse(n_sub == n_bus, 0,mat[n_sub,n_bus]))
    n_bus <- n_bus +1 
  }  
  n_sub <- n_sub + 1
}

df <- data.frame(x = paste(sapply(subtypes, function(x) rep(x,length(subtypes)))), y = rep(subtypes, length(subtypes)),
                            vals = log(inter_sub))

mat <- as.data.frame(mat)
colnames(mat) <- subtypes
rownames(mat) <- subtypes

p1 <- ggplot(df, aes(x =x, y = y)) + 
      geom_tile(aes(fill = vals),color = "white")  + scale_fill_gradient(low = "white",
                                                     high = "steelblue")
p2 <- ggplot(data.frame(inter = inter_tumor, subtype = subtypes , total = sapply(subtype_gene_sets,function(x) length(x))), aes(x = subtypes, y = inter_tumor)) +
      geom_bar(stat ="identity",)

grid.arrange(p1,p2, ncol = 2)
