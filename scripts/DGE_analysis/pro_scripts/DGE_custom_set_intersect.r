#!/usr/bin/Rscript
sh <- suppressPackageStartupMessages
sh(library(org.Hs.eg.db))
sh(library(AnnotationDbi))
sh(library(ggplot2))
sh(library(optparse))
sh(library(RColorBrewer))

fisher_test <- function(g1,g2,complete_set) {
    ct <-matrix(c(
    length(setdiff(complete_set,union(g1,g2))),
    length(setdiff(g1,g2)),
    length(setdiff(g2,g1)),
    length(intersect(g1,g2))),
    nrow=2)
  res <- fisher.test(ct,alternative = "greater")
  pval <- res$p.value
  return(res$p.value)
}

parser <- OptionParser()

parser <- add_option(parser,
                     c("-f","--file_list"),
                     default = NULL,
                     help =""
)

args <- parse_args(parser)
 

AIMS <- "/home/alma/ST-2018/CNNp/DGE/data/custom_gene_sets/aims.dat"
HPA <- "/home/alma/ST-2018/CNNp/DGE/data/custom_gene_sets/HPA.prognostic.breast.dat"
mammaprint <- "/home/alma/ST-2018/CNNp/DGE/data/custom_gene_sets/mammaprint.dat"
oncotypedx <- "/home/alma/ST-2018/CNNp/DGE/data/custom_gene_sets/oncotype.dat"
pam50 <- "/home/alma/ST-2018/CNNp/DGE/data/custom_gene_sets/pam50.dat"
allgenes_pth <- "/home/alma/ST-2018/CNNp/DGE/data/custom_gene_sets/allgenes.in.dataset.symbols.dat"

custom_gene_sets_pth <- list(AIMS =AIMS,
                         HPA.Progonsitc = HPA,
                         MammaPrint = mammaprint,
                         OncoTypeDX = oncotypedx,
                         PAM50 = pam50)

all_genes <- read.delim(allgenes_pth, sep = '\n', col.names = c('genes'))


custom_gene_sets <- list()

for (geneset in names(custom_gene_sets_pth)) {
  custom_gene_sets[[geneset]] <- read.delim(custom_gene_sets_pth[[geneset]], sep = '\n', col.names = c('genes'))
}

dge_results_paths <- read.delim(args$file_list, sep = '\n', col.names = 'files', stringsAsFactors = F)
dge_results_paths <- dge_results_paths$files

dge_fisher_mat <- data.frame(matrix(0, ncol = length(dge_results_paths), nrow = length(custom_gene_sets)))
colnames(dge_fisher_mat) <- as.character(sapply(dge_results_paths, function(x) basename(x)))
rownames(dge_fisher_mat) <- names(custom_gene_sets)

for (dgeres in dge_results_paths) {
  
  dgegenes <- read.table(dgeres, sep = ',',
                         header =1,
                         row.names = 1,
                         stringsAsFactors = F)
  
  if (!('symbol' %in% colnames(dgegenes))) {
    dgegenes$symbol <- mapIds(org.Hs.eg.db,
                            keys=as.character(dgeres$genes),
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
  }
  dge_fisher <- sapply(custom_gene_sets, function(x) fisher_test(g1 = as.character(dgegenes$symbol),g2 = as.character(x$genes), complete_set = as.character(all_genes$genes)))
  dge_fisher_mat[,basename(dgeres)] <- dge_fisher  
}

dge_fisher_mat[,'nullref'] <- rep(0.01, length(custom_gene_sets))
coul = colorRampPalette(brewer.pal(8, "YlOrRd"))(25)
png('custom_set_intersect_res.png', width = 1000, height = 1000)
heatmap(-log(as.matrix(dge_fisher_mat)),
        Colv = NA,
        Rowv = NA,
        scale = "none",
        col = coul,
        margins = c(30,10),
        cexRow = 1,
        cexCol = 1)
dev.off()