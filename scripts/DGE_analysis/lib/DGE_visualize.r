#!/usr/bin/env Rscript
sh <- suppressPackageStartupMessages

sh(library(ggplot2))
sh(library(optparse))
sh(library(gridExtra))
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))

LIB_PATH <- "~/ST-2018/CNNp/DGE/scripts/DGE_analysis/lib"

source(paste(c(LIB_PATH,"modOptparse.r"),collapse ="/"))

# Parser ----------------------------


parser <- OptionParser()

parser <- add_option(parser,c("-d","--dge_result"),
                     type = "character",
                     help = paste(c("Differential gene Expression",
                                    "file either generated using",
                                    "either generated using DESeq",
                                    "or edgeR."),
                                  collapse =" ")
                     )

parser <- add_option(parser,c("-f","--feature_file"),
                     type = "character",
                     help = paste(c("feature file containing",
                                    "tumor annotations and",
                                    "potential other usefull",
                                    "information"),
                                  collapse =" ")
)


parser <- add_option(parser,c("-c","--count_file"),
                     type = "character",
                     help = paste(c("count file containing",
                                    "for given sections"),
                                  collapse =" ")
)


parser <- add_option(parser,c("-m","--method"),
                     type = "character",
                     help = paste(c("method by which the",
                                    "DGE file was generated",
                                    'either "deseq2" or "edger"'),
                                  collapse =" ")
)

parser <- add_option(parser,c("-t","--topgenes"),
                     type = "integer",
                     default = 20,
                     help = paste(c("number of top/bottom genes",
                                    "to use for relative frequencies",
                                    "default is set to 20"),
                                  collapse =" ")
)


parser <- add_option(parser,c("-o","--output"),
                     type = "character",
                     default = NA,
                     help = paste(c("number of top/bottom genes",
                                    "to use for relative frequencies",
                                    "default is set to 20"),
                                  collapse =" ")
)


args <- splitMultipleArgs(parse_args(parser,
                                     args = allowMultipleArgs()))


# Main ------------------------------

fancy_file <- args$dge_result
n_sections <- length(args$feature_file)
args$method <- tolower(args$method)

args$feature_file <- sort(args$feature_file)
args$count_file <- sort(args$count_file)


if (args$method == "deseq") {
  pval_col <- "padj"
  logfc <- "log2FoldChange"
} else if (args$method == "edger") {
  pval_col <- "FDR"
  logfc <- "logFC"
}


dgm_ori <- read.csv(args$dge_result, sep = ",",header = TRUE, row.names = 1)
print(head(dgm_ori))
dgm_ori <- dgm_ori[order(dgm_ori[pval_col],decreasing = FALSE),]

plot_dge_list <- list()
plot_list <- list()
#par(mfrow = c(n_sections,2))

for (section in c(1:n_sections)) {
  
  feature_file <- args$feature_file[section]
  count_file <- args$count_file[section]
  
  fm <- read.csv(file = feature_file, header = TRUE, sep = "\t", row.names = 1)
  cmat <- read.csv(count_file, sep = "\t",header = TRUE, row.names = 1)
  
  inter <- intersect(rownames(fm),rownames(cmat))
  
  fm <- fm[inter,]
  cmat <- cmat[inter,]
  dgm <- dgm_ori[intersect(rownames(dgm_ori),colnames(cmat)),]
  
  print(length(intersect(rownames(dgm_ori),colnames(cmat))))
  
  xlab <- "xcoord"
  ylab <- "ycoord"
  feature <- "tumor"
  
  fm[xlab] <- round(fm[xlab])
  fm[ylab] <- round(fm[ylab])
    
  fm$tumor <- relevel(fm$tumor, "tumor")
  cmap <- c("red","green")
  names(cmap) <- levels(fm$tumor)
  
  N <- args$topgenes
  
  pos_idx <- which(dgm[logfc] > 0 )
  neg_idx <- which(dgm[logfc] < 0 )
  
  print(length(pos_idx))
  print(length(neg_idx))
  
  if (length(pos_idx) > 0) {
    pos_idx <- pos_idx[1:min(length(pos_idx),N)]
    genes_pos <- rownames(dgm)[pos_idx]
    rel.freq.pos <- rowSums(cmat[,genes_pos])
    rel.freq.pos <- rel.freq.pos/sum(rel.freq.pos)
    rel.freq.pos <- rel.freq.pos / max(rel.freq.pos)
  } else {
    rel.freq.pos <- rep(0,dim(fm)[1])    
  }
  
  if (length(neg_idx) > 0) {
    neg_idx <- neg_idx[1:min(length(neg_idx),N)]  
    genes_neg <- rownames(dgm)[neg_idx]
    rel.freq.neg <- rowSums(cmat[,genes_neg])
    rel.freq.neg <- rel.freq.neg/sum(rel.freq.neg)
    rel.freq.neg <- rel.freq.neg / max(rel.freq.neg)
  } else {
    rel.freq.neg <- rep(0,dim(fm)[1])
  }
  
  print(colnames(fm))
  print(fm[[xlab]][1:10])
  print(n_sections)
    
  plot_list[[section]] <- ggplot(data = fm, aes_string(x = xlab,y = ylab)) + 
            geom_tile(aes_string(x = xlab, y = ylab), 
              fill=rgb(r=rel.freq.pos,b=0.0,g=rel.freq.neg, alpha = 0.8), interpolate = T, color = "black") + 
            ggtitle(paste(c(
                    basename(args$count_file[section]),
                    "\n",
                    length(genes_pos),
                    "top.sig overexpressed\n",
                    length(genes_neg),
                    "top.sig underexpressed"),
                    collapse = " "))

  plot_list[[n_sections + section]] <- ggplot(data = fm, aes_string(x = xlab, y = ylab)) + 
            geom_tile(aes_string(x = xlab, y = ylab, fill = feature), color = "black") + 
            #geom_point(size = s.size, color = 'black', alpha = 0.2) +
            scale_color_manual(values = cmap) +
            theme(legend.position = "none") +
            ggtitle(paste(c(basename(args$feature_file[section]),
                          "\n",
                          "Pathologists annotations\n"),
                          collapse = " "))
  
}

#all_plots <- append(plot_dge_list, plot_ano_list)
g <- grid.arrange(grobs = plot_list, ncol=n_sections, nrow = 2)
if (is.na(args$output)) {
  X11()
  plot(g)
  cat("Press [enter] to continue ...")
  fil <- readLines(con="stdin", 1)
  cat(fil, "\n")

  } else {
  
  if (grepl("\\.pdf",args$output)) {
    output_name <- args$output
  } else {
    output_name <- paste(c(args$output,".pdf"),collapse = "")
  } 
  ggsave(output_name, g, width = 297, height = 210, units = "mm")

}
