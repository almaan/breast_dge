library(ggplot2)
library(latticeDensity)
library(spatgraphs)
library(network)
library(igraph)
library(gridExtra)
library(sp)
source("/home/alma/ST-2018/CNNp/DGE/scripts/DGE_analysis/lib/zone_generation.r")
load <- TRUE
if (load) {
  cpth <- "/home/alma/ST-2018/CNNp/DGE/data/count_data/her2lum/count_data-23287_C2.tsv"
  pth <- "/home/alma/ST-2018/CNNp/DGE/data/curated_feature_files/her2lum/23287_C2.feature_file.tsv"
  dgpth <- "/home/alma/ST-2018/CNNp/DGE/res/DGEresults/section_within_patient/DGE_analysis.2019-01-31-21-39-11.71906.fancy.tsv"
  fm <- read.csv(file = pth, header = TRUE, sep = "\t", row.names = 1)
  cmat <- read.csv(cpth, sep = "\t",header = TRUE, row.names = 1)
  inter <- intersect(rownames(fm),rownames(cmat))
  fm <- fm[inter,]
  cmat <- cmat[inter,]
  dgm <- read.csv(dgpth, sep = ",",header = TRUE, row.names = 1)
  dgm <- dgm[order(dgm$FDR),]
  xlab <- "xcoord"
  ylab <- "ycoord"
  feature <- "tumor"
}

#G <- spatgraph(fm[c(xlab,ylab)], "geometric", par = 1.5)
#net <- graph_from_adj_list(G$edges)
#net <- as.undirected(net)
#net$weight <- gene.val
#plot.igraph(net, layout = as.matrix(fm[c(xlab,ylab)]), vertex.size = net$weight, edge.width = 1, vertex.label = NA)
#fm$zones <- make_zones(fm[c(xlab,ylab)],labels = fm$tumor, zone_method = "mult_levels", ulim = 3, llim = 0)

#spy1 <- fm[c(xlab,ylab)]
#spy1 <- spy1[chull(spy1),] 

fm$tumor <- relevel(fm$tumor, "tumor")
cmap <- c("red","green")
names(cmap) <- levels(fm$tumor)

N = 50

#idx_pos <- intersect(c((nrow(dgm)-N+1):nrow(dgm)),which(dgm$logFC > 0))
#idx_neg <- intersect(c(1:N),which(dgm$logFC < 0))
idx_pos <- which(rownames(dgm) == "ENSG00000123358")
idx_neg <- which(rownames(dgm) == "ENSG00000123358")

#idx_neg <- which(dgm$logFC < quantile(dgm$logFC[idx_neg],0.2))

#alpha.neg <- rowSums(cmat[,rownames(dgm)[idx_neg]])
#alpha.neg <- alpha.neg/max(alpha.neg)
#alpha.pos <- rowSums(cmat[,rownames(dgm)[idx_pos]])
#alpha.pos <- alpha.pos/max(alpha.pos)
#alpha.pos[is.nan(alpha.pos)] <- 0
#alpha.neg[is.nan(alpha.neg)] <- 0
#alpha.pos.size <- 4
alpha.pos <- cmat[["ENSG00000123358"]]
alpha.pos[is.na(alpha.pos)] <- 0
alpha.pos <- alpha.pos/max(alpha.pos)
alpha.pos[is.na(alpha.pos)] <- 0
alpha.neg <- alpha.pos * 0


s.size <- 4

p1 <- ggplot(data = fm, aes(x = xcoord,y = ycoord)) + 
  #            geom_polygon(data = spy1, fill = NA, color = "black", alpha = 0.2 ) +
              geom_point(size = s.size, color = rgb(r=alpha.pos,g=0,b=alpha.neg, alpha = 0.8)) + 
             ggtitle(paste0(c(length(idx_pos),"most overexpressed and",length(idx_neg), "\n most underexpressed genes"),collapse = " "))
#p2 <- ggplot(data = fm, aes(x = xcoord,y = ycoord)) +



p3 <- ggplot(data = fm, aes(x = xcoord, y = ycoord)) + 
      geom_point(size = s.size,aes(color = tumor)) +
      scale_color_manual(values = cmap) +
      theme(legend.position = "none")+
      ggtitle("Pathologists annotations")
#
  
  
grid.arrange(p1, p3, ncol=2)
g <- arrangeGrob(p1,p3, ncol =2)
m <- unlist(regexec("[0-9]{5}.[A-Z][0-9]", cpth))
id <- substr(cpth,m,m+7)
opth <- "/home/alma/ST-2018/CNNp/DGE/res/DGEresults/zone_within_subtype/analysis/"
ggsave(filename=paste(c(opth,id,'dge.res.png'),collapse =""),g)
