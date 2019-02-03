library(org.Hs.eg.db)
library(igraph)
library(cluster)
library(anocva)
load("/home/alma/ST-2018/CNNp/DGE/data/Robjects/matrices.r")

DGE_PTH <- "/home/alma/ST-2018/CNNp/DGE/res/DGEresults/section_and_tumor_within_all/DGE_analysis.2019-02-01-14-25-38.40231.fancy.tsv"

logfc <- "logFC"
pval <- "FDR"
N <- 1000


dgm <- read.csv(DGE_PTH, sep = ",",header = TRUE, row.names = 1)
dgm <- dgm[order(dgm[logfc],decreasing = FALSE),]

pos_idx <- which(dgm[logfc] > 0 )
pos_idx <- pos_idx[1:min(length(pos_idx),N)]

neg_idx <- which(dgm[logfc] < 0 )
neg_idx <- neg_idx[1:min(length(neg_idx),N)]

genes_pos <- rownames(dgm)[pos_idx]
genes_neg <- rownames(dgm)[neg_idx]

keep_genes <- c(genes_pos, genes_neg)

matrices$count_matrix <- matrices$count_matrix[,keep_genes]

keep.samples <- rowSums(matrices$count_matrix) > 0


matrices$count_matrix <- matrices$count_matrix[keep.samples,]
matrices$feature_matrix <- matrices$feature_matrix[keep.samples,]
matrices$feature_matrix[] <-lapply(matrices$feature_matrix,factor)
matrices$count_matrix <- as.data.frame(matrices$count_matrix)



rownames(matrices$feature_matrix) <- rownames(matrices$count_matrix)

adj <- cor(matrices$count_matrix, method = "pearson")
over_thrs <- abs(adj) > 0.8
adj[over_thrs] <- 1
adj[!over_thrs] <- 0
diag(adj) <- 0
adj <- adj[-which(rowSums(adj) == 0),-which(colSums(adj) == 0)]

gr <- graph_from_adjacency_matrix(adj, mode ="undirected")
plot(gr, vertex.size=5, vertex.label=NA)
clu <- components(gr)
grp <- groups(clu)
grp <- grp[sapply(grp, length) > 5]

mapped <- mapIds(org.Hs.eg.db, keys=unlist(grp),
                 column="SYMBOL", keytype="ENSEMBL", multiVals="first")

gr <- induced_subgraph(gr,unlist(grp))
clu <- components(gr)

plot.igraph(gr,
            layout =layout_with_fr,
            vertex.size=20, 
            vertex.label=mapped,
            vertex.label.dist = 3,
            vertex.color = clu$membership)


