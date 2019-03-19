library(gplots)

pth <- "/home/alma/ST-2018/CNNp/DGE/res/DGEresults_all/logfccounts2.tsv"
pth_map <- "/home/alma/ST-2018/CNNp/DGE/res/DGEresults_all/results_summary.tsv"

sets <- read.csv(pth,sep='\t', header = T, row.names = 1 )
ressum <- read.csv(pth_map,sep='\t', header = T, row.names = 1, stringsAsFactors = F)
sep1 <- rep("T:",dim(ressum)[1])
sep2 <- rep("A",dim(ressum)[1])
sep3 <- rep("C",dim(ressum)[1])

new_names <- sapply(c(1:dim(ressum)[1]), function(x)
                      paste(c("T",ressum[x,'tested_for'],"A",ressum[x,'accounted_for'],"C",ressum[x,"conditioned_on"]),collapse  = '_'))
new_names <- as.data.frame(new_names,row.names = rownames(ressum))

intermat <- matrix(0, nrow = dim(sets)[1], ncol = dim(sets)[1])


for (xx in 1:(dim(sets)[1]-1)){
  print(xx)
  intermat[xx,xx] <- length(sets$all_genes[xx])
  s <- 1+xx
  for (jj in s:dim(sets)[1]){
    g1 <-strsplit(as.character(sets$all_genes[xx]),",")[[1]]
    g2 <- strsplit(as.character(sets$all_genes[jj]),",")[[1]]
    inter <- length(intersect(g1,g2))
    print(c(xx,jj,inter))
    intermat[xx,jj] <- inter
    intermat[jj,xx] <- inter
  }
}

cn <-as.character(sapply(rownames(sets), function (x) 
    paste(tail(strsplit(x,"/")[[1]],ifelse(2,1,length(strsplit(x,"/")[[1]]) > 2),collapse = "."))))

cn <-  gsub("\\.fancy\\.tsv","",cn)
colnames(intermat) <- new_names[cn,1]
rownames(intermat) <- new_names[cn,1]

df <- as.data.frame(intermat)
pdf("/home/alma/ST-2018/CNNp/DGE/res/heatmaps/heatmap_all_gens.pdf",width = 40, height = 40)
heatmap.2(log1p(as.matrix(df)),scale = "none", Rowv = T, Colv = T, dendrogram = "none", key = F,
        margins = c(60,60),trace = "none",cexRow = 2.3, cexCol = 2.3, offsetRow = 1.0)

dev.off()


new_names_down <- sapply(c(1:dim(ressum)[1]), function(x)
  paste(c("T",ressum[x,'tested_for'],"A",ressum[x,'accounted_for'],"C",ressum[x,"conditioned_on"],"D"),collapse  = '_'))

new_names_down <- as.data.frame(new_names_down,row.names = rownames(ressum))


new_names_up <- sapply(c(1:dim(ressum)[1]), function(x)
  paste(c("T",ressum[x,'tested_for'],"A",ressum[x,'accounted_for'],"C",ressum[x,"conditioned_on"],"U"),collapse  = '_'))

new_names_up <- as.data.frame(new_names_up,row.names = rownames(ressum))

new_names_down <- as.character(new_names_down[cn,1])
new_names_up <- as.character(new_names_up[cn,1])


dirdf <- as.data.frame(matrix(0,nrow = length(new_names_up)*2),ncol = 1,row.names = as.character(c(new_names_down,new_names_up)))
dirdf[1:dim(sets)[1],] <- as.character(sets$neg_genes)
dirdf[(1+dim(sets)[1]):(2*dim(sets)[1]),] <- as.character(sets$pos_genes)



intermat_lfc <- matrix(0, nrow = dim(dirdf)[1], ncol = dim(dirdf)[1])

for (xx in 1:(dim(dirdf)[1]-1)){
  intermat_lfc[xx,xx] <- length(dirdf[xx,1])
  s <- 1+xx
  for (jj in s:dim(dirdf)[1]){
    g1 <-strsplit(as.character(dirdf[xx,1]),",")[[1]]
    g2 <- strsplit(as.character(dirdf[jj,1]),",")[[1]]
   
    inter <- length(intersect(g1,g2))
    print(c(xx,jj,inter))
    intermat_lfc[xx,jj] <- inter
    intermat_lfc[jj,xx] <- inter
  }
}

rownames(intermat_lfc) <- rownames(dirdf)
colnames(intermat_lfc) <- rownames(dirdf)

pdf("/home/alma/ST-2018/CNNp/DGE/res/heatmaps/heatmap_all_genes_lfc.pdf",width = 65, height = 65)
heatmap.2(log1p(intermat_lfc),scale = "none", Rowv = T, Colv = T, dendrogram = "none", key = F,
          margins = c(60,60),trace = "none",cexRow = 2.3, cexCol = 2.3, offsetRow = 1.0)

dev.off()
