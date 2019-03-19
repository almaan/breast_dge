#!/usr/bin/Rscript

library(gplots)
library(viridis)
library(optparse)
library(heatmaply)
  
make_heatmap_res <- function(x,
                             stem_name,
                             Rowv = T,
                             Colv = T,
                             colors_dir= 1,
                             dendrogram = "none",
                             margins = c(60,60)) {
  
  #' Takes an interaction matrix as 
  
  X <- log1p(x)
  
  pdf(paste(c(stem_name,"pdf"),collapse = "."),
      width = 65,
      height = 65)
  
  
  
  heatm <- heatmap.2(X,
                  scale = "none",
                  Rowv = Rowv,
                  Colv = Colv,
                  key = F,
                  hclustfun = function(w){hclust(dist(w),method="ward.D2")},
                  dendrogram = dendrogram,
                  margins = margins,
                  trace = "none",
                  col = viridis(100, direction = colors_dir),  
                  na.color = 'black',
                  cexRow = 2.3,
                  cexCol = 2.3,
                  offsetRow = 1.0)

  title(stem_name)
  dev.off()
  
  sv <- X[rev(heatm$rowInd),heatm$colInd]
  hmly <-   heatmaply(sv,
              margins = c(60,60,60,60),
              Rowv = NA,
              Colv = NA,
              col = viridis(100, direction = colors_dir),
              file = paste(c(stem_name,".html"),collapse =""))
  rm(hmly)
  
  return(heatm)
}

mi_sub = function(x) {
  p = prop.table(x)
  rs = rowSums(p)
  cs = colSums(p)
  indep = rs %*% t(cs)
  logratios = log2(p / indep)
  logratios[x==0] = 0
  return(p * logratios)
  
}

mi = function(x) {
  return(sum(mi_sub(x)))  
  
}

make_interact_matrix <- function(set_list,
                                 set_names,
                                 stat_test = F,
                                 complete_set = NA) {
  n_sets <- length(set_list)
  intermat <- matrix(0, nrow = n_sets, ncol = n_sets)
  for (xx in 1:(n_sets-1)){
    intermat[xx,xx] <- length(set_list[xx])
    s <- 1+xx
    for (jj in s:n_sets){
      g1 <- strsplit(as.character(set_list[xx]),",")[[1]]
      g2 <- strsplit(as.character(set_list[jj]),",")[[1]]
      if (any(c(g1,g2) == "none")) {
        inter <- NA
      } else {
        inter <- length(intersect(g1,g2))
        
        if (is.character(stat_test)) {
          ct <-matrix(c(
                        length(setdiff(complete_set,union(g1,g2))),
                        length(setdiff(g1,g2)),
                        length(setdiff(g2,g1)),
                        inter),
                      nrow=2)
          
          if (stat_test == 'fisher') {
            inter <- fisher.test(ct,alternative = "greater")
            inter <- -log(inter$p.value)
          } else if (stat_test == 'mi') {
            inter = -log(mi(ct))
          }
        }  
      }
      intermat[xx,jj] <- inter
      intermat[jj,xx] <- inter
      
    }
  }
  
  intermat[is.infinite(intermat)] <- max(intermat[!((is.infinite(intermat)) | (is.na(intermat)))])
  rownames(intermat) <- set_names
  colnames(intermat) <- set_names
  
  return(intermat)
}


prs <- OptionParser()

prs <- add_option(prs, c("-rd","--res_dir"))
prs <- add_option(prs, c("-i","--include"))
prs <- add_option(prs, c("-rf","--res_file"))
prs <- add_option(prs, c("-o","--output_dir"))
prs <- add_option(prs, c("-t","--tag"),default = '')
prs <- add_option(prs, c("-s","--complete_set"))
prs <- add_option(prs, c("-c", "--conditioned_on"),default='')
prs <- add_option(prs, c("-ts", "--tested_for"),default='')


args <- parse_args(prs)

dir_pth <- args$res_dir
dge_pth <- args$include
res_pth <- args$res_file
comp_pth <- args$complete_set
opath <- args$output_dir
tag <- args$tag



ressum <- read.csv(res_pth,sep='\t', header = T, row.names = 1, stringsAsFactors = F)
filenames <- read.delim(dge_pth,sep = "\n", header = F, stringsAsFactors = F)
complete_set <- as.character(read.delim(comp_pth,sep = ",", header = F, stringsAsFactors = F))
complete_set <- unique(as.character(sapply(complete_set,
                                          function(x) strsplit(x,',')[[1]])))


colnames(filenames) <- "names"
inter <- intersect(as.character(rownames(ressum)),as.character(filenames$names))
filenames <- filenames[which(filenames$names %in% inter),]
ressum <- ressum[inter,]

if (grepl("\\w",args$conditioned_on,perl=T)) {
  to_include <- grepl(paste(c(args$conditioned_on,"|none"),collapse=""),ressum$conditioned_on,perl=T)
  ressum <- ressum[to_include,]
  filenames <- filenames[to_include]
  print('adjusted for conditioning')
}


if (grepl("\\w",args$tested_for,perl=T)) {
  to_include <- grepl(paste(c(args$tested_for,"|none"),collapse=""),ressum$tested_for,perl=T)
  ressum <- ressum[to_include,]
  filenames <- filenames[to_include]
  print('adjusted for testing')
}

if (!all(rownames(ressum) == filenames)) {
  print("ERROR: Non-matching indices")
  quit()
}

# Make Interpretable names -------------------

new_names <- sapply(c(1:dim(ressum)[1]), function(x)
  paste(c("T",ressum[x,'tested_for'],"A",ressum[x,'accounted_for'],"C",ressum[x,"conditioned_on"]),collapse  = '_'))

new_names <- as.character(new_names)

new_names_down <- sapply(c(1:dim(ressum)[1]), function(x)
  paste(c("T",ressum[x,'tested_for'],"A",ressum[x,'accounted_for'],"C",ressum[x,"conditioned_on"],"D"),collapse  = '_'))

new_names_up <- sapply(c(1:dim(ressum)[1]), function(x)
  paste(c("T",ressum[x,'tested_for'],"A",ressum[x,'accounted_for'],"C",ressum[x,"conditioned_on"],"U"),collapse  = '_'))

new_names_lfc <- as.character(c(new_names_down,new_names_up))

cnames <- c("total","all_genes","negative","positive","neg_genes","pos_genes")
statmat <- data.frame(matrix(0,length(filenames), length(cnames)), row.names = filenames)
colnames(statmat) <- cnames

print('initate stats matrix generation')
# Generate Dataframe with analysis information -----------------
for (ii in filenames){
  file <- paste(c(dir_pth,ii,list.files(paste(c(dir_pth,ii),collapse = "/"),pattern ="fancy.tsv")[1]),collapse = "/")
  a <- try(res <- read.csv(file, sep = ",", header = T, row.names = 1))
  if(!(class(a) == "try-error")) {
    statmat[ii,"total"] <- dim(res)[1]
    statmat[ii,"all_genes"] <- paste(rownames(res),collapse =",")
    if ('logFC' %in% colnames(res)) {
      statmat[ii,"negative"] <- sum(res$logFC < 0)
      statmat[ii,"positive"] <- sum(res$logFC > 0)
      statmat[ii,"neg_genes"] <- paste(rownames(res[res$logFC < 0,]),collapse =",")
      statmat[ii,"pos_genes"] <- paste(rownames(res[res$logFC > 0,]),collapse =",")
    } else {
      statmat[ii,3:4] <- rep(0,2)
      statmat[ii,5:length(cnames)] <- 'none'
    }
  } else {
    statmat[ii,3:4] <- rep(0,2)
    statmat[ii,5:length(cnames)] <- 'none'
  }
}
print("generated stats matrix")
write.table(statmat, file = paste(c(opath,"DGE_heatmap_data.tsv"),collapse = "/"),
          sep = "\t", quote = F, row.names = T, col.names = T)
print("saved stats matrix")

#complete_set <- unique(as.character(sapply(paste(statmat$all_genes,collapse =","),
#                                          function(x) strsplit(x,',')[[1]])))


# Generate interaction Matrices -------------------

sets_all <- statmat['all_genes']
sets_lfc <- as.data.frame(matrix(0,nrow = length(new_names_up)*2),
                                ncol = 1,
                                row.names = cnames)

print("created empty sets frames")

sets_lfc[1:dim(statmat)[1],] <- as.character(statmat$neg_genes)
sets_lfc[(1+dim(statmat)[1]):(2*dim(statmat)[1]),] <- as.character(statmat$pos_genes)
colnames(sets_lfc) <- "lfc_genes"

print('generated gene set information')

imat_list <- list()
name_list <- list()
tests <- c("none","fisher","mi")
for (tst in c(1:length(tests))) {
  imat_list[[2*tst -1]] <- make_interact_matrix(sets_all$all_genes,new_names,stat_test = tests[tst], complete_set = complete_set)
  imat_list[[2*tst]] <- make_interact_matrix(sets_lfc$lfc_genes,new_names_lfc,stat_test = tests[tst], complete_set = complete_set)
  name_list[[2*tst-1]] <- paste(c(opath,paste(c(tag,paste(c("DGE_heatmap_non_directional",tests[tst]),collapse="_")),collapse =".")),collapse="/")
  name_list[[2*tst]] <- paste(c(opath,paste(c(tag,paste(c("DGE_heatmap_directional",tests[tst]),collapse="_")),collapse =".")),collapse="/")
  
}

print('constructed name lists and interaction matrices')
# Extremely Ugly Plottting ---------------


use_dendro <- rep("none",length(name_list))
use_dendro[1:2] <- 'both'
use_cols <- rep(NA,length(name_list))
use_cols[1:2] <- T
color_dir <- rep(1,length(name_list))
color_dir[(length(name_list)-1):length(name_list)] <- -1
  
hm <- list()

for (ii in 1:length(name_list)) {
  x <- imat_list[[ii]]
  ridx <- !(rowSums(is.na(x)) >= (dim(x)[2] -1))
  cidx <- !(colSums(is.na(x)) >= (dim(x)[1] -1))
  x <- x[ridx,cidx]
  
  if (ii > 2){
    ridx <- rev(hm[[2-(ii %% 2)]]$rowInd)
    cidx <- hm[[2-(ii %% 2)]]$colInd
    x <- x[ridx,cidx]
  }  
  
  hm[[ii]] <- make_heatmap_res(x, 
                               name_list[[ii]],
                               Rowv = use_cols[ii],
                               Colv = use_cols[ii],
                               colors_dir = color_dir[ii],
                               dendrogram = use_dendro[ii]
                         )
}