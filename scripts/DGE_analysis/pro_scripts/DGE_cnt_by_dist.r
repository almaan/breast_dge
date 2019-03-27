g#!/usr/bin/Rscript

library(spatstat)
library(viridis)
library(optparse)
library(org.Hs.eg.db)
library(AnnotationDbi)


# Wrapper for multiple argument parsing -----------------

allowMultipleArgs <- function(){
  
  #' Modify trailing arguments passed such that space
  #' separated arguments to same flag becomes joined by
  #' commas; a format supported by optparse, and which later 
  #' easily can be split into separate parts again
  
  
  oriArgs <- commandArgs(trailingOnly = TRUE)
  flags.pos <- which(sapply(oriArgs, function(x) '-' == substr(x,1,1)))
  newArgs <- c()
  
  if (length(flags.pos) > 1) {
    for (i in 1:(length(flags.pos)-1))
    {
      if ((flags.pos[i] + 1) != flags.pos[i+1]) {
        pos <- c((flags.pos[i]+1):(flags.pos[i+1]-1))
        newArgs <- c(newArgs,oriArgs[flags.pos[i]], paste(oriArgs[pos],collapse=','))
      } else {
        newArgs <- c(newArgs,oriArgs[flags.pos[i]])
      }
    }
  }
  
  if (length(oriArgs) > tail(flags.pos,n=1)) {
    pos <- c((flags.pos[length(flags.pos)]+1):length(oriArgs))
    newArgs <- c(newArgs, oriArgs[tail(flags.pos,n=1)],paste(oriArgs[pos],collapse=','))
  } else {
    newArgs <- c(newArgs, oriArgs[tail(flags.pos,n=1)])
  }
  return(newArgs)
}


splitMultipleArgs <- function(optArgs) {
  
  #' Use in combination with allowMultipleArgs
  #' will split all commaseparated arguments
  #' into individual elements in list
  
  for (i in 1:length(optArgs)) {
    if (grepl(",",optArgs[[i]])) {
      optArgs[[i]] <- unlist(strsplit(optArgs[[i]],','))
    }
  }
  
  return(optArgs)
}

# Parser ------------------------------

parser <- OptionParser()

parser <- add_option(parser,
                     c("-d","--dge_res"),
                     default = NULL,
                     help =""
)


parser <- add_option(parser,
                     c("-c","--count_files"),
                     default = NULL,
                     help = ""
)

parser <- add_option(parser,
                     c("-f","--feature_files"),
                     default = NULL,
                     help = ""
)

parser <- add_option(parser,
                     c("-p","--polydeg"),
                     default = 5,
                     help = ""
)

parser <- add_option(parser,
                     c("-t","--title"),
                     default = NULL,
                     help = ""
)


# Variables -----------------

if (interactive()) {
  # example set when running in interactive mode to test performance
  # ajdust paths to make compatible with local version
  
  st_cnt_files <- c('/home/alma/ST-2018/CNNp/DGE/data/count_data/her2nonlum/count_data-23567_D2.tsv',
                    '/home/alma/ST-2018/CNNp/DGE/data/count_data/her2nonlum/count_data-23567_E2.tsv',
                    '/home/alma/ST-2018/CNNp/DGE/data/count_data/her2nonlum/count_data-23567_E1.tsv')
  
  feat_files <- c('/home/alma/ST-2018/CNNp/DGE/data/curated_feature_files/her2nonlum/23567_D2.feature_file.tsv',
                  '/home/alma/ST-2018/CNNp/DGE/data/curated_feature_files/her2nonlum/23567_E2.feature_file.tsv',
                  '/home/alma/ST-2018/CNNp/DGE/data/curated_feature_files/her2nonlum/23567_E1.feature_file.tsv')

  genes_pth <- '/home/alma/ST-2018/CNNp/DGE/res/DGEresults2/section_and_zone_within_subtype_new_d2/her2nonlum/DGE_analysis.2019-03-19-15-11-50.64798.tme_vs_tumor.fancy.tsv'
  
  polydeg <- 5
  
  main_title <- "test"

} else {
  # For CLI usage
  args <- splitMultipleArgs(parse_args(parser, args = allowMultipleArgs())) 
  st_cnt_files <- args$count_files
  feat_files <- args$feature_files
  genes_pth <- args$dge_res
  polydeg <- args$polydeg
  main_title <- args$title
}

# Validate input -----------------------

# check for consitency in number of files
if(!(length(st_cnt_files) == length(feat_files))) {
  print('Unequal number of count and feature files')
  quit(status = 1)
}

# makes sure count and feature fiels are matched
for (filenum in c(1:length(feat_files))) {
  pat <- "[0-9]{5}_[A-Z][0-9]" # section id pattern
  # grep for section specific pattern in filenames
  frgx <- regexpr(pat,feat_files[filenum],perl =T)
  crgx <- regexpr(pat,st_cnt_files[filenum],perl =T)
  # extract match
  len <- attr(frgx, 'match.length')
  fsec <- substr(feat_files[filenum],frgx[1],frgx[1] + len)
  csec <- substr(st_cnt_files[filenum],crgx[1],crgx[1] + len)
  
  # exit if count and feature file section id does not match  
  if (!(fsec == csec)) {
    print('Not properly mathced files')
    quit(status =1)
  }
}

print('>> Matching feature files and count files')
print('>> Will Initiate analysis')

# Main -----------------------

# load genes to be used in analysis 
genes <- read.table(genes_pth, sep = ',',
                    header =1,
                    row.names = 1,
                    stringsAsFactors = F)

genes <- genes$genes # get ENSEMBL ids
genes <- genes[!(is.na(genes))] # remove any potential NA

mx <- 50 #  maximum distance from distance from tumor spot
dr <- 0.1 #  distance increment
lims <- seq(0,mx+dr,dr) #  lower limits 
meanpnt <- lims + dr/2.0 #  mean point within each interval

# dataframe to hold mean expression values for each distance value
data <- data.frame(matrix(0,ncol = length(genes), nrow = length(lims)),
                   row.names = paste(lims,(lims + dr),sep = '-'))
colnames(data) <- genes

# get gene symbol names  
symbol <- mapIds(org.Hs.eg.db,
                 keys=as.character(genes),
                 column="SYMBOL",
                 keytype="ENSEMBL",
                 multiVals="first")

# if symbols are not found use ENSEMBL id
symbol[is.na(symbol)] <- genes[is.na(symbol)]

# iterate over all files
for (filenum in c(1:length(st_cnt_files))){
  
  st_cnt_pth <- st_cnt_files[filenum]
  feat_pth <- feat_files[filenum]
  
  feat <- read.table(feat_pth, sep = '\t', header = 1, row.names = 1, stringsAsFactors = F)
  cnt_raw <- read.table(st_cnt_pth, sep = '\t', header = 1, row.names = 1, stringsAsFactors = F)
  cnt_raw <- cnt_raw / rowMeans(cnt_raw)
  
  # extract spots present in feature and count file
  interspt <- as.character(intersect(rownames(feat),rownames(cnt_raw)))
  cnt_raw <- cnt_raw[interspt,]
  feat <- feat[interspt,]
  
  # create matrix with all as columns genes
  cnt <- data.frame(matrix(0, nrow = nrow(cnt_raw),ncol = length(genes)))
  colnames(cnt) <- genes
  rownames(cnt) <- interspt
  
  # interection of genelist and genes in st-count data
  inter <- as.character(intersect(genes,colnames(cnt_raw)))
  
  # fill count matrix with scaled st-counts
  cnt[,inter] <- cnt_raw[,inter]
  remove(cnt_raw) # to free up space
  
  cnt <- cnt[,inter]
  
  # get tumor and non-tumor spot coordinates
  tmrspts <- feat[feat$tumor == 'tumor',c('xcoord','ycoord')]
  nonspot <- feat[feat$tumor == 'non',c('xcoord','ycoord')]
  
  # create distance matrix (euclidian distance)
  dmat = crossdist.default(X = nonspot[,1], Y = nonspot[,2], x2 = tmrspts[,1], y2 = tmrspts[,2])
  rownames(dmat) = rownames(nonspot)
  colnames(dmat) = rownames(tmrspts)
  
  # get minimum distance to a tumor spot for each non-tumor spot
  mindist <- apply(dmat, 1, min)
  
  for (ii in c(1:length(lims))) {
    # find spots with mindist within I = [lim,lim+dr)
    idx_within <- (mindist >= lims[ii]) & (mindist < lims[ii] + dr) 
    # get count matrix index for spots within I = [lim,lim+dr)
    spt_within <- which(rownames(cnt) %in% rownames(dmat)[idx_within]) #  
    # if I is non-empty
    if (length(spt_within) > 0) {
      # add mean of observed counts taken over points within I
      data[ii,] <- data[ii,] + colMeans(cnt[spt_within,]) 
    }
  }

}

# adjust for multiple sections having been used
dataf <- data / length(st_cnt_files)
# remove genes and distances where no observetions are made
dataf <- dataf[rowSums(dataf) > 0,colSums(dataf) > 0]
# remove distances where no spots are present
meanpntf <- meanpnt[rowSums(dataf) > 0]
# name of analysis
bname <- gsub("\\.fancy\\.tsv","",basename(genes_pth))
print(paste(c(dirname(genes_pth),paste(c(bname,"count_by_distance.png"),collapse = '.')),collapse = '/'))
# save result to png-file
png(paste(c(dirname(genes_pth),paste(c(bname,"count_by_distance.png"),collapse = '.')),collapse = '/'),
    width = 1000,
    height = 500)

# panel for plot and legend
par(mfrow=c(1,2))
# generate colormap
cmap <- viridis(dim(data)[2])

# get plot title
if (is.null(main_title)) {
  plt_title <- main_title
}

# Plot Results -------------------------

# initiate plot
plot(meanpntf,dataf[,1],
     ylim = c(0,max(dataf + 1)),
     xlab = 'distance to nearest tumor spot',
     ylab = 'mean expression',
     col = cmap[1],
     main = gsub('DGE_analysis\\.','',bname),
     cex.main = 0.9,
     pch = 19
     )

# add points for all genes
for (ii in 2:(ncol(dataf))) {
  points(x= meanpntf, dataf[,ii], col = cmap[ii], pch = 19)
}

# add fitted curves for each gene using same color as points
fitted <- list()
for (ii in c(1:ncol(dataf))) {
  # using polynomial of degree "polydeg" to fit data
  fitted[[genes[ii]]] <- lm(formula = y ~ poly(x, polydeg, raw=TRUE),
                            data = data.frame(x = meanpntf, # distance is independent variable
                                              y = dataf[,ii]) # mean expression is dependent variable
                            )
  
  lines(meanpntf, predict(fitted[[genes[ii]]], data.frame(x=meanpntf)), col=cmap[ii], lwd = 5)
}

# plot legend in different subplot for readability
plot(x=NULL, xlim = c(0,2), ylim = c(0,1), xlab = '', ylab = '')
legend(0,1,legend= c(symbol), fill = cmap, cex = 0.5, bty = 'n', ncol = 5  )
dev.off()


