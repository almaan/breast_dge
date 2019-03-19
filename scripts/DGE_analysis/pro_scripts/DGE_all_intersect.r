pth <- "/home/alma/ST-2018/CNNp/DGE/res/DGEresults_all/dgefiles.dat"
files <- read.delim(pth, sep ="\n",header = F, stringsAsFactors = F)
files <- as.character(files$V1)

genelist <- list()
testnames <- c()
for (ii in 1:length(files)) {
  genes <- c()
  suc <- try(genes <- read.csv(files[ii], row.names = 1, header = T, stringsAsFactors = F))
  id <- gsub('\\.fancy\\.tsv','',basename(files[ii]))

  if (class(suc) == "try-error") {
    suc2 <- try(genes <- read.csv(files[ii], row.names = 1, header = T, stringsAsFactors = F, sep = "\t"))
    
    if (!(class(suc2) == "try-error")) {
      genes <- rownames(genes)
      genelist[[id]] <- genes
      testnames <- c(testnames,id)
    } else {
        next
    }
  } else {
    genes <- rownames(genes)
    genelist[[id]] <- genes
    testnames <- c(testnames,id)
  }
  
}

intermat <- matrix(0, nrow = length(genelist), ncol = length(genelist))

for (ii in 1:(length(genelist)-1)){
  intermat[ii,ii] <- length(genelist[[ii]])
  s <- 1+ii
  for (jj in s:length(genelist)){
    inter <- length(intersect(genelist[[ii]],genelist[[jj]]))
    print(c(ii,jj,inter))
    intermat[ii,jj] <- inter
    intermat[jj,ii] <- inter
  }
}

colnames(intermat) <- names(genelist)
rownames(intermat) <- names(genelist)
df <- as.data.frame(intermat)
