mpth <- "/home/alma/ST-2018/CNNp/DGE/res/DGEresults_all/"
pth <- "/home/alma/ST-2018/CNNp/DGE/res/DGEresults2/ffiles.dat"
filenames <- read.delim(pth,sep = "\n", header = F, stringsAsFactors = F, col.names = 'files')
filenames <- filenames$files
#filenames <- t(filenames)
#filenames <- as.vector(filenames)


df <- data.frame(matrix(0,length(filenames), 6), row.names = unlist(filenames))
colnames(df) <- c("total",
                  "negative",
                  "positive"
                  #"all_genes",
                  #"neg_genes",
                  #"pos_genes"
                  )

for (ii in c(1:length(filenames))){
  #file <- paste(c(mpth,ii,list.files(paste(c(mpth,ii),collapse = "/"),pattern ="fancy.tsv")[1]),collapse = "/")
  #file <- paste(c(mpth,ii),collapse = "")
  a <- try(res <- read.csv(filenames[ii], sep = ",", header = T, row.names = 1, stringsAsFactors = F))
  if(!(class(a) == "try-error")) {
      df[ii,"total"] <- dim(res)[1]
      df[ii,"negative"] <- sum(res$logFC < 0)
      df[ii,"positive"] <- sum(res$logFC > 0)
 #     df[ii,"all_genes"] <- paste(rownames(res),collapse =",")
  #    df[ii,"neg_genes"] <- paste(rownames(res[res$logFC < 0,]),collapse =",")
   #   df[ii,"pos_genes"] <- paste(rownames(res[res$logFC > 0,]),collapse =",")
  } else {
    df[ii,] <- rep(0,length(colnames(df)))
  }
}

rownames(df) <- as.vector(sapply(rownames(df),function(x) basename(x)))
#df <- df[order(rownames(df),decreasing = T),]
write.table(df,"/home/alma/ST-2018/CNNp/DGE/res/DGEresults2/logfccounts2.tsv", sep ="\t",row.names = T, col.names = T, quote = F)

