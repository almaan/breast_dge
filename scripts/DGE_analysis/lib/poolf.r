library(dbscan)
#NOTE: Written as having the genes as rows and samples as columns
get_idx <- function(data,
                    niter,
                    k_neighbors,
                    lim){
  
  idx_list <- matrix(0,nrow = niter, ncol = k_neighbors + 1 )
  nn <- kNN(as.matrix(data),
            k = k_neighbors)
  
  if (dim(nn$id)[1] >= niter) {
    choice <- sample(c(1:dim(nn$id)[1]),niter)
  } else {
    choice <- c()
    while(length(choice) < niter) {
      choice <- c(choice,sample(c(1:dim(nn$id)[1]),
                                min(niter-dim(nn$id)[1],dim(nn$id)[1])))
    }
  }
  
  for (ii in 1:length(choice)){
    inlim <- as.vector(nn$dist[choice[ii],] <= lim)
    n_inlim <- sum(inlim)
    
    if (length(n_inlim) == 0) {
      idx <- replicate(k_neighbors,c(choice[ii]))
    } else if(n_inlim < k_neighbors ) {
      idx <- c(choice[ii],as.vector(nn$id[choice[ii],inlim]))
      idx <- c(idx,sample(idx, k_neighbors - n_inlim, replace =TRUE))
    } else {
      idx <- c(choice[ii],nn$id[choice[ii],])
    }
    idx_list[ii,] <- idx
  }
  
  return(idx_list)
}

get_ratio <- function(n_samples, tabb){
  ratio <- sapply(names(tabb), function(x) round(n_samples*tabb[[x]]/sum(tabb)))
  if (sum(ratio) > n_samples) {
    ratio[sample(1:length(tabb))] = ratio[sample(1:length(tabb))] -1
  } else if (sum(ratio) < n_samples) {
    ratio[sample(1:length(tabb))] = ratio[sample(1:length(tabb))] +1
  }
  return(as.vector(ratio))
}

make_pseudo <- function(cnt, ft, select, k_neighbors,lim, n_samples, transpose = TRUE){
  if(transpose){
    cnt <- t(cnt)
  }
  
  tab <- sort(table(ft[[select]]),decreasing = TRUE)
  features <- names(tab)
  
  ratio <- get_ratio(n_samples = n_samples,tab)
  pseudo_mat <- data.frame(matrix(0,nrow = dim(cnt)[1], ncol = n_samples))
  rownames(pseudo_mat) <- rownames(cnt)
  colnames(pseudo_mat) <-c(1:n_samples)
  pseudo_feat <-data.frame(0,row = n_samples, ncol = 1)
  pseudo_feat <- data.frame(unlist(as.vector(sapply(1:length(features),
                            function(x)  replicate(ratio[x],features[x])))))
  
  colnames(pseudo_feat) <- select
  bot = 0
  for (jj in 1:length(features)){
    pos <-ft[select] == features[jj]
    idx <- get_idx(cbind(ft$xcoord[pos],ft$ycoord[pos]),
                   niter = ratio[jj],
                   k_neighbors = k_neighbors,
                   lim = lim)
    
    for (ii in c(1:ratio[jj])) {
      pseudo_mat[,bot + ii] <-rowSums(cnt[,pos][,idx[ii,]])
    }
    bot <- ratio[jj]
  }
  
  if(transpose){
    pseudo_mat <- data.frame(t(pseudo_mat))
  }
  return(list(pseudo_cnt = pseudo_mat, pseudo_feat = pseudo_feat))
}
