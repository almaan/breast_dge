sh <- suppressPackageStartupMessages
sh(library(dbscan))
sh(library(futile.logger))
sh(library(docstring))


#NOTE: Written as having the genes as rows and samples as columns
get_idx <- function(data,
                    n_samples,
                    k_members,
                    lim){
  #'Function to get indices of samples to pool from provided
  #'data. Indices will be given according to position in the
  #'data object.
  #'
  #'Arguments:
  #' data - matrix or vector with dimension (n_points, n_features)
  #' n_samples - number of pseudo-samples to generate
  #' k_members - number of members to include in each pseudo-sample
  #' lim - limit within which neighbours must reside. Negative for infinity.
  #'Returns:
  #' idx_mat - matrix with dimensions (n_samples, k_members) where indices
  #' of spots to be pooled in each sample are listed along second axis 
  
  idx_mat <- c()
  #if number of members are larger than the number
  #of data points adjust for this in KD-tree construction
  #one less neighbour than k_members is used to include spot
  adj_k_members <- min(k_members, dim(data)[1])
  nn <- kNN(as.matrix(data),
            k = adj_k_members-1)
  
  #choose n_samples random points, replace is less spots than n_samples
  doreplace <- (dim(nn$id)[1]) < n_samples
  if (doreplace) {
    flog.debug("More requested samples than spots with feature > Using Replacement")
  
  } else {
    flog.debug("Less requested samples than spots with feature > No Replacement")
  }
  choice <- sample(1:(dim(nn$id)[1]),size = n_samples, replace = doreplace)
  #TODO add alternative for purely random sampling not based on distance
  #get neighbour(s) of random points
  for (spot in choice){
   #find how many neighbours within limit that spot have
   inlim <- nn$dist[spot,] <= lim
   n_inlim <- sum(inlim)
   
   if (n_inlim == 0){
     #CASE: no neighbours are within limit 
     flog.debug("No neighbours within threshold > Replicating single spot")
     idx <- replicate(k_members,spot)
   
   } else if (inlim < k_members -1) {
     #CASE: has neighbours but less than specifed
     flog.debug(sprintf("%d neighbour(s) within threshold. > Sampling from neighbours",
               n_inlim,k_members))
     
     idx <- c(spot,nn$id[spot,inlim])
     idx <- c(idx, sample(idx, k_members - length(idx),replace = TRUE))
   
   } else {
     #CASE: has k_members within limit
     flog.debug('%d neighbours where within threshold > No Resampling required')
     idx <- c(spot,nn$id[spot,inlim])
   }
   idx_mat <- rbind(idx_mat,idx)
  }
  flog.debug(sprintf("Dimensions of index list are > (%d,%d) ", dim(idx_mat)[1], dim(idx_mat)[2]))
  return(idx_mat)
}

get_ratio <- function(n_samples, tabb){
  #'Function to scale number of datapoints to n_samples whilst
  #'maintaining ratio between different features
  #'Arguments:
  #' n_samples - number of pseudo-samples to be generated
  #' feat_tab - a table of features in dataset
  #' Returns:
  #' ratio - a vector with absolute number of samples to be drawn for each feature 

  #compute scaled values
  ratio <- sapply(names(tabb), function(x) round(n_samples*tabb[[x]]/sum(tabb)))
  
  #adjust for evenutal abundance/deficiency in samples from rounding
  if (sum(ratio) > n_samples) {
    ratio[sample(1:length(tabb))] = ratio[sample(1:length(tabb))] -1
  } else if (sum(ratio) < n_samples) {
    ratio[sample(1:length(tabb))] = ratio[sample(1:length(tabb))] +1
  }
  #adjust for zero assigned members in low abundance cases
  #TODO: make sure that this works for multiple features not just one
  if (any(ratio == 0)){
    pos <-(ratio== 0)
    ratio[pos] = ratio[pos] + 1
    ratio[!(pos)] = ratio[!(pos)] - 1
  }
  flog.debug(sprintf("the ratio between categories are %s", 
                     paste(rbind(ratio,names(ratio)), collapse = ":")))
  
  return(as.vector(ratio))
}

make_pseudo <- function(cnt, ft, select, k_members,lim, n_samples, transpose = TRUE){
  #'Function to make pseudo-count and feature matrices by pooling spots together
  #'Arguments:
  #' cnt - count matrix (assumed to be on form n_sample x n_genes)
  #' ft - feature matrix (must be aligned with cnt)
  #' select - feature to pool w.r.t.
  #' k_members - number of spots to be used in each pooling
  #' lim - distance within which neighbours of spots must fall in pooling
  #' n_samples - number of pseudo-samples to generate
  #' transpose - to adjust for orientation of count matrix
  #'Returns:
  #' pseudo_mat - list with two elements "psuedo_mat" and "pseudo_feat"  
  if(transpose){
    flog.debug('Count matrix is taken as n_samples x n_genes. Using Temporary transposition')
    cnt <- t(cnt)
  }
  
  tab <- sort(table(ft[[select]]),decreasing = TRUE)
  features <- names(tab)
  ratio <- get_ratio(n_samples = n_samples,tab)
  pseudo_mat <- data.frame(matrix(0,nrow = dim(cnt)[1], ncol = n_samples))
  pseudo_feat <-data.frame(matrix(0,nrow = n_samples, ncol = 1))
  rownames(pseudo_mat) <- rownames(cnt)
  rownames(pseudo_feat) <- colnames(pseudo_mat)
  colnames(pseudo_feat) <- select
  bot = 0

  for (jj in 1:length(features)){
    #get spots with specific feature annotation
    pos <-ft[[select]] == features[jj]
    flog.debug(sprintf("Sampling from spots with feature %s", features[jj]))
    
    #get indices of spots to pool
    idx <- get_idx(cbind(ft[['xcoord']][pos],ft[['ycoord']][pos]),
                   n_samples = ratio[jj],
                   k_members = k_members,
                   lim = lim)
    
    for (ii in c(1:ratio[jj])) {
      pseudo_mat[,bot + ii] <-rowSums(cnt[,pos][,idx[ii,]])
      pseudo_feat[[select]][bot + ii] <- features[jj]
      }
    bot <- bot + ratio[jj]
  }
  
  if(transpose){
    pseudo_mat <- data.frame(t(pseudo_mat))
  }
  return(list(pseudo_cnt = pseudo_mat, pseudo_feat = pseudo_feat))
}
