library(ggplot2)

make_zones <- function(fm,
                       dr = 1.0, 
                       method = "euclidian",
                       feature_name = "tumor",
                       foci_label = "tumor",
                       env_label = "non",
                       xlab = 'xcoord',
                       ylab = 'ycoord')  {
  
  idx_foci <- fmat[[feature_name]] == foci_label
  idx_env <- fmat[[feature_name]] == env_label
  
  zones <- matrix(0,length(fmat[[feature_name]]))
  dm <- as.matrix(dist(fmat[c(xlab,ylab)]))
  max.dist <- ceiling(max(dm[idx_foci,idx_env]))
  withinCutoff <- 1
  reduceCutoff <- 0.0
  
  while (sum(withinCutoff) > 0  & (max.dist - reduceCutoff) > 0) {
    withinCutoff <- rowSums(dm[,idx_foci] < (max.dist - reduceCutoff)) >= 1
    if (sum(withinCutoff) > 0) {
      print(max.dist - reduceCutoff)
      zones[idx_env & withinCutoff] <- (max.dist - reduceCutoff) 
    }
    reduceCutoff <- reduceCutoff + dr
  }
  return(zones)
}

visualize_zones <- function(fm,
                            feature_name = "tumor",
                            foci_label = "tumor",
                            env_label = "non",
                            xlab = 'xcoord',
                            ylab = 'ycoord') {
  
  alph <- as.numeric(fm[[feature_name]] == foci_label)
  s.size <- 5.0

  plot.obj <- ggplot(fm, aes(x = fm[[xlab]], y = fm[[ylab]])) + 
    
    geom_point(size = s.size*1.2, alpha = alph, color = "black") +
    geom_point(aes(color = zones), size = s.size) +
    scale_color_gradient(low = "red", high = "green") + 
    xlab("x-coordinates") + 
    ylab("y-coordinates")
  
  return(plot.obj)
  
}
