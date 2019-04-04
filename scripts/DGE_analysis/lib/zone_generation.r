library(ggplot2)

make_zones <- function(crd,
                       labels,
                       zone_method = "mult_levels",
                       dist_method = "euclidian",
                       dr = 1.0, 
                       foci_label = "tumor",
                       llim = 0,
                       ulim  = NULL,
                       verbose = TRUE)  {
  
#' generates discrete zones of spots based on their
#' distance to the nearest spot of a given label. Useful for
#' selection of microenvironment.
#' 
#' Two methods are available, "mult_levels" and "three_levels".
#' 
#' mult_levels: 
#' assigns labels to spots within all dicrete intervals of size dr. All spots
#' not labeled with foci_label will be given a labeled based on their distance to
#' nearest foci_label spot.
#' 
#' three_levels:
#' will create only three levels "0: foci_label", "1:interface" and "2:peripheral"
#' all spots with a distance to nearest foci_label spit within the
#' interval (llim,ulim] are labeled as interface whilst remaining non foci_label
#' spot are labeled as peripheral   
#' 
#' @param : crd - N x 2 vector or matrix containing spot coordinates
#' @param : labels - N x 1 vector containg labels of feature of interest
#' @param : zone_method - string, choose from "mult_levels" or "three_levels". See above for
#'          explanation
#' @param : dist_method  - string, method used to calculate distance between spots. Default is
#'          euclidian. Available are also "maximum", "manhattan", "canberra", "binary" or 
#'          "minkowski.
#' @param : dr - numerical, increment for distance when using method mult_levels
#' @param : foci_label - label of foci to which distance should be computed to
#' @param : llim - numerical, lower limit to be used in "three_levels" method
#' @param : ulim - numerical, upper limit to be used in "three_levels" method
#' @param : verbose - bool, activates verbose mode. Default True.
#' 
#' @returns : zones - N x 1 vector with zone labels
  
  # print method
  if (verbose) {print(paste(c("using method >>",zone_method), collapse = " "))}
  
  labels <- unlist(labels)
  
  # make sure limits are provided
  if (zone_method == "three_levels" & is.null(ulim)) {
      if (verbose){print(paste(c("upper and lower limit needs to be specified when using",
                    'method "three_levels". Do this using arguments llim and',
                    'ulim for lower respective upper limit'), collapse = " "))}
      return(NA)
  
  } else {
  
    # get bool indices of foci spots
    idx_foci <- labels == foci_label
    # prepare vector
    zones <- matrix("tumor",nrow = length(labels))
    # compute distance matrix between spots
    dm <- as.matrix(dist(crd), method = dist_method , diag = TRUE) #  diag true for proper size
    
    # output based on method
    if (zone_method == "three_levels") {
      zones[!(idx_foci)] <- "distal"
      idx_inter <- (!(idx_foci) & (rowSums(dm[,idx_foci] <= ulim) >=1) & !(rowSums(dm[,idx_foci] < llim) >=1))
      idx_below_inter <- (!(idx_foci) & (rowSums(dm[,idx_foci] < llim) >=1))
      # only assign if intermediary zones
      if (any(idx_inter)) {
        zones[idx_inter] <- "micro"
      }
      
      if (any(idx_below_inter)){
        # assign dummy index to region between intermediary and tumor zone
        zones[idx_below_inter] <- "inter"
      }
    
    if (zone_method == "tme_vs_all") {
      zones[,] <- "distal_and_tumor"
      idx_inter <- (!(idx_foci) & (rowSums(dm[,idx_foci] <= ulim) >=1) & !(rowSums(dm[,idx_foci] < llim) >=1))
      idx_below_inter <- (!(idx_foci) & (rowSums(dm[,idx_foci] < llim) >=1))
      # only assign if intermediary zones
      if (any(idx_inter)) {
        zones[idx_inter] <- "micro"
      }
      
      if (any(idx_below_inter)){
        # assign dummy index to region between intermediary and tumor zone
        zones[idx_below_inter] <- "inter"
      }
        
    } else if (zone_method == "mult_levels"){
      
      max.dist <- ceiling(max(dm[idx_foci,!(idx_foci)]))
      withinCutoff <- 1
      reduceCutoff <- 0.0
      
      while (sum(withinCutoff) > 0  & (max.dist - reduceCutoff) > 0) {
        withinCutoff <- rowSums(dm[,idx_foci] < (max.dist - reduceCutoff)) >= 1
        if (sum(withinCutoff) > 0) {
          zones[!(idx_foci) & withinCutoff] <- (max.dist - reduceCutoff) 
        }
        reduceCutoff <- reduceCutoff + dr
      }
      
    } else {
      if(verbose){print('ERROR: Specified method is not valid')}
      zones <- NA
    }
    return(zones)
  }
}

visualize_zones <- function(crd,
                            labels,
                            zones,
                            foci_label = "tumor"
                            ) {
  
  #' Visualize the generated zones from the function 'generate_zones'
  #' Will plot all spots with colors corresponding to distance to
  #' nearest foci_label  spot, on the spectrum between red to green.
  #' Spots annotated as foci_label will be marked with black border
  #' 
  #' @param : crd - N x 2 vector or matrix containing spot coordinates
  #' @param : labels - N x 1 vector containg labels of feature of interest
  #' @param : foci_label - label of foci to which distance should be computed to
  #' 
  #' @return : plot.obj - a grob (ggplot) object 
  
  # set border to zero for non foci_label spots  
  alph <- as.numeric(labels == foci_label)
  
  # genrate plot
  s.size <- 5.0 
  plot.obj <- ggplot(crd, aes(x = crd[,1], y = crd[,2])) + 
    geom_point(size = s.size*1.2, alpha = alph, color = "black") + #  black border on feature_label spots
    geom_point(aes(color = zones), size = s.size) + #  color by zone
    scale_color_gradient(low = "red", high = "green") + #  set color scheme
    xlab("x-coordinates") +
    ylab("y-coordinates")
  
  return(plot.obj)
  
}
