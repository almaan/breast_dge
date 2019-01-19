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
    print(pos)
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
