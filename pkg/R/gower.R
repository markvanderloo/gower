


#' Gower's distace
#'
#' Recycling over the data.frame with the least records.
#'
#' @param x data.frame
#' @param y data.frame
#' 
#' @export
gower_dist <- function(x, y, ranges){
  ilog = which(sapply(x,is.logical))
  icat = which(sapply(x,is.factor))
  inum = which(sapply(x,is.numeric))
  
  if (missing(ranges)){
    ranges = sapply(inum, function(i){ 
      r <- diff(range(c(range(x[,i],na.rm=TRUE),range(y[,i],na.rm=TRUE))))
      if (r>0) 
        r 
      else 
        pmax( max(x[,i],na.rm=TRUE), max(y[,i],na.rm=TRUE), na.rm=TRUE )
    })
  }
  .Call(R_gower,x,y,ilog,icat, inum, as.numeric(ranges))
}
