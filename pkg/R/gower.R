


#' Gower's distace
#'
#' Recycling over the data.frame with the least records.
#'
#' @param x data.frame
#' @param y data.frame
#' 
#' @export
gower_dist <- function(x, y){
  ilog = which(sapply(x,is.logical))
  icat = which(sapply(x,is.factor))
  inum = which(sapply(x,is.numeric))
  
  .Call(R_gower,x,y,ilog,icat, inum)
}
