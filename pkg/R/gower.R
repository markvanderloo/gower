


#' Gower's distace
#'
#' Recycling over the data.frame with the least records.
#'
#' @param x data.frame
#' @param y data.frame
#' 
#' @export
gower_dist <- function(x, y){

  pair <- match(names(x),names(y),nomatch = 0L)
  factor_pair <- sapply(x,is.factor)
  .Call(R_gower,x,y,pair, as.integer(factor_pair))
  
}
