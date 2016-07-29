


#' Gower's distance
#'
#' Compute Gower's distance, pairwise between records in two data sets \code{x} 
#' and \code{y}. Records from the smallest data set are recycled over.
#' 
#' @section Details:
#' There are three ways to specify which columns of \code{x} should be compared
#' with what columns of \code{y}. The first option is do give no specification. 
#' In that case columns with matching names will be used. The second option
#' is to use only the \code{pairs_y} argument, specifying for each column in \code{x}
#' in order, which column in \code{y} must be used to pair it with (use \code{0} 
#' to skip a column in \code{x}). The third option is to explicitly specify the
#' columns to be matched using \code{pair_x} and \code{pair_y}.
#' 
#' 
#' @section Note:
#' Gower (1971) originally defined a similarity measure (\eqn{s}, say)
#' with values ranging from 0 (completely dissimilar) to 1 (completely similar).
#' The distance returned here equals \eqn{1-s}.
#'
#'
#' @param x \code{[data.frame]}
#' @param y \code{[data.frame]}
#' @param pair_x \code{[numeric|character] (optional)} Columns in \code{x} used for comparison. 
#'    See Details below.
#' @param pair_y \code{[numeric|character] (optional)} Columns in \code{y} used for comparison. 
#'    See Details below.
#' @param eps = \code{[numeric] (optional)} Computed numbers (variable ranges) 
#'    smaller than \code{eps} are treated as zero. 
#' 
#' 
#' @return
#'   A \code{numeric} vector of length \code{max(nrow(x),nrow(y))}.
#' 
#' @references 
#' Gower, John C. "A general coefficient of similarity and some of its 
#' properties." Biometrics (1971): 857-871.
#' 
#' @export
gower_dist <- function(x, y, pair_x=NULL, pair_y=NULL, eps = 1e-8){
  gower_work(x=x,y=y,pair_x=pair_x,pair_y=pair_y,n=NULL,eps=eps)
}

#' Find the top-n matches
#' 
#' @description
#' 
#' Find the top-n matches in \code{y} for each record in \code{x}.
#' 
#' @inheritParams gower_dist
#' @param n The top-n indices and distances to return.
#' 
#' @seealso \code{\link{gower_dist}}
#' 
#' @return 
#'  A \code{list} with two array elements: \code{index}
#'  and \code{distance}. Both have size \code{n X nrow(x)}. Each ith column 
#'  corresponds to the top-n best matches of \code{x} with rows in \code{y}.
#' 
#' @examples 
#' # find the top 4 best matches in the iris data set with itself.
#' x <- iris[1:3,]
#' lookup <- iris[1:10,]
#' gower_topn(x=x,y=lookup,n=4)
#' 
#' 
#' @export
gower_topn <- function(x, y, pair_x=NULL, pair_y = NULL, n=5, eps=1e-8){
  gower_work(x=x,y=y,pair_x=pair_x,pair_y=pair_y,n=n,eps=eps)
}



gower_work <- function(x, y, pair_x=NULL, pair_y = NULL, n=NULL, eps=1e-8){
  stopifnot(is.numeric(eps), eps>0)
  
  
  if (is.null(pair_x) & is.null(pair_y)){
    pair <- match(names(x),names(y),nomatch = 0L)
  } else if (is.null(pair_x)){
    pair <- pair_y
  } else {
    if (is.character(pair_x) & is.character(pair_y)){
      m <- match(names(x),pair_x,nomatch=0)
      pair_x <- pair_x[m]
      pair_y <- pair_y[m]
    }
    pair <- numeric(ncol(x))
    pair[pair_x] <- pair_y
  }
  factor_pair <- as.integer(sapply(x,is.factor))
  eps <- as.double(eps)
  pair <- as.integer(pair-1L)
  if (is.null(n)){
    .Call(R_gower, x, y , pair, factor_pair, eps)
    
  } else {
    L <- .Call(R_gower_topn, x, y, pair, factor_pair, as.integer(n), eps)
    names(L) <- c("index","distance")
    dim(L$index) <- c(n,nrow(x))
    dim(L$distance) <- dim(L$index)
    dimnames(L$index) <- list(topn=NULL,row=NULL)
    dimnames(L$distance) <- dimnames(L$index)
    L
  }
    
}











