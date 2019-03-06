
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
#' @param eps \code{[numeric] (optional)} Computed numbers (variable ranges) 
#'    smaller than \code{eps} are treated as zero.
#' @param weights \code{[numeric] (optional)} A vector of weights of length \code{ncol(x)}
#'    that defines the weight applied to each component of the gower distance.
#' @param ignore_case \code{[logical]} Toggle ignore case when neither \code{pair_x}
#'    nor \code{pair_y} are user-defined. 
#' @param nthread Number of threads to use for parallelization. By default,
#'    for a dual-core machine, 2 threads are used. For any other machine 
#'    n-1 cores are used so your machine doesn't freeze during a big computation. 
#'    The maximum nr of threads are determined using \code{omp_get_max_threads} at C level.
#' 
#' 
#' @return
#'   A \code{numeric} vector of length \code{max(nrow(x),nrow(y))}.
#'   When there are no columns to compare, a message is printed and both
#'   \code{numeric(0)} is returned invisibly.
#' 
#' @seealso \code{\link{gower_topn}}
#' 
#' @references 
#' Gower, John C. "A general coefficient of similarity and some of its 
#' properties." Biometrics (1971): 857-871.
#' 
#' @export
gower_dist <- function(x, y, pair_x=NULL, pair_y=NULL, eps = 1e-8, weights=NULL,
                       ignore_case=FALSE, nthread=getOption("gd_num_thread")){
  check_recycling(nrow(x),nrow(y))
  gower_work(x=x,y=y,pair_x=pair_x,pair_y=pair_y
    , n=NULL, eps=eps, weights=weights, ignore_case=ignore_case, nthread=nthread)
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
#'  When there are no columns to compare, a message is printed and both
#'  \code{distance} and \code{index} will be empty matrices; the list is
#'  then returned invisibly.
#'
#' @examples 
#' # find the top 4 best matches in the iris data set with itself.
#' x <- iris[1:3,]
#' lookup <- iris[1:10,]
#' gower_topn(x=x,y=lookup,n=4)
#' 
#' 
#' @export
gower_topn <- function(x, y, pair_x=NULL, pair_y = NULL, n=5, eps=1e-8, weights = NULL,
                       ignore_case=FALSE, nthread=getOption("gd_num_thread")){
  gower_work(x=x,y=y,pair_x=pair_x,pair_y=pair_y
  , n=n, eps=eps, weights=weights, ignore_case=ignore_case, nthread=nthread)
}



gower_work <- function(x, y, pair_x, pair_y, n, eps, weights, ignore_case, nthread){
  stopifnot(is.numeric(eps), eps>0) 
  
  if (is.null(pair_x) & is.null(pair_y)){
    xnames <- if(ignore_case) toupper(names(x)) else names(x)
    ynames <- if(ignore_case) toupper(names(y)) else names(y)
    pair <- match(xnames, ynames, nomatch = 0L)
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
  if ( !any(pair > 0) ){
    message("Nothing to compare")
		return( if (is.null(n)){ # gower_dist
				invisible(numeric(0))
			} else { # gower_topn
				invisible(list(distance=matrix(0)[0,0],index=matrix(0L)[0,0]))
			}
		)
	}


  if ( !is.null(weights) && ( any( weights < 0 ) || !all(is.finite(weights)) ) ){
    stop("At least one element of 'weights' is not a finite nonnegative number"
       , call. = FALSE)
  }
  if ( !is.null(weights) && length(weights) < length(pair) ){
    msg <- sprintf("%d weights specified, expected %d"
                 , length(weights), length(pair))
    stop(msg, call. = FALSE)
  }
  # If the user didn't pass any weights, then weight all components of the 
  # distance equally.
  if (is.null(weights))
      weights <- rep(1, ncol(x))

  # check column classes

  nthread <- as.integer(nthread)
  ranges <- numeric(length(pair))
  for ( i in seq_along(pair)){
    if (pair[i] == 0 ) next
    ranges[i] <- .Call("R_get_xy_range",x[[i]],y[[pair[i]]],nthread)
  }

  factor_x <- sapply(x,is.factor)
  factor_y <- sapply(y,is.factor)
  
  for ( i in seq_along(pair) ){
    if ( pair[i] == 0 ) next
    iy <- pair[i]
    if (!factor_x[i] & !factor_y[iy]) next

		if ( factor_x[i] && !factor_y[iy] ){
			stop("Column ", i, " of x is of class factor while matching column "
          , pair[i]," of y is of class ", class(y[[iy]]) )
		}
		if ( !factor_x[i] && factor_y[iy] ){
			stop("Column ", i, " of x is of class ", class(x[[i]]), " while matching column "
          , pair[i]," of y is of class factor" )
		}
    if (factor_x[i] && factor_y[iy]){
			if ( !isTRUE( all.equal( levels(x[[i]]), levels(y[[iy]]) ) ) ){
        stop("Levels in column ", i, " of x do  not match those of column ", 
        pair[i], " in y.")

			}
		}
	}


  factor_pair <- as.integer(factor_x)
  
  eps <- as.double(eps)
  # translate to C-indices (base-0).
  pair <- as.integer(pair-1L)
  if (is.null(n)){
    .Call("R_gower", x, y , ranges, pair, factor_pair, eps, weights, nthread)
    
  } else {
    L <- .Call("R_gower_topn", x, y, ranges, pair, factor_pair, as.integer(n), eps, weights, nthread)
    names(L) <- c("index","distance")
    dim(L$index) <- c(n,nrow(x))
    dim(L$distance) <- dim(L$index)
    dimnames(L$index) <- list(topn=NULL,row=NULL)
    dimnames(L$distance) <- dimnames(L$index)
    L
  }
    
}

RECYCLEWARNING <- tryCatch(1:3+1:2,warning=function(e) e$message)

check_recycling <- function(nx,ny){
  mx <- max(nx,ny)
  mn <- min(nx,ny)
  if ((mx %% mn) != 0) warning(RECYCLEWARNING, call.=FALSE)
}








