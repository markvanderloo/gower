#' Gower's distance/similarity measure.
#'
#' A C-based implementation of Gower's distance.
#' 
#' @name gower-package
#' @docType package
#' @useDynLib gower, .registration=TRUE
#' 
{}


.onLoad <- function(libname, pkgname){
  max_threads <- 1L
  max_threads <- .Call("R_get_max_threads",PACKAGE="gower")
  thread_limit <- .Call("R_get_thread_limit",PACKAGE="gower")
  max_threads <- min(max_threads, thread_limit)
  # leave one core for the user to control the machine.
  if (max_threads > 2) max_threads <- max_threads - 1L
  options(gd_num_thread = as.integer(max_threads))
}
