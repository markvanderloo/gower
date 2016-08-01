#' Gower's distance/similarity measure.
#'
#' A C-based implementation of Gower's distance.
#' 
#' @name gower-package
#' @docType package
#' @useDynLib gower R_gower R_gower_topn
#' 
{}


.onLoad <- function(libname, pkgname){
  options(gd_num_thread = parallel::detectCores())
}