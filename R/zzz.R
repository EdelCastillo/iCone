#' iCone
#' 
#'  MSI data peak picking
#' 
#' @name rPPGAS
#' @aliases iCone-package 
#' @author Esteban del Castillo
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib iCone
NULL  

.onUnload <- function (libpath)
{
  library.dynam.unload("iCone", libpath)
}