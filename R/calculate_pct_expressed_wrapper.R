#' Calculate Percent Expression
#'
#' @param data A numeric matrix of gene expression data.
#' @param cellTypes A character vector of cell types.
#' @param genes A character vector of gene names.
#'
#' @return A numeric matrix of percent expression values.
#' @import Rcpp
#' @useDynLib resrect, .registration = TRUE
#' @export
calculatePctExp <- function(data, cellTypes, genes) {
  .Call('_resrect_calculatePctExp', PACKAGE = 'resrect', data, cellTypes, genes)
}

# Register the C++ function with Rcpp
Rcpp::sourceCpp('./../RESRECT/src/calculate_pct_exp.cpp')
