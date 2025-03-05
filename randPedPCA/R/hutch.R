
#' Compute the number of vectors to use for Hutchinson trace estimation
#'
#' Follows Skorski, M. (2021). Modern Analysis of Hutchinson’s Trace Estimator.
#' 2021 55th Annual Conference on Information Sciences and Systems (CISS), 1–5.
#' https://doi.org/10.1109/CISS50987.2021.9400306
#'
#' @param e A \code{numeric} denoting the relative error margin
#' @param d A \code{numeric}. 1-d is the probability the the relative error is bounded
#' by e.
#'
#' @return a scalar
getNumVectorsHutchinson <- function(e, d){
  2*(2+8*sqrt(2)/3*e)*log(2/d)/e^2
}


#' Trace estimation for sparse L inverse matrices
#'
#' Using Hutchinson's method
#'
#' @param L A pedigree's L inverse matrix
#' @param numVectors, an \code{integer} specifying how many random vectors to use
#'
#' If you do not have a good reason to do otherwise, use the function \code{hutchpp} instead.
#'
#' The higher \code{numVectors}, the higher the accuracy and the longer the running time.
#' Accuracy can be estimated with the function \code{getNumVectorsHutchinson}.
#'
#' @return a scalar
randTraceHutchinson <- function(L, numVectors){
  dim <- nrow(L)
  testVectors <- rnorm(n = dim * numVectors)
  testMatrix <- matrix(testVectors, nrow = dim, ncol = numVectors)
  testMatrixColNorms <- apply(testMatrix, 2, function(col) sqrt(sum(col^2)))
  normTestMatrix <- sweep(testMatrix, 2, testMatrixColNorms, FUN="/") * sqrt(dim)

  Y <- spam::backsolve(t(L), normTestMatrix)
  Y <- spam::forwardsolve(L, Y)
  Ests <- colSums(Y * normTestMatrix)
  return(mean(Ests))
}
