



#' Generat range matrix for SVD
#'
#' @param L a pedigree's L inverse matrix in sparse 'spam' format
#' @param rank  `integer` how many principal components to return
#' @param depth `integer` number of iterations for generating the range matrix
#' @param numVectors `integer` > `rank` to specify the oversampling for the
#'
#' @return
#' @export
#'
#' @importFrom spam backsolve
#' @importFrom spam forwardsolve
#' @importFrom stats rnorm
#' @examples
randRangeFinder <- function(L, rank, depth, numVectors){
  dim <- nrow(L)
  testVectors <- rnorm(n = dim * numVectors)
  testMatrix <- matrix(testVectors, nrow = dim, ncol = numVectors)
  Q <- testMatrix
  for (i in 1:depth){
    qrObject <- base::qr(Q)
    Q <- qr.Q(qrObject)
    Q <- spam::backsolve(t(L),Q)
    Q <- spam::forwardsolve(L,Q)
  }
  qrObject <- qr(Q)
  Q <- qr.Q(qrObject)
  return(Q[,1:rank])
}

#' Singular value decomposition in sparse triangular matrix
#'
#' Uses random linear algebra, see Halko et al. (2010). Singular value
#' decomposition (SVD) decomposes a matrix \eqn{X=U\Sigma W^T}
#'
#' @param L a pedigree's L inverse matrix in sparse 'spam' format
#' @param rank  `integer` how many principal components to return
#' @param depth `integer` number of iterations for generating the range matrix
#' @param numVectors `integer` > `rank` to specify the oversampling for the
#'
#' @return A list of three: u (=U), d (=Sigma), and v (=W^T)
#' @export
#'
#' @importFrom spam backsolve
#' @examples
randSVD <- function(L, rank, depth, numVectors){
  # L: lower cholesky factor of animal matrix
  # rank: number of PCs
  # depth: power iteration, higher -> more accurate approximation, ~5 is usually sufficient
  # numVectors: usually rank + 5~10, must be larger than rank
  dim <- nrow(L)
  Q <- randRangeFinder(L, rank, depth, numVectors)
  C <- t(backsolve(t(L), Q))
  svdObject <- svd(C)
  U <- Q %*% svdObject$u
  D <- svdObject$d ** 2
  V <- svdObject$v
  return(list(u=U[,1:rank], d=D[1:rank], v=V[1:rank,]))
}



#' PCA of pedigree L inverse sparse matrix
#'
#' @param L a pedigree's L inverse matrix in sparse `spam` format
#' @param method `string` only randSVD (the default) is implemented
#' @param rank  `integer` how many principal components to return
#' @param depth `integer` number of iterations for generating the range matrix
#' @param numVectors `integer` > `rank` to specify the oversampling for the
#' range matrix
#'
#' @return A `list` of one element, `scores`, the principal components
#' @export
#'
#' @examples
#'
rppca <- function(L, method="randSVD", rank=10, depth=3, numVectors=15){
  #check L is the right kind of sparse matrix

    if(method=="randSVD"){
    rsvd = randSVD(L, rank=rank, depth=depth, numVectors=numVectors)
    scores = rsvd$u %*% diag(rsvd$d)
    dimnames(scores) <- list(NULL, paste0("PC", 1:rank))
    return(list(scores= scores,
                d=rsvd$d
                # Would be good to also return variance proportions.
                # PCA usually returns loadings but not sure this makes sense
                # here where we only have the loading of the range matrix?
                )
           )

  } else {
    stop(paste0("Method ", method," not implemented"))
  }
}

