



#' Generat range matrix for SVD
#'
#' @param L a pedigree's L inverse matrix in sparse 'spam' format
#' @param rank  `integer` how many principal components to return
#' @param depth `integer` number of iterations for generating the range matrix
#' @param numVectors `integer` > `rank` to specify the oversampling for the
#'
#' @return The range matrix for `randSVD`
#' @export
#'
#' @importFrom spam backsolve
#' @importFrom spam forwardsolve
#' @importFrom stats rnorm
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

#' Compute the number of vectors used for trace estimation
#'
#' Follows Skorski, M. (2021). Modern Analysis of Hutchinson’s Trace Estimator.
#' 2021 55th Annual Conference on Information Sciences and Systems (CISS), 1–5.
#' https://doi.org/10.1109/CISS50987.2021.9400306
#'
#' @param e `numeric` denoting the relative error margin
#' @param d `numeric`. 1-d is the probability the the relative error is bounded
#' by e.
#'
#' @return a scalar
#' @export
getNumVectorsHutchinson <- function(e, d){
  2*(2+8*sqrt(2)/3*e)*log(2/d)/e^2
}


#' Trace estimation for sparse L inverse matrices
#'
#' Using Hutchinson's method
#'
#' @param L A pedigree's L inverse matrix
#' @param numVectors, an `integer` specifying how many random vectors to use
#'
#' The higher `numVectors`, the higher the accuracy and the longer the runtime.
#' Accuracy can be estimated with the function `getNumVectorsHutchinson`.
#'
#' @return a scalar
#' @export
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

#' PCA of pedigree L inverse sparse matrix
#'
#'
#'
#' @param L a pedigree's L inverse matrix in sparse `spam` format
#' @param method `string` only randSVD (the default) is implemented
#' @param rank  `integer` how many principal components to return
#' @param depth `integer` number of iterations for generating the range matrix
#' @param numVectors `integer` > `rank` to specify the oversampling for the
#' range matrix
#' @param nVecTraceEst `integer` number of random vectors to be used for
#' relationship matrix trace estimation
#' @param returnRotation `logical`whether or not to returen the rotation
#' matrix, i.e. the right singulare values of the relationship matrix. FALSE by
#' default.
#'
#' The output slots are named like those if R's built in `prcomp` function.
#' Rotation is not returned by default as it is the transpose of the PC scores,
#' which are returned in `x`. `scale` and `center` are set to `FALSE`.
#'
#' @return A `list` containing:
#' * `x`, the principal components,
#' * `sdev`, the variance components of each PC. Note that the total variance is
#' not known per se and this these components cannot be used to compute the
#' proportion of the total variance accounted for by each PC. However, if
#' `nVecTraceEst` is specified, `rppca` will estimate the total variance and
#' return variance proportions.
#' * `vProp` the estimated variance proportions accounted for by each PC.
#' Only returned if `nVecTraceEst` is set.
#' * `scale` always `FALSE`
#' * `center` always `FALSE`
#' * `rotation` the right singular values of the (implicit) relationship matrix.
#' Only returned if `returnRotation == TRUE`
#'
#' @export
#'
#' @examples rppca(pedLinv)
#'
rppca <- function(L,
                  method="randSVD",
                  rank=10,
                  depth=3,
                  numVectors=15,
                  nVecTraceEst,
                  returnRotation=FALSE){
  #check L is the right kind of sparse matrix
  nn <- dim(L)[1]
  if(method=="randSVD"){
    rsvd = randSVD(L, rank=rank, depth=depth, numVectors=numVectors)
    scores = rsvd$u %*% diag(rsvd$d)
    dimnames(scores) <- list(NULL, paste0("PC", 1:rank))
    if(returnRotation){
      pc <- list(x= scores,
                 sdev=rsvd$d / sqrt(max(1, nn-1)),
                 center=FALSE,
                 scale=FALSE,
                 rotation=t(x)
      )
    } else {
      pc <- list(x= scores,
                 sdev=rsvd$d / sqrt(max(1, nn-1)),
                 center=FALSE,
                 scale=FALSE
      )
    }

    class(pc) <- "rppca"
    return(pc)

  } else {
    stop(paste0("Method ", method," not implemented"))
  }
}

