library(Matrix)
setwd("~/git_repos/RandPedigreePCA/")
randRangeFinder <- function(L, rank, depth, numVectors){
  dim <- nrow(L)
  testVectors <- rnorm(n = dim * numVectors)
  testMatrix <- matrix(testVectors, nrow = dim, ncol = numVectors)
  Q <- testMatrix
  for (i in 1:depth){
    qrObject <- Matrix::qr(Q)
    Q <- qr.Q(qrObject)
    print("fine up to here")
    print(typeof(L))
    print(typeof(t(L)))
    Q <- Matrix::solve(t(L),Q) 
    Q <- Matrix::solve(L,Q)
  }
  qrObject <- qr(Q)
  Q <- qr.Q(qrObject)
  return(Q[,1:rank])
}

randSVD <- function(L, rank, depth, numVectors){
  # L: lower cholesky factor of animal matrix
  # rank: number of PCs
  # depth: power iteration, higher -> more accurate approximation, ~5 is usually sufficient
  # numVectors: usually rank + 5~10, must be larger than rank
  dim <- nrow(L)
  Q <- randRangeFinder(L, rank, depth, numVectors)
  print("Range finder done")
  C <- t(Matrix::solve(t(L), Q))
  print("Last backsolve done")
  svdObject <- svd(C)
  print("svd done")
  U <- Q %*% svdObject$u
  D <- svdObject$d ** 2
  V <- svdObject$v
  return(list(u=U[,1:rank], d=D[1:rank], v=V[1:rank,]))
}

#L <- spam::read.MM('~/temp/pedLInv.mtx')
L <- Matrix::readMM('~/temp/pedLInv.mtx')
#K <- as(L, "CsparseMatrix")


#L <- Matrix::readMM('datasets/pedLInv50.mtx')
t0 <- Sys.time()
svdObject <- randSVD(L, 10, 2, 15)
Sys.time() -t0
# 33s on hux
# 0.1s for ped50

smoothScatter(svdObject$u[,1:2])

# for sanity check
#invAdense <- solve(as.matrix(A))
#svdObjectDense <- svd(invAdense)

