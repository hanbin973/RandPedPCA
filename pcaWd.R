library(spam)
library(Matrix)
setwd("~/git_repos/RandPedigreePCA/")
randRangeFinder <- function(L, rank, depth, numVectors){
  dim <- nrow(L)
  testVectors <- rnorm(n = dim * numVectors)
  testMatrix <- matrix(testVectors, nrow = dim, ncol = numVectors)
  Q <- testMatrix
  print("Test matrix done")
  for (i in 1:depth){
    print("Iteration")
    print(i)
    qrObject <- Matrix::qr(Q)
    
    Q <- qr.Q(qrObject)
    print("QR decomp. done")
    Q <- spam::backsolve(t(L),Q) 
    print("backsolve done")
    Q <- spam::forwardsolve(L,Q)
    print("forwardsolve done")
    #print(summary(L))
    #Q <- Matrix::solve(t(L), Q)#, system="LDLt")
    
    #Q <- backsolve(L, forwardsolve(L, Q, transpose=T, upper.tri=T))
    print("Solve worked.")
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
  print("Going to run range finder")
  Q <- randRangeFinder(L, rank, depth, numVectors)
  print("Range finder done")
  C <- t(spam::backsolve(t(L), Q))
  svdObject <- svd(C)
  U <- Q %*% svdObject$u
  D <- svdObject$d ** 2
  V <- svdObject$v
  return(list(u=U[,1:rank], d=D[1:rank], v=V[1:rank,]))
}


t0 <- Sys.time()
t0
#L <- spam::read.MM('~/temp/pedLInv.mtx')
Sys.time() - t0

#L <- Matrix::readMM('datasets/pedLInv.mtx')

L <- Matrix::readMM('~/temp/pedLInv.mtx')


# does not work
#K <- as.spam(L)
svdObject <- randSVD(L, 5, 10, 15)



typeof(L)
as(L, "CsparseMatrix")
K <- plot(svdObject$u[,1:2])

# for sanity check
#invAdense <- solve(as.matrix(A))
#svdObjectDense <- svd(invAdense)


# other stuff -------------------------------------------------------------



a <- matrix(rnorm(12), 3, 4)
A <- a %*% t(a)
chl <- Cholesky(A)
typeof(chl)
str(chl)

Matrix::solve(chl, c(1,2,3))


Matrix()


?Cholesky
