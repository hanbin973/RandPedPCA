library(randPedPCA)

ped <- pedigree(sire=pedMeta$fid,
                dam=pedMeta$mid,
                label=pedMeta$id)
li <- getLInv(ped)
ll <- getL(ped)

image(ll)
image(li)
show(ll)
show(li)

# returns the value of A * G (but taking L^-1 as input)
oraculum <- function(Li, G){
  Y <- spam::backsolve(t(Li), G)
  spam::forwardsolve(Li, Y)
}


hutchinson <- function(matVecOracle, num_queries, dimension = -1, hutch_dist = function(m, n) 2*matrix(sample(1:2, size = m*n, replace = TRUE, prob = c(0.5, 0.5)), nrow = m, ncol = n) - 3){
  if(!is.function(matVecOracle)){
    dimension <- nrow(matVecOracle)
    matVecOracle2 <- function(z) matVecOracle %*% z
  } else if(dimension == -1){
    stop('Passed MatVec handle but not the matrix dimension. "dimension" is the third (optional) input argument.')
  }
  G <- hutch_dist(dimension, num_queries)
  trace_est <- sum(diag(t(G) %*% matVecOracle2(G))) / num_queries
  return(trace_est)
}

aa <- getA(ped)

tv <- sum(diag(aa))

t0 <- Sys.time()
hutchinson(aa, 100)
Sys.time() - t0




hutchplusplus <- function(matVecOracle, num_queries, dimension = -1, hutch_dist = function(m, n) 2*matrix(sample(1:2, size = m*n, replace = TRUE, prob = c(0.5, 0.5)), nrow = m, ncol = n) - 3, sketch_dist = function(m, n) 2*matrix(sample(1:2, size = m*n, replace = TRUE, prob = c(0.5, 0.5)), nrow = m, ncol = n) - 3, sketch_frac = 2/3){
  if(!is.function(matVecOracle)){
    dimension <- nrow(matVecOracle)
    matVecOracle2 <- function(z) matVecOracle %*% z
  } else if(dimension == -1){
    stop('Passed MatVec handle but not the matrix dimension. "dimension" is the third (optional) input argument.')
  }
  S_num_queries <- round(num_queries * sketch_frac / 2)
  Hutch_num_queries <- num_queries - S_num_queries
  S <- sketch_dist(dimension, S_num_queries)
  Q <- qr.Q(qr(matVecOracle2(S), 0))
  G <- hutch_dist(dimension, Hutch_num_queries)
  G <- G - Q %*% t(Q) %*% G
  trace_est <- sum(diag(t(Q) %*% matVecOracle2(Q))) + sum(diag(t(G) %*% matVecOracle2(G))) / Hutch_num_queries
  return(trace_est)
}

t0 <- Sys.time()
hutchplusplus(aa, 200)
Sys.time() - t0
h <- sapply(1:100, function(x) hutchinson(aa, 10))
hpp <- sapply(1:100, function(x) hutchplusplus(aa, 10))

hist(h)
hist(hpp, col =2, add=T)


hutchinsonLi <- function(Li,
                         num_queries,
                         hutch_dist = function(m, n) 2*matrix(sample(1:2, size = m*n, replace = TRUE, prob = c(0.5, 0.5)), nrow = m, ncol = n) - 3
                         ){
  dimension <- dim(Li)[1]
  G <- hutch_dist(dimension, num_queries)
  trace_est <- sum(diag(t(G) %*% oraculum(Li, G))) / num_queries
  return(trace_est)
}
hutch_dist  = function(m, n) 2*matrix(sample(1:2, size = m*n, replace = TRUE, prob = c(0.5, 0.5)), nrow = m, ncol = n) - 3
sketch_dist = function(m, n) 2*matrix(sample(1:2, size = m*n, replace = TRUE, prob = c(0.5, 0.5)), nrow = m, ncol = n) - 3
hutchplusplusLi <- function(Li,
                            num_queries,
                            sketch_frac = 2/3){
  dimension <- dim(Li)[1]
  print("Dim done.")
  S_num_queries <- round(num_queries * sketch_frac / 2)

  Hutch_num_queries <- num_queries - S_num_queries
  print("Num queries done.")
  S <- sketch_dist(dimension, S_num_queries)
  print("S done.")
  Q <- qr.Q(Matrix::qr(oraculum(Li, S), 0))
  print(class(Q))
  print(dim(Q))
  print("Q done.")
  G <- hutch_dist(dimension, Hutch_num_queries)
  print("G done.")
  print(dim(G))
  G <- G - Q %*% (t(Q) %*% G)
  print("G done 2.")
  trace_est <- sum(diag(t(Q) %*% oraculum(Li, Q))) + sum(diag(t(G) %*% oraculum(Li, G))) / Hutch_num_queries
  return(trace_est)
}
class(dogli)
?qr
t0 <- Sys.time()
hutchplusplus(aa, 100)
Sys.time() - t0

t0 <- Sys.time()
hutchplusplusLi(pedLInv, 100)
Sys.time() - t0
traceback()




image(dogli)
dogli <- importLinv("~/Dropbox/Code/PedPCA/pedLInv.mtx")
t0 <- Sys.time()
#h01 <- hutchinsonLi(dogli, 10) # OK
h01 <- hutchplusplusLi(dogli, 10) # computer says no
Sys.time() - t0
traceback()


S <- matrix(rnorm(dim(dogli)[1] * 10), nrow=dim(dogli)[1])
# check return value of oraculum
or <- oraculum(dogli, S)
class(or)
class(dogli)

dogliEst <- sapply(1:100, function(x) hutchinsonLi(dogli, 10))
dogliEstPP <- sapply(1:100, function(x) hutchplusplusLi(dogli, 10))


hist(dogliEstPP, main="Trace of A (1.3M Labrador Retrievers)\n100 replicate estimates",
     col="#FF000040", breaks=seq(1350000, 1650000, by=10000))
hist(dogliEst, col="#40404040", breaks=seq(1350000, 1650000, by=10000),
     add=T)
legend("topright",
       fill=c("#40404040","#FF000040"),
       legend=c("Hutchinson", "Hutch++")
       )
# Checking types ----------------------------------------------------------

S <- sketch_dist(dim(pedLInv)[1], 100)
class(S)

class(oraculum(pedLInv, S))

class(spam::backsolve(t(pedLInv), S))
image(spam::backsolve(t(pedLInv), S))

Q <- qr.Q(Matrix::qr(oraculum(Li, S), 0))


# Matrix multiplication dimension -----------------------------------------


B <- matrix(1:9, nrow=1)
Bt <- t(B)
B %*% Bt
Bt %*% B
