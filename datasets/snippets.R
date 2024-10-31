library(pedigreeTools)

pp <- pedigree(c(NA, NA, 1, 1, 3, 4, 6, 5, 6, 5),
               c(NA, NA, 2, 2, 4, 3, 5, 6, 5, 6),
               1:10)
pp
getA(pp)
getL(pp)
getLInv(pp)

# empty matrix of passed on ancestry
pMatInv <- matrix(0.0, nrow = 10, ncol = 10)

# fill matrix
for(i in 1:nrow(ped)){
  pMatInv[i,i] <- 1
  if(!is.na(ped[i,2])){
    pMatInv[i, ped[i,2]] <- -1/2
  }
  if(!is.na(ped[i,3])){
    pMatInv[i, ped[i,3]] <- -1/2
  }
}
pMatInv
solve(pMatInv)


i4 <- as.matrix(iris[,1:4])

i4s <- apply(i4,1,scale, scale = F)
typeof(i4s)
prd <- t(i4s) %*% i4s
cvmat <- cov(i4s)
image(cvmat)
image(prd)


