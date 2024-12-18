

# get data
i4 <- as.matrix(iris[,1:4])
i4cent <- scale(i4, scale = F, center = T) # traits centered
i4tcent <- scale(t(i4), scale = F, center = T) # inds centered



# Prcomp ------------------------------------------------------------------

# prcompom on centered i

pci4cent <- prcomp(i4cent, center = F)
# centered by prcomp
pci4 <- prcomp(i4cent, center = F)

# identical as expected
plot(pci4cent$x[,1:2], col=iris$Species)
plot(pci4$x[,1:2], col=iris$Species)



# Cov matrices ------------------------------------------------------------

cprTrait <- crossprod(i4)
cprTraitCent <- crossprod(i4cent)
cprIndTraitCent <- tcrossprod(i4cent)
cprInd <- tcrossprod(i4)


# 1st and 2nd are 4x4 cov mats
svdList <- lapply(list(cprTrait, cprTraitCent, cprIndTraitCent, cprInd), svd)

plot((i4 %*% svdList[[1]]$v)[,1:2], col=iris$Species)
plot((i4cent %*% svdList[[2]]$v)[,1:2], col=iris$Species)

plot((svdList[[3]]$u %*% diag(svdList[[3]]$d))[,1:2], col=iris$Species)
plot((svdList[[4]]$u %*% diag(svdList[[4]]$d))[,1:2], col=iris$Species)
