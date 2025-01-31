


ped1 <- pedigree(pedMeta$fid,
                 pedMeta$mid,
                 pedMeta$id)
AA <-  getA(ped1)
LI1 <- getLInv(ped1)
L1 <- getL(ped1)
image(sparse2spam(LI1))
image(pedLInv)



svdAA <- svd(AA)

svdLI <- randSVD(pedLInv, 10, 3, 15)
svdAA$d[1:10]
svdLI$d[1:10]

svdAA$u[100:110,1]
svdLI$u[100:110,1]


pcL1 <- prcomp(t(L1), center = F)
pcL1c <- prcomp(t(L1), center = T)
pcLI <- rppca(pedLInv)
pcLIc <- rppca(pedLInv, cent=T)

pcL1$sdev[1:10]
pcL1c$sdev[1:10]
pcLI$sdev[1:10]
pcLIc$sdev[1:10]

plot(pcL1$x[,1:2])
plot(pcL1c$x[,1:2], main="Pedigree, with centered L")
plot(pcLI$x[,1:2], main="Pedigree PCA (default)")
plot(pcLIc$x[,1:2], main="Pedigree PCA (centred)")

pcGcent <- prcomp(pedGeno, center=T)
pcGnocent <- prcomp(pedGeno, center=F)
pcGcent$sdev[1:10]
plot(pcGnocent$x[,1:2], main="SNP pCA, not centered")
plot(pcGcent$x[,1:2], main="SNP PCA, centered")

plot(pcGcent$x[,1], pcL1c$x[,1], main='Correlation PCs1')
plot(pcGcent$x[,2], pcL1c$x[,2], main='Correlation PCs2')



# 2nd scenario - diverged founders ----------------------------------------

ped2 <- pedigree(pedMeta2$fid,
                 pedMeta2$mid,
                 pedMeta2$id)
AA2 <-  getA(ped2)
L2 <- getL(ped2)

pcL2 <- prcomp(t(L2), center = F)
pcL2c <- prcomp(t(L2), center = T)
pcLI2 <- rppca(pedLInv2)
pcLI2c <- rppca(pedLInv2, cent=T)


#pcG2cent <- prcomp(pedGeno2, center=T) # not included in the package anymore as too large
#pcG2nocent <- prcomp(pedGeno2, center=F) # not included in the package anymore as too large

# plotting
plot(pcL2c$x[,1:2], main="Pedigree2, with centered L")
plot(pcL2$x[,1:2], main="Pedigree2, with un-centered L")
plot(pcLI2$x[,1:2], main="Pedigree2 PCA (default)")
plot(pcLI2c$x[,1:2], main="Pedigree2 PCA (centred)")

plot(pcG2nocent$x[,1:2], main="SNP2 PCA, not centered")
plot(pcG2cent$x[,1:2], main="SNP2 PCA, centered")

plot(pcL2c$x[,1], pcG2cent$x[,1])
plot(pcL2c$x[,2], pcG2cent$x[,2])
