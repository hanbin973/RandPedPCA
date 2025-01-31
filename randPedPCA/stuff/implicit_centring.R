# code to demonstrate the expected outcome of centring L^(-1)

# generate pedigree from the metadata of the 2nd test dataset (two diverged populations)
ped2 <- pedigree(pedMeta2$fid,
                 pedMeta2$mid,
                 pedMeta2$id)
# get L
L2 <- getL(ped2)

# run naive PCA
pcL2 <- prcomp(t(L2), center = F, scale.=F)
# run naive PCA, this time centring L (which is equivalent to centring A)
pcL2c <- prcomp(t(L2), center = T, scale.=F)

# now rppca w/o centring
pcLI2 <- rppca(pedLInv2)
# rppca with centring
pcLI2c <- rppca(pedLInv2, center=T)

pcPed2 <- rppca(ped2)
pcPed2c <- rppca(ped2, center=T)

plot(pcL2$x[,1:2], main="Naive, not centred")
plot(pcL2c$x[,1:2], main="Naive, centred")
plot(pcLI2$x[,1:2], main="rppca, not centred") # same as naive, not centred as expected
plot(pcLI2c$x[,1:2], main="rppca, centred") # same pattern as centred naive as expected


plot(pcL2c$x[,1], pcLI2c$x[,1])
plot(pcL2c$x[,2], pcLI2c$x[,2])

summary(lm(pcL2c$x[,1]~pcLI2c$x[,1]))
summary(lm(pcL2c$x[,2]~pcLI2c$x[,2]))
plot(pcL2c$x[,1], pcPed2c$x[,1])
plot(pcL2c$x[,2], pcPed2c$x[,2])

plot(pcLI2c$x[,1], pcPed2c$x[,1])
plot(pcLI2c$x[,2], pcPed2c$x[,2])

pcLI2$sdev
pcLI2$sdev^2

plot(pcL2c$x[,1], pcLI2c$x[,1]/sqrt(2650))
abline(0,1)


plot(pcL2c$x[,2], pcLI2c$x[,2]/sqrt(2650))
abline(0,1)


summary(pcPed2c)
summary(pcPed2)

summary(pcLI2)
summary(pcLI2c)
