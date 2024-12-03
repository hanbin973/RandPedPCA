# Example to show that PCs 1 and 2 are swapped for a small example that can be
# analysed in a naive way with R's built-in prcomp()


# Load data ---------------------------------------------------------------


# smaller simulation

linv <- importLinv("../datasets/pedLInv.mtx")
#linvMM <- Matrix::readMM("../datasets/pedLInv.mtx")
dim(linv)
pLab <- factor(read.table("../datasets/popLabel.csv", header = T)[,1])
class(linv)
dim(linv)
pc <- rppca(linv)
class(pc)
object.size(pc$x)
plot(pc$x[,1:2], col=pLab)
pc
summary(pc)

randTraceHutchinson(linv, numVectors = 40000)




# try  naive approach, gat a inv matrix
ainv <- readMM("../datasets/pedAInv.mtx")
ainv[1:10, 1:10]

# get a matrix
aa <- solve(ainv)
traa <- matrixcalc::matrix.trace(as.matrix(aa))


pcTr <- rppca(linv,totVar=traa)

# get a inv matrix also by multiplying L^I^T and L^I
ltl <- t(linv) %*% linv
# get a matrix
ltli <- solve(ltl)

matrixcalc::matrix.trace(as.matrix(ltli))

# compare both A matrices
all(zapsmall(ltli - aa , digits = 10)==0) # true
all(zapsmall(ltli - aa , digits = 13)==0) # false, identical down to 12 digits



# compute built-in PCA on A matrix
# pcaaNoscale <- prcomp(aa)
# pcaaNocenter <- prcomp(aa, center = F)
# summary(pcaaNocenter)$importance[,1:10]
# summary(pcaaNocenter)$importance[,1:4]
#
# pclltiNoscale <- prcomp(ltli)
#
# plot(sseq,sqrt(2)/(sseq-1))
#
# plot(pcaaNoscale$x[,1:2]) # 1 and 2 swapped
# plot(pcaaNoscale$x[,2:1]) # alright if PCs 1 and 2 swapped
#
# plot(pcaaNocenter$x[,1:2], col=pLab)
#
#
# plot(pcltliNoscale$x[,1:2]) # 1 and 2 swapped
# plot(pcltliNoscale$x[,2:1]) # alright if PCs 1 and 2 swapped

plot(pc$x[,1:2], col=pLab) # RandPed PCA
grid()
legend("topright",
       col=1:3,
       pch=1,
       legend=c("A", "AB", "B"))
# look into correlations between PCs computed by both methods
pcCopy <- pc$scores # a copy of the RandPedPCA scores (to rename PCs to rPCs)
dimnames(pcCopy)[[2]] <- paste0("r", dimnames(pcCopy)[[2]])

# compute correlation matrix
ccors <- as.matrix(cor(pcaaNoscale$x[,1:10], pcCopy))
ccors2 <- as.matrix(cor(pcaaNocenter$x[,1:10], pcCopy))
# visualise correlation matrix
ComplexHeatmap::Heatmap(as.matrix(ccors),
                        column_order = NULL,
                        row_order = NULL,
                        row_title = "prcomp (R built-in)",
                        column_title = "RandPedPCA", name = "Corr. coef.")

# this is looking better, centering turned off in prcomp()
ComplexHeatmap::Heatmap(as.matrix(ccors2),
                        column_order = NULL,
                        row_order = NULL,
                        row_title = "prcomp (R built-in)",
                        column_title = "RandPedPCA", name = "Corr. coef.")

#ll <- solve(linv) # does not work, go via Cholesky decomposition of A

# prcompL <- prcomp(linv)
# plot(prcompLinv$x[,1:2], col=pLab)
# crLinv <- solve(crossprod(linv))
# image(aa)
# image(linv)
# image(crLinv)
# image(aa - crLinv)
# all(zapsmall(aa - crLinv)==0)

chA <- chol(aa)
chA
image(chA)
#prcompCholA <- prcomp(chA, center = F)
prcompCholAT <- prcomp(t(chA), center=F, scale. = F)
# plot(prcompCholA$x[,1:2],col=pLab)
# grid()
# legend("topleft",
#        col=1:3,
#        pch=1,
#        legend=c("A", "AB", "B"))
#
#
# plot(prcompCholA$rotation[,1:2],col=pLab)
# summary(prcompCholA)$importance[,1:10]
summary(prcompCholAT)$importance[,1:10]
summary(prcompCholAT)$importance[,1:4]



# this is the right one!! (chA is the upper factor, must be transposed)
plot(prcompCholAT$x[,1:2],col=pLab)
grid()
legend("topright",
       col=1:3,
       pch=1,
       legend=c("A", "AB", "B"))
prcompCholAT$sdev[1:10]
prcompCholAT$sdev[1:10]^2
plot(prcompCholA$rotation[,1:2],col=pLab)
sqrt(pc$sdev * sqrt(2599)) / sqrt(2599)
zapsmall(summary(prcompCholA)$importance[1,1:10] - pc$sdev)
sum(prcompCholA$sdev)
sum(prcompCholAT$sdev)
# both 35.14031

sum(prcompCholAT$sdev^2)
# 1.32626

pc$sdev

prcompCholAT$sdev[1:10]
sqrt(2599)
# 50.98039
sum(diag(aa))
# 50.98039

hist(diag(aa))


svdchA <- svd(t(chA))
svdchA$d[1:10] / sqrt(2599)
sum(svdchA$d^2) - sum(diag(aa))


# What is the equivalent to X? --------------------------------------------


linv
dim(linv)
pcli <- prcomp(linv)

sum(pcli$sdev)
sum(pcli$sdev^2)

levels(pLab)
plot(pcli$x[,1:2], col = pLab)
plot(pcli$x[,2:3], col = pLab)
plot(pcli$rotation[,1:2], col = pLab)
dim(pcli$rotation)
dim(pcli$x)
summary(pcli)$importance[,1:10]


eaa <- eigen(aa)
length(eaa$values)
dim(eaa$vectors)
plot(eaa$vectors[,1:2], col =pLab)


# SVD on L
svdchA <- svd(t(chA))

svdchA$u[1:5, 1:5]
svdchA$v[1:5, 1:5]
svdchA$d[1:5]
# eigen decomp on A
eA <- eigen(aa)
eA$vectors[1:5, 1:5]
eA$values[1:5]
sqrt(eA$values[1:5])
all(zapsmall(svdchA$u - eA$vectors, digits = 6) == 0)

image(svdchA$u)
image(eA$vectors)

ccc <- cor(svdchA$u, eA$vectors)
image(ccc)





# More equivalences -------------------------------------------------------
# individuals in rows, trait in cols
# we are interested in relationships between individuals

# iris dataset, dropping the labels
i4 <- as.matrix(iris[,1:4])
i4cent <- scale(i4, scale = F, center = T) # traits centered
i4tcent <- scale(t(i4), scale = F, center = T) # inds centered
# sds of the data
sds <- apply(i4, 2, sd)
# Sepal.Length  Sepal.Width Petal.Length  Petal.Width
#    0.8280661    0.4358663    1.7652982    0.7622377
sum(sds)
# 3.791468
# vars of the data cols
vars <- apply(i4, 2, var)
# Sepal.Length  Sepal.Width Petal.Length  Petal.Width
#    0.6856935    0.1899794    3.1162779    0.5810063
sum(vars)
# 4.572957
pci <- prcomp(i4)
pciNN <- prcomp(i4, center = F, scale. = F)
pci$center
pci$scale
pciNN$center
pciNN$scale
prci <- princomp(i4)
summary(pci)
sum(pci$sdev)
sum(pci$sdev^2)
# identical to orig data
# 4.572957

plot(pci$x[,1:2], col=iris$Species, main="prcomp") # good, for reference


sum(prci$sdev^2)
# again, sum is basically the same
# 4.542471

#plot(prci$scores[,1:2], col=iris$Species)
plot(prci$scores[,1], -1* prci$scores[,2],
     col=iris$Species,
     main="princomp")                             # good too
prci$center
prci$scale


# doing it "by hand"
# we want to plot individuals and thus summarise traits we thus need a
# covariance matrix of the traits
# covariance matrices
cMat1 <- crossprod(i4tcent) # cov mat of the inds (large)
cMat2 <- crossprod(i4cent) # cov mat of the traits (small)


# covariance matrices as used by princomp()
cwt1 <- cov.wt(i4cent) # cov mat of traits (cols) thus small
cwt2 <- cov.wt(i4tcent) # cov mat of inds (rows) thus larger
dim(cwt1$cov)
dim(cwt2$cov)
cov.wt

cMat1 / cwt2$cov

cMat2 / cwt1$cov
# the above a scaled by 1/n-1


image(cMat1)
image(cwt2$cov)
e1 <- eigen(cMat1)
e1wt <- eigen(cwt2$cov)
dim(e1$vectors)

plot(e1$vectors[,1:2], col=iris$Species,
     main="Eigenvectors of individual cov mat") # strange


# tried to scale cov mat as is done inside princomp
plot(e1wt$vectors[,1], -1*e1wt$vectors[,2], col=iris$Species,
     main="not sure")


e2 <- eigen(cMat2)
e2wt <- eigen(cwt1$cov)
plot((i4cent %*% e2$vectors)[,1:2], col=iris$Species,
     main="Centered iris data %*% ED od trait cov mat") # scale ok
plot((i4cent %*% e2wt$vectors)[,1:2], col=iris$Species,
     main="Centered iris data %*% ED od trait cov mat") # scale ok



chol()

ewt1 <- eigen(cwt1$cov)
ewt2 <- eigen(cwt2$cov)
ewt2$values[1:10]
plot(ewt2$vectors[,1:2])
plot((i4cent %*% ewt1$vectors)[,1:2], col=iris$Species)
plot((t(i4cent) %*% ewt2$vectors)[,1:2])
plot((t(i4cent) %*% ewt2$vectors)[1:2,])
plot((ewt2$vectors %*% i4cent)[,1:2])
plot((ewt2$vectors %*% i4cent)[1:2,])
plot(ewt2$vectors[,1:2], col=iris$Species)

# SVD
isvd <- svd(i4cent)
isvd$d
sqrt(isvd$d)

plot(isvd$u[,1:2], col=iris$Species, main="SVD (u) of centered iris data") # scale wrong, but identical to ED of individual cov mat

isvd2 <- svd(cMat1)
plot(isvd2$u[,1:2], col=iris$Species, main="SVD (u) of centered iris individual cov mat") # scale wrong, but identical to ED of individual cov mat
plot(isvd2$v[,1:2])
#ch <- chol(cMat1) # does not work






# Hutchinson --------------------------------------------------------------

getNumVectorsHutchinson(.2, .95)

# Doubled pedigree (two unrelated pedigrees) ------------------------------


L10 <- importLinv("../datasets/LInv10.mtx")

L1010 <- importLinv("../datasets/LInv1010.mtx")


pc10 <- rppca(L10)
pc1010 <- rppca(L1010)



plot(pc10$scores[,1:2])
plot(pc1010$scores[,1:2])
plot(pc1010$scores[,c(1,3)])
plot(pc1010$scores[,c(1,4)])
plot(pc1010$scores[,c(1,5)])
plot(pc1010$scores[,c(1,6)])


# Comparing variance components -------------------------------------------




# variance components (built-in method)
summary(pcaaNoscale)$importance[,1:10]


# no great resemblance between these values of
#  SD from prcomp and d from random svd
summary(pcaaNocenter)$importance[,1:10]
#  PC1      PC2      PC3      PC4       PC5       PC6     PC7       PC8
#  Standard deviation     17.43141 14.79962 1.417906 1.330322 0.6352784 0.5078248 0.44014 0.3378062
#  Proportion of Variance  0.57397  0.41374 0.003800 0.003340 0.0007600 0.0004900 0.00037 0.0002200
#  Cumulative Proportion   0.57397  0.98770 0.991500 0.994840 0.9956100 0.9960900 0.99646 0.9966800
#  PC9      PC10
#  Standard deviation     0.2848558 0.2750638
#  Proportion of Variance 0.0001500 0.0001400
#  Cumulative Proportion  0.9968300 0.9969700
pc$sdev
# [1] 17.4314082 14.7996244  1.4178247  1.3302139  0.6318616  0.4974759  0.4243559  0.3244110  0.2563754  0.2230309


sdSum <- sum(summary(pcaaNocenter)$importance[1,])
# 67.61325
varSum <- sum(summary(pcaaNocenter)$importance[1,]^2)
# 529.3925

summary(pcaaNocenter)$importance[1,1:10]^2/varSum
# as prcomp

# same again
round(pc$sdev^2/varSum, 3)

sum(diag(aa))
# 3446.95
sqrt(sum(diag(aa)))
# 58.71073

sum(diag(aa)) / sqrt(dim(aa)[1]-1) # same as sum of SDs from prcomp
# 67.61325



summary(pcaaNocenter)$importance[1,1:10]^2/ sum(diag(aa)) / sqrt(dim(aa)[1]-1)
# far too low, is the trace not what we need?
sum(diag(aa))
# 3446.95

3446.95 / sqrt(2599)
# same as total sd from prcomp of A!
# [1] 67.61326



round(pc$sdev^2 / (sum(diag(aa)) / sqrt(2599))^2, 3)
round(pc$sdev / (sum(diag(aa)) / sqrt(2599)), 3)
pc$sdev^2 / sum(diag(aa)^2) / 2599

aaColVar <- apply(aa, 2, var)
plot(aaColVar)
sum(aaColVar)
image(aa)





# Comparing pedigree matrices ---------------------------------------------

# Matrix::writeMM(obj = pedA,    file = "pedA.mtx")
# Matrix::writeMM(obj = pedAInv, file = "pedAInv.mtx")
# Matrix::writeMM(obj = pedLInv, file = "pedLInv.mtx")



AAi <- Matrix::readMM("../datasets/pedAInv.mtx")
AA <- solve(AAi)
Li <- Matrix::readMM("../datasets/pedLInv.mtx")
image(AAi)
image(AA)
image(Li)
sum(diag(AA))
zapsmall(AA[1000:1010, 1000:1010])



# Computing the total variance from a pedigree ----------------------------
library(pedigreeTools)

metaData <- read.table("../datasets/pedMeta.csv",
                       header=T,
                       sep=",")
head(metaData)

ped <- pedigree(sire  = metaData$fid,
                dam   = metaData$mid,
                label = metaData$id)

ped
inbr <- inbreeding(ped)
pedA <- getA(ped)
plot( inbr, diag(as.matrix(pedA)))
abline(1,1, col=2)
sum(inbr+1)


pc2 <- rppca(ped)
plot(pc2$x[,1:2], col=factor(metaData$population))
pc2$sdev^2

pc2$sdev
pc$sdev
pc2$varProps
sum(diag(aa))
# 3549.629
svdChA <- svd(chA)
sum(svdChA$d)
# 1797.833
sum(svdChA$d^2)
# 3549.629 (as diag sum of A)

(pc2$sdev * sqrt(2649))^2 / sum(svdChA$d^2)
pc2$sdev^2 * 2649 / sum(svdChA$d^2)

summary(prcompCholAT)$importance[,1:10]

pcaaNocenter$sdev[1:10]
#summary(prcompCholAT)$importance[,1:10] # dont want that!, wrong input (A insted of L)
sum(prcompCholAT$sdev^2) * sqrt(2599)


# Matrix conversions ------------------------------------------------------

llii <- getLInv(ped) # dtCMatrix
lliiCsparse <- as(llii, "CsparseMatrix")
class(llii)
linvMM # readMM returns dgTMatrix
class(linvMM)
class(lliiCsparse)



as(llii, "CsparseMatrix") # returns dtCMatrix
as(linvMM, "CsparseMatrix") # returns dgCMatrix
as.spam.dgCMatrix(as(llii, "dgCMatrix"))
as.spam.dgCMatrix(as(llii, "CsparseMatrix"))
as.spam.dgCMatrix(as(linvMM, "dgCMatrix"))
as.spam.dgCMatrix(as(linvMM, "CsparseMatrix"))
?as.spam.dgCMatrix
?as
showMethods(coerce)
