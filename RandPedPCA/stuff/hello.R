# Example to show that PCs 1 and 2 are swapped for a small example that can be
# analysed in a naive way with R's built-in prcomp()

# smaller simulation

linv <- as.spam.dgCMatrix(as(readMM("../datasets/pedLInv.mtx"), "CsparseMatrix"))
pc <- rppca(linv)
plot(pc$scores[,1:2])

# try  naive approach, gat a inv matrix
ainv <- readMM("../datasets/pedAInv.mtx")
ainv[1:10, 1:10]

# get a matrix
aa <- solve(ainv)

# get a inv matrix also by multiplying L^I^T and L^I
ltl <- t(linv) %*% linv
# get a matrix
ltli <- solve(llt)

# compare both A matrices
all(zapsmall(ltli - aa , digits = 10)==0) # true
all(zapsmall(ltli - aa , digits = 13)==0) # false, identical down to 12 digits



# compute built-in PCA on A matrix
pcaaNoscale <- prcomp(aa)
pcaaNocenter <- prcomp(aa, center = F)
pclltiNoscale <- prcomp(ltli)



plot(pcaaNoscale$x[,1:2]) # 1 and 2 swapped
plot(pcaaNoscale$x[,2:1]) # alright if PCs 1 and 2 swapped

plot(pcltliNoscale$x[,1:2]) # 1 and 2 swapped
plot(pcltliNoscale$x[,2:1]) # alright if PCs 1 and 2 swapped

plot(pc$scores[,1:2]) # RandPed PCA

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
pc$d
#  [1] 888.65996 754.49060  72.27821  67.81399  32.01035  25.14347  21.42995  14.56372  12.50480
#  [10]  12.15989
sqrt(pc$d)
#  [1] 29.810400 27.467992  8.501659  8.234925  5.657769  5.014326  4.629250  3.816244  3.536213
#  [10]  3.487103
