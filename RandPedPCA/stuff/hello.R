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

# visualise correlation matrix
ComplexHeatmap::Heatmap(as.matrix(ccors),
                        column_order = NULL,
                        row_order = NULL,
                        row_title = "prcomp (R built-in)",
                        column_title = "RandPedPCA", name = "Corr. coef.")

# variance components (built-in method)
summary(pcaaNoscale)$importance[,1:10]



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
