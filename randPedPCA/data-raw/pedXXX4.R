# ---- DESCRIPTION --------------------------------------------------------

# Simulation of three populations, of which two (A and B) are separate and
# the third (AB) is crossbred/admixed with continuous migration (from A and B).
# Populations A and B are selected on two different traits, while population AB
# is selected on an index of these two traits (just to build-up differences).

# The simulation can use one common founder population or two founder
# populations that have diverged some time ago.

# The simulation saves:
# * pedigree information (=expected Identity By Descent (IBD), eIBD)
# * IBD haplotypes since the founders (=realised IBD, rIBD)
# * Identity By State (IBS) genotypes
# * Trait genetic and phenotypic values

# We would only really need eIBD, but we save more for comparison and pedagogy.

# ---- INSTALLATION ------------------------------------------------------

#pkgs <- c("AlphaSimR", "dplyr", "pedigreeTools", "SIMplyBee")
#install.packages(pkg = pkgs)

# ---- SETUP -------------------------------------------------------------

library(package = "AlphaSimR")
library(package = "dplyr")
library(package = "pedigreeTools")

# ---- SIMULATION - ONE OR TWO FOUNDER POPULATION ------------------------
set.seed(987765543)
# Simulate founder genomes - one common founder population
founderGenomes <- runMacs(nInd = 200, nChr = 10, segSites = 1100,
                          species = "GENERIC")

# Simulate founder genomes - two founder populations that have a common origin
#founderGenomes2 <- runMacs(nInd = 100, nChr = 10, segSites = 1100, split = 100,
#                           species = "GENERIC")
# ... uncomment the line below and run the code to use the two-founder population
# founderGenomes <- founderGenomes2

# Initiate simulation parameters
SP <- SimParam$new(founderPop = founderGenomes)
# ... track global pedigree and recombinations
# (recobinations will enable tracking IBD)
SP$setTrackPed(isTrackPed = TRUE)
SP$setTrackRec(isTrackRec = TRUE)
# ... add two complex traits
varG <- matrix(data = c( 1.0, -0.3,
                        -0.3,  1.0), byrow = TRUE, nrow = 2)
SP$addTraitA(nQtlPerChr = 100, mean = c(0, 0), var = diag(varG), corA = varG)
varE <- matrix(data = c(2.0, 0.0,
                        0.0, 2.0), byrow = TRUE, nrow = 2)

# Monitor function
collectData <- function(pop, data = NULL, population, generation) {
  remove <- FALSE
  if (is.null(data)) {
    remove <- TRUE
    data <- vector(mode = "list", length = 3)
    names(data) <- c("pedigree", "haploIBD", "genoIBS")
    data$pedigree <- data.frame(id = NA, population = NA, generation = NA,
                                mid = NA, fid = NA,
                                gv1 = NA, pv1 = NA,
                                gv2 = NA, pv2 = NA)
    data$haploIBD <- matrix(data = NA, ncol = sum(pop@nLoci))
    data$genoIBS <- matrix(data = NA, ncol = sum(pop@nLoci))
  }
  data$pedigree <- rbind(data$pedigree,
                         data.frame(id = pop@id,
                                    population = population,
                                    generation = generation,
                                    mid = pop@mother,
                                    fid = pop@father,
                                    gv1 = pop@gv[, 1],
                                    pv1 = pop@pheno[, 1],
                                    gv2 = pop@gv[, 2],
                                    pv2 = pop@pheno[, 2]))
  data$haploIBD <- rbind(data$haploIBD,
                         pullIbdHaplo(pop = pop))
  data$genoIBS <- rbind(data$genoIBS,
                        pullSegSiteGeno(pop = pop))
  if (remove) {
    data$pedigree <- data$pedigree[-1, ]
    data$haploIBD <- data$haploIBD[-1, ]
    data$genoIBS <- data$genoIBS[-1, ]
  }
  return(data)
}

# Founder population & split
founders <- newPop(rawPop = founderGenomes)
founders <- setPheno(pop = founders, varE = diag(varE))
#popA <- founders[1:50]
#popB <- founders[51:100]
#popC <- founders[101:150]
#popD <-founders[151:200]
pop0 <- founders

# generation 0
gg <- 0
data <- collectData(pop=pop0, data=NULL, population="ABCD", generation = gg)
gg <- gg +1
for(i in 1:10){
  pop0 <- randCross(pop = pop0, nCrosses = 200)
  data <- collectData(pop=pop0, data=data, population="ABCD", generation = gg)
  gg <- gg +1
}
pop0 <- randCross(pop = pop0, nCrosses = 200)
pop0 <- setPheno(pop = pop0, varE = diag(varE))
popA <- pop0[1:50]
popB <- pop0[51:100]
popC <- pop0[101:150]
popD <- pop0[151:200]
data <- collectData(pop = popA, data = data, population = "A", generation = gg)
data <- collectData(pop = popB, data = data, population = "B", generation = gg)
data <- collectData(pop = popC, data = data, population = "C", generation = gg)
data <- collectData(pop = popD, data = data, population = "D", generation = gg)
gg <- gg +1
## GO ON HERE
# Select on each trait and keep the populations separate
ss <- 0 # sample size from other pops
for (generation in 12:20) {
  parentsA <- c(selectInd(pop = popA, nInd = 10, trait = 1),
                sample(popB, ss),
                sample(popC, ss),
                sample(popD, ss)
                )
  parentsB <- c(selectInd(pop = popB, nInd = 10, trait = 2),
                sample(popA, ss),
                sample(popC, ss),
                sample(popD, ss)
  )
  parentsC <- c(selectInd(pop = popC, nInd = 10, trait = 1),
                sample(popA, ss),
                sample(popB, ss),
                sample(popD, ss)
  )
  parentsD <- c(selectInd(pop = popD, nInd = 10, trait = 2),
                sample(popA, ss),
                sample(popB, ss),
                sample(popC, ss)
  )
  popA <- randCross(pop = parentsA, nCrosses = 50)
  popB <- randCross(pop = parentsB, nCrosses = 50)
  popC <- randCross(pop = parentsC, nCrosses = 50)
  popD <- randCross(pop = parentsD, nCrosses = 50)
  popA <- setPheno(pop = popA, varE = diag(varE))
  popB <- setPheno(pop = popB, varE = diag(varE))
  popC <- setPheno(pop = popC, varE = diag(varE))
  popD <- setPheno(pop = popD, varE = diag(varE))
  data <- collectData(pop = popA, data = data, population = "A", generation = generation)
  data <- collectData(pop = popB, data = data, population = "B", generation = generation)
  data <- collectData(pop = popC, data = data, population = "C", generation = generation)
  data <- collectData(pop = popD, data = data, population = "D", generation = generation)
  gg <- gg +1
}


#factor(data$pedigree$population)
table(data$pedigree$population)
data$pedigree$population <- factor(data$pedigree$population)
summary(data$pedigree$population)



# ---- eIBD COVARIANCE & PRECISION FACTOR --------------------------------

ped <- pedigree(sire = factor(data$pedigree$fid),
                dam = factor(data$pedigree$mid),
                label = factor(data$pedigree$id))

# We really need only A and Linv, but saving other objects for pedagogy.

# Relatedness - covariance (limited scalability)
pedA <- getA(ped = ped)
# ... factorisation
pedL <- getL(ped = ped)
pedT <- getT(ped = ped)
pedD <- getD(ped = ped)

# Relatedness - precision (very sparse so very scalable)
pedAInv <- getAInv(ped = ped)
# ... factorisation
pedLInv4 <- sparse2spam(getLInv(ped = ped))
pedTInv <- getTInv(ped = ped)
pedDInv <- getDInv(ped = ped)



#Matrix::writeMM(obj = pedLInv, file = "pedLInv4pop.mtx")
#write.csv(x = data$pedigree$population, file = "popLabel4pop.csv", row.names = FALSE)
#write.csv(x = data$pedigree, file = "pedMeta4pop.csv", row.names = FALSE)

pedMeta4 <- data$pedigree
pedGeno4 <- data$genoIBS
usethis::use_data(pedLInv4, overwrite = TRUE)
usethis::use_data(pedMeta4, overwrite = TRUE)
#usethis::use_data(pedGeno2, overwrite = TRUE)

#saveRDS(pedGeno4, "~/temp/pedGeno4.rds")
