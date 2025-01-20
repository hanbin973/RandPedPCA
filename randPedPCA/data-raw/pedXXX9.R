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
founderGenomes <- runMacs(nInd = 270, nChr = 10, segSites = 1100,
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
  pop0 <- randCross(pop = pop0, nCrosses = 270)
  data <- collectData(pop=pop0, data=data, population="ABCD", generation = gg)
  gg <- gg +1
}
pop0 <- randCross(pop = pop0, nCrosses = 270)
pop0 <- setPheno(pop = pop0, varE = diag(varE))
pop1 <- pop0[1:30]
pop2 <- pop0[31:60]
pop3 <- pop0[61:90]
pop4 <- pop0[91:120]
pop5 <- pop0[121:150]
pop6 <- pop0[151:180]
pop7 <- pop0[181:210]
pop8 <- pop0[211:240]
pop9 <- pop0[241:270]
data <- collectData(pop = pop1, data = data, population = "p1", generation = gg)
data <- collectData(pop = pop2, data = data, population = "p2", generation = gg)
data <- collectData(pop = pop3, data = data, population = "p3", generation = gg)
data <- collectData(pop = pop4, data = data, population = "p4", generation = gg)
data <- collectData(pop = pop5, data = data, population = "p5", generation = gg)
data <- collectData(pop = pop6, data = data, population = "p6", generation = gg)
data <- collectData(pop = pop7, data = data, population = "p7", generation = gg)
data <- collectData(pop = pop8, data = data, population = "p8", generation = gg)
data <- collectData(pop = pop9, data = data, population = "p9", generation = gg)
gg <- gg +1
## GO ON HERE
# Select on each trait and keep the populations separate
ss <- 0 # sample size from other pops
for (generation in 12:20) {
  parents1 <- c(selectInd(pop = pop1, nInd = 20, trait = 1),
                sample(pop2, ss),
                sample(pop4, ss),
                sample(pop5, ss)
                )
  parents2 <- c(selectInd(pop = pop2, nInd = 20, trait = 1),
                sample(pop1, ss),
                sample(pop3, ss),
                sample(pop5, ss)
  )
  parents3 <- c(selectInd(pop = pop3, nInd = 20, trait = 1),
                sample(pop2, ss),
                sample(pop5, ss),
                sample(pop6, ss)
  )
  parents4 <- c(selectInd(pop = pop4, nInd = 20, trait = 1),
                sample(pop1, ss),
                sample(pop5, ss),
                sample(pop7, ss)
  )
  parents5 <- c(selectInd(pop = pop5, nInd = 20, trait = 1),
                sample(pop2, ss),
                sample(pop4, ss),
                sample(pop6, ss),
                sample(pop8, ss)
  )
  parents6 <- c(selectInd(pop = pop6, nInd = 20, trait = 1),
                sample(pop3, ss),
                sample(pop5, ss),
                sample(pop9, ss)
  )
  parents7 <- c(selectInd(pop = pop7, nInd = 20, trait = 1),
                sample(pop4, ss),
                sample(pop5, ss),
                sample(pop6, ss)
  )
  parents8 <- c(selectInd(pop = pop8, nInd = 20, trait = 1),
                sample(pop5, ss),
                sample(pop7, ss),
                sample(pop9, ss)
  )
  parents9 <- c(selectInd(pop = pop9, nInd = 20, trait = 1),
                sample(pop5, ss),
                sample(pop6, ss),
                sample(pop8, ss)
  )
  pop1 <- randCross(pop = parents1, nCrosses = 30)
  pop2 <- randCross(pop = parents2, nCrosses = 30)
  pop3 <- randCross(pop = parents3, nCrosses = 30)
  pop4 <- randCross(pop = parents4, nCrosses = 30)
  pop5 <- randCross(pop = parents5, nCrosses = 30)
  pop6 <- randCross(pop = parents6, nCrosses = 30)
  pop7 <- randCross(pop = parents7, nCrosses = 30)
  pop8 <- randCross(pop = parents8, nCrosses = 30)
  pop9 <- randCross(pop = parents9, nCrosses = 30)

  pop1 <- setPheno(pop = pop1, varE = diag(varE))
  pop2 <- setPheno(pop = pop2, varE = diag(varE))
  pop3 <- setPheno(pop = pop3, varE = diag(varE))
  pop4 <- setPheno(pop = pop4, varE = diag(varE))
  pop5 <- setPheno(pop = pop5, varE = diag(varE))
  pop6 <- setPheno(pop = pop6, varE = diag(varE))
  pop7 <- setPheno(pop = pop7, varE = diag(varE))
  pop8 <- setPheno(pop = pop8, varE = diag(varE))
  pop9 <- setPheno(pop = pop9, varE = diag(varE))

  data <- collectData(pop = pop1, data = data, population = "p1", generation = generation)
  data <- collectData(pop = pop2, data = data, population = "p2", generation = generation)
  data <- collectData(pop = pop3, data = data, population = "p3", generation = generation)
  data <- collectData(pop = pop4, data = data, population = "p4", generation = generation)
  data <- collectData(pop = pop5, data = data, population = "p5", generation = generation)
  data <- collectData(pop = pop6, data = data, population = "p6", generation = generation)
  data <- collectData(pop = pop7, data = data, population = "p7", generation = generation)
  data <- collectData(pop = pop8, data = data, population = "p8", generation = generation)
  data <- collectData(pop = pop9, data = data, population = "p9", generation = generation)
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
#pedL <- getL(ped = ped)
#pedT <- getT(ped = ped)
#pedD <- getD(ped = ped)

# Relatedness - precision (very sparse so very scalable)
#pedAInv <- getAInv(ped = ped)
# ... factorisation
pedLInv9 <- getLInv(ped = ped)
#pedTInv <- getTInv(ped = ped)
#pedDInv <- getDInv(ped = ped)



#Matrix::writeMM(obj = pedLInv, file = "pedLInv4pop.mtx")
#write.csv(x = data$pedigree$population, file = "popLabel4pop.csv", row.names = FALSE)
#write.csv(x = data$pedigree, file = "pedMeta4pop.csv", row.names = FALSE)

pedMeta9 <- data$pedigree
pedGeno9 <- data$genoIBS
usethis::use_data(pedLInv9, overwrite = TRUE)
usethis::use_data(pedMeta9, overwrite = TRUE)
#usethis::use_data(pedGeno2, overwrite = TRUE)

#saveRDS(pedGeno9, "~/temp/pedGeno9.rds")
