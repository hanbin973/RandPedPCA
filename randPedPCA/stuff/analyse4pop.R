dat4 <- read.table("pedMeta4pop.csv",
                   sep=",",
                   header=T)
library(randPedPCA)
head(dat4)
table(dat4$population)
ped4 <- pedigree(sire = dat4$fid,
                 dam = dat4$mid,
                 label = dat4$id
                  )
pc4 <- rppca(ped4, center=F)
pc5 <- rppca(ped4, center=T)
summary(pc4)
range(dat4$generation)
plot(pc4$x[,1:2],col=factor(dat4$population))
plot(pc4$x[,c(1,3)],col=factor(dat4$population))
plot(pc4$x[,c(1,4)],col=factor(dat4$population))
plot(pc4$x[,c(1,5)],col=factor(dat4$population))

plot(pc4$x[,2:3],col=factor(dat4$population))
plot(pc4$x[,3:4],col=factor(dat4$population))


plot(pc4$x[,1:2],col=heat.colors(21)[dat4$generation+1])
plot(pc4$x[,2:3],col=heat.colors(21)[dat4$generation+1])
plot(pc4$x[,3:4],col=heat.colors(21)[dat4$generation+1])
plot(pc4$x[,4:5],col=heat.colors(21)[dat4$generation+1])
plot(0:20, col=heat.colors(20)[0:20])
length(heat.colors(20)[0:20])
pairs(pc4$x[,1:4],
      col=factor(dat4$population),
      upper.panel = NULL)
up <- function(x,y){
  points(pc4[,x], pc4[,y], col=factor(dat4$population))
}

pairs(pc5$x[,1:4],
      col=factor(dat4$population),
      upper.panel = NULL)
traceback()

rbound <- rbind(pc4$x, pc5$x)
halfPairs <- function(rbound, metaData) {
  pairs(rbound[,1:4],
      upper.panel =function(x, y) {
        points(x[1:(length(x)/2)], y[1:(length(y)/2)], col=factor(metaData$population))
        grid()
        },
      lower.panel = function(x, y) {
        points(x[(length(x)/2 + 1):length(x)], y[(length(y)/2+1):length(y)], col=factor(metaData$population))
        grid()
      }
      )
}

halfPairs(rbind(pc4$x, pc4nc$x), pedMeta4)
# compare to SNPs ---------------------------------------------------------

library(rsvd)

ped <- pedigree(sire =pedMeta$fid,
                dam=pedMeta$mid,
                label=pedMeta$id)
ti <- getTInv(ped)
image(ti)
ti[100:180, 1:30]


?readRDS
load("data/pedMeta4.rda")
ped4 <- pedigree(sire=pedMeta4$fid,
                 dam=pedMeta4$mid,
                 label=pedMeta4$id)
pc4 <- rppca(ped4, center = T)
pc4nc <- rppca(ped4, center = F)
plot(pc4$x[,1:2], col=factor(pedMeta4$population))
plot(pc4$x[,2:3], col=factor(pedMeta4$population))

plot(pc4nc$x[,1:2], col=factor(pedMeta4$population))
plot(pc4nc$x[,2:3], col=factor(pedMeta4$population))

# down-sample SNPs to speed up PCA
g4 <- readRDS("~/temp/pedGeno4.rds")[,sample(1:11000, 1000)]
dim(g4)
pcg4 <- rpca(g4, center=T, scale=F)
pcg4nc <- rpca(g4, center=F, scale=F)
plot(pcg4$x[,1:2], col=pedMeta4$population)
plot(pcg4$x[,2:3], col=pedMeta4$population)
plot(pcg4nc$x[,1:2], col=pedMeta4$population)


halfPairs(rbind(pcg4$x, pcg4nc$x), pedMeta4)
