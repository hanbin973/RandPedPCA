
library(randPedPCA)
head(pedMeta9)
table(pedMeta9$population)
ped9 <- pedigree(sire = pedMeta9$fid,
                 dam = pedMeta9$mid,
                 label = pedMeta9$id
                  )
pc9 <- rppca(ped9, center=T)
pc9nc <- rppca(ped9, center=F)


halfPairs <- function(rbound, col) {
  pairs(rbound[,1:4],
      upper.panel =function(x, y) {
        points(x[1:(length(x)/2)], y[1:(length(y)/2)], col=col)
        grid()
        },
      lower.panel = function(x, y) {
        points(x[(length(x)/2 + 1):length(x)], y[(length(y)/2+1):length(y)], col=col)
        grid()
      }
      )
}

halfPairs(rbind(pc9$x, pc9nc$x), pedMeta9$population)
halfPairs(rbind(pc9$x, pc9nc$x), rainbow(10)[as.numeric(pedMeta9$population)])
# compare to SNPs ---------------------------------------------------------

library(rsvd)


# down-sample SNPs to speed up PCA
g9 <- readRDS("~/temp/pedGeno9.rds")
dim(g9)
indivIndex <- (5670-2700):5670
g9 <- readRDS("~/temp/pedGeno9.rds")[indivIndex,sample(1:11000, 1000)]
dim(g9)
pcg9 <- rpca(g9, center=T, scale=F)
pcg9nc <- rpca(g9, center=F, scale=F)


halfPairs(rbind(pcg9$x, pcg9nc$x), col=rainbow(10)[as.numeric(pedMeta9$population)][indivIndex])
