
# pedigree from example data
ped <- pedigree(pedMeta$fid,
                pedMeta$mid,
                pedMeta$id)
# get additive relationship mat
AA <- getA(ped)
LL <- t(getL(ped))

# run naive PCA with and w/o centring
pcNaive <- prcomp(LL, center = F)
pcNaiveC <- prcomp(LL, center = T)
summary(pcNaive)$importance[,1:5]

summary(pcNaiveC)$importance[,1:5]

# get L and centre
LLc <- apply(t(getL(ped)), 2, function(x) x-mean(x))

image(LLc)

AAcent <- LLc %*% t(LLc)
sum(diag(AA))
sum(diag(AAcent))

pcNaiveLLcent <- prcomp(LLc, center = F)
summary(pcNaiveLLcent)$importance[,1:5]

sum(diag(AA))
sum(diag(AAcent))
sum(summary(pcNaive)$importance[1,]^2)*2649
sum(summary(pcNaiveC)$importance[1,]^2)*2649


# comparing plots
plot(pcNaive$x[,1:2])
plot(pcNaiveC$x[,1:2])
plot(pcNaiveAAcent$x[,1:2])

summary(pcNaive)$importance[,1:5]
summary(pcNaiveC)$importance[,1:5]
summary(pcNaiveLLcent)$importance[,1:5]


rpc01 <- rppca(ped, center = F)
rpc01c <- rppca(ped, center = T)
summary(rpc01) # matches naive approach
summary(rpc01c) # Stdev values match but proportions don't.
# This is because totVar is off - it does not apply after centring.


summary(rppca(ped, center = F))
summary(rppca(ped, center = F, totVar=3000))
summary(rppca(ped, center = T, totVar=3000))
summary(rppca(ped, center = T))
