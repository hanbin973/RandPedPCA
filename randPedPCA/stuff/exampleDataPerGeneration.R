
pc <- rppca(pedLinv)
head(pedMeta)
summary(pedMeta$generation)


# get rang of generation values
genRan <- range(pedMeta$generation)

# turn into sequence
genSeq <- do.call(seq,as.list(genRan))


for(i in genSeq){
  print(i)
  png(filename = sprintf("Rplot%03d.png", i))
  xrange = range(pc$scores[,1])
  yrange = range(pc$scores[,2])

  plot(pc$scores[pedMeta$generation==i, 1],
       pc$scores[pedMeta$generation==i, 2],
       main=paste0("Generation ", i),
       xlim = xrange,
       ylim=yrange,
       xlab="PC1",
       ylab="PC2",
       col=factor(pedMeta$population)[pedMeta$generation==i]
  )
  grid()
  legend("topleft",
         pch=1,
         col=1:3,
         legend=levels(factor(pedMeta$population)),
         title = "Population"
         )


  dev.off()
}
# create gif using ImageMagick:
# e.g. magick -delay 10 *png ped.gif
