


#' @method plot rppca
#' @export
plot.rppca <- function(x, dims=c(1,2), ...){
  ss <- summary(x)$importance
  if(dim(ss)[1] == 3){
    plot(
      x$x[,dims],
      xlab=paste0("PC", dims[1], " (", signif(100*ss[2,dims[1]], 3), "%)"),
      ylab=paste0("PC", dims[2], " (", signif(100*ss[2,dims[2]], 3), "%)")
    )
  } else {
    plot(
      x$x[,dims]
    )
  }
}
