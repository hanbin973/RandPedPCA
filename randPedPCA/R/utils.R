
# Utility functions

# The oracle function that returns the value of A * G (but taking L^-1 as input)
oraculumLi <- function(Li, G){
  Y <- spam::backsolve(t(Li), G)
  return(spam::forwardsolve(Li, Y))
}
