combinepmp <-
function(Xp, xmp)
 {
  # combine a subspace with one element belonging to the remaining space 
  xmp <- as.matrix(xmp)
  Xmp <- matrix( xmp, nrow(Xp), length(xmp), byrow=TRUE, dimnames=list(NULL, colnames(xmp)))
  X <- cbind(Xp, Xmp)
  X <- X[, order(colnames(X))]
  X
}
