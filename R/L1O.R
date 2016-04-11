L1O <-
function( X, Y ,lambda, ... )
 {
    if  ( ! is.matrix(Y)) Y=t(t(Y))
    if  ( ! is.matrix(X)) X=t(t(X))

    n <- nrow(X)
    out <- as.data.frame(t(sapply(seq(1:n), function(i)
     {
       Xi <- X[-i,]
       Yi <- Y[-i,]
       Ei <- GP_C ( Xi, Yi, lambda, ... )
       x = X[i,, drop=FALSE]
       OUT = GP_P (Ei, x)
       o_mean <- OUT$yp 
       o_sd   <- as.numeric(sqrt(OUT$Sp_diag))
       MD  <- (  (o_mean -  Y[i,] ) /  o_sd  ) ^2
       c(Xi=X[i,], Yi=Y[i,], mean=o_mean, sd = o_sd, mahalanobis = MD )
     } )))
     attr(out, "class") <- 'L1O'
     return(out)

 }

plot.L1O <- function(L1O, ...)
{
   with(L1O, errbar (Yi, mean, mean-sd, mean+sd , ... ) )
   a1  <- min(L1O$mean)
   a2  <- max(L1O$mean)
   lines (c(a1,a2), c(a1,a2), lty=2 ) 
}
