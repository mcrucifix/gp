L1O <-
function( X, Y ,lambda, ... )
 {
    if  ( ! is.matrix(Y)) Y=t(t(Y))
    if  ( ! is.matrix(X)) X=t(t(X))

    n <- nrow(X)
    as.data.frame(t(sapply(seq(1:n), function(i)
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


 }
