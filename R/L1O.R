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
       c(Xi=X[i,], Yi=as.numeric(Y[i,]), mean=o_mean, sd = o_sd, mahalanobis = MD )
     } )))
     attr(out, "class") <- 'L1O'
     attr(out, "scaled:center") = attr(Y,"scaled:center")   
     attr(out, "scaled:scale")  = attr(Y,"scaled:scale")   
     return(out)

 }

plot.L1O <- function(L1O, extended=FALSE, rescale = FALSE, ilog = FALSE, ...)
{
   Yi <- L1O$Yi 
   sd <- L1O$sd
   mean <- L1O$mean


   if (rescale)
   {
    Yi <- Yi * attr(L1O, "scaled:scale") + attr(L1O, "scaled:center")
    sd <- sd * attr(L1O, "scaled:scale") 
    mean <- mean * attr(L1O, "scaled:scale") + attr(L1O, "scaled:center")
   }

   if (ilog)
   {
    Yi <- exp(Yi)
    mean <- exp(mean + sd^2/2)
    sd   <- mean * sqrt ( exp ( sd^2 ) - 1 ) 

   }


   if (extended)
   {
     cmin <- min(Yi, mean-sd)
     cmax <- max(Yi, mean+sd)
     dc   <- cmax - cmin

     errbar(mean, Yi, Yi-sd, Yi+sd, main = "Predicted vs true")
     k <-  ( ks.test ( (mean - Yi) / sd , pnorm ) )
     lines(c(cmin,cmax), c(cmin,cmax))

     qqnorm( (mean - Yi) / sd)    # should be a straight line
     lines(c(cmin,cmax), c(cmin,cmax), lty=2)

     plot(mean, (mean - Yi) / sd , main = "standard error", ...)
     abline(h=0, lty=2, col='blue')
     text (cmin+dc*0.6, cmin+dc*0.9, sprintf('p = %.2f',k$p.value), pos = 4)
   } else 
   {
     cmin <- min(Yi, mean-sd)
     cmax <- max(Yi, mean+sd)
     errbar (Yi, mean, mean-sd, mean+sd , ... ) 
     lines(c(cmin,cmax), c(cmin,cmax))
   }
}
