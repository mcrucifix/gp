sintegral_p <- function(rhoarray, X, cgen, p)
{
  # computes marginal integrals over rho_p

genmatrix <- function(cmatrix, ip=1)
 ## this is the service function computing the simple integral
 ## over -p|p
 {
     if (ip == 0)
      { rhomatrix <- 1 }
     else 
      { rhomatrix  <- as.numeric(rho[[ip]]$mpp) }

     if (class(cmatrix) == "list")
       # if cmatrix is a list, it must be named ! 
       # use the above rather than 'is.list' because
       # a sparseMatrix is also a list
       { 
         sapply ( cmatrix, function(cmatrix) 
           {  
            sum(cmatrix * rhomatrix)  
           },
           USE.NAMES=TRUE)
       }
       else 
       {
         sum(cmatrix * rhomatrix)  
       }
  }

n <- dim(rhoarray)
m <- length(dim(rhoarray))

rho <- MarginalDens(rhoarray,p)
Xgrid <- rho$grid
rho <- rho$marginals

sapply(seq(nrow(Xgrid$p)), 
    function(ip) 
  {  
    marginal <- Xgrid$p[ip,]
    names(marginal) <- colnames(Xgrid$p)
    grid <- convert(combinepmp(Xgrid$mp, Xgrid$p[ip,, drop=FALSE]),n)
    cmatrix   <- cgen(X[na.fail(grid),, drop=FALSE])    
    # generates matrix
    genmatrix(cmatrix, ip)
   }
 ) 

}

