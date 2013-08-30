cintegral <-
function(rhoarray, X, cgen, p)
{
 # X is the expanded grid filling space (assuming equal weights)
 # cgen is the outer generator 
 # rhoarray is an array of scalars
 n <- dim(rhoarray)
 m <- length(dim(rhoarray))
 
 lsum <- function(x) 
  # used to do the third integral, over x_p
  # considers the possibility that genmatrix has
  # in fact returned a list
   { inames <- names(x[[1]])
     shapes <- sapply(x[[1]], dim, simplify=FALSE, USE.NAMES=TRUE)
     if (is.null(inames)) 
        {
         rsum(simplify(x))
        } else 
        {
        sapply(inames, 
         function(iname) 
          { 
            out <- rsum( sapply(x, function(x) {x[[iname]] } ) ) 
            ishape <- shapes[[iname]]
            if (! is.null (ishape)) array(out, ishape) else out 
          } , 
         simplify=FALSE, USE.NAMES=TRUE) 
        }
  }

 genmatrix <- function(ccmatrix, ip=1)
   ## this is the service function computing the double integral
   ## over -p|p
   ## which will then need  to be intgrated over rho_p
   {
    rhorhomatrix  <- outer(as.numeric(rho[[ip]]$mpp), 
                            as.numeric(rho[[ip]]$mpp))   # matrix
    if (is.list(ccmatrix))
       # if ccmatrix is a list, it must be named ! 
       { 
       sapply ( ccmatrix, function(ccmatrix) 
               {  
                  isum(sweep(ccmatrix, c(1,2), rhorhomatrix,"*") ) *rho[[ip]]$p 
               },
              USE.NAMES=TRUE, simplify=FALSE)
       }
       else 
       {
          isum(sweep(ccmatrix, c(1,2), rhorhomatrix,"*") ) *rho[[ip]]$p 
       }
   }

 if (is.null(p))
 {
    # the user requests p = the entire domain. 
    # just the one provided. 
   rho <- list(list(mpp= rhoarray / sum(rhoarray), p=1))
   ccmatrix   <- cgen(X)                 # generates matrix
   genmatrix(ccmatrix)

 } else 
 {

   rho <- MarginalDens(rhoarray,p)
   Xgrid <- rho$grid
   rho <- rho$marginals
   lsum   ( 
   #that's the double integral for all elements of space p
   # should result in a vector of scalars
    lapply(seq(nrow(Xgrid$p)), 
    function(ip) 
         {  
           marginal <- Xgrid$p[ip,]
           names(marginal) <- colnames(Xgrid$p)
           grid <- convert(combinepmp(Xgrid$mp, Xgrid$p[ip,, drop=FALSE]),n)
           ccmatrix   <- cgen(X[na.fail(grid),, drop=FALSE])                 # generates matrix
           genmatrix(ccmatrix, ip)

   } ) ) 
 }
}
