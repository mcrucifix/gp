MarginalDens <-
function(rhoarray,p)
{
 n <- dim(rhoarray)
 m<- length(dim(rhoarray))
 Xgrid <-gridpmp(m,p,n)
 
 if (length(Xgrid$mp) == 0)
 {
  out<-lapply(rhoarray, function(i) list(mpp=1, p=i))
 } else
 {
 out <- lapply(seq(nrow(Xgrid$p)), 
        function(i) 
        normalise(na.fail(rhoarray[convert(combinepmp(Xgrid$mp, 
                  Xgrid$p[i,, drop=FALSE]),n)])))
 }
 return(list(grid=Xgrid, marginals=out))
}
