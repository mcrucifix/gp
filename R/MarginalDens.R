MarginalDens <-
function(rhoarray,p)
{
 # we suppose that the number of grid point is the same along all directions = n
 n <- dim(rhoarray)[1]
 m<- length(dim(rhoarray))
 Xgrid <-gridpmp(m,p,n)
 print('Xgrid...')
 print(Xgrid)

 out <- lapply(seq(nrow(Xgrid$p)), 
        function(i) 
        normalise(rhoarray[convert(combinepmp(Xgrid$mp, Xgrid$p[i,, drop=FALSE]),n)]))
 return(list(grid=Xgrid, marginals=out))
}
