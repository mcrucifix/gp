l3d <-
function(x,h)
# 3d laplacien
{
 r <- lapply(dim(x), function(i) {seq(2,i-1)})
 names(r) <- c('x','y','z')
 return(( x[r$x+1, r$y, r$z] +
        x[r$x, r$y+1, r$z] +
        x[r$x, r$y, r$z+1] +
        x[r$x-1, r$y, r$z] +
        x[r$x, r$y-1, r$z] +
        x[r$x, r$y, r$z-1] -
       6*x[r$x, r$y, r$z]  ) / (h^2) )
}
