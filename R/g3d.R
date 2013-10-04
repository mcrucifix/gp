g3d <-
function(x,y, h, omega)
# 3d laplacien
{
 r <- lapply(dim(y), function(i) {seq(1,i-1)})

 names(r) <- c('x','y','z')
 midvalue <- function(x) {(x[-length(x)] + x[-1])/2.}
 newgrid <- as.data.frame(apply(x,2,midvalue))
 new_egrid <- lapply(expand.grid(newgrid), function(i) array(i,dim(y)-1))


# change caused by precession : y * d/dx * omega_x - x * d/dy * omega_y
# change caused by obliquity  : (1 - z^2) * d/dz * omega_z

# maximum possible change : sum of the absolute value of both

 return( 
    list(grid=newgrid, var = 

       #precession effects
       ( omega$pre *  abs ( 
       ( new_egrid$ecosw    *  ( y[r$x+1, r$y, r$z] - y[r$x, r$y, r$z] ) )     -
       ( new_egrid$esinw   *  ( y[r$x, r$y+1, r$z] - y[r$x, r$y, r$z] )  ) )   +  

       #obliquity effects
       ( omega$eps *  abs ( sqrt ( 1 - new_egrid$eps^2  ) * 
       ( y[r$x, r$y, r$z+1] - y[r$x, r$y, r$z] )   )  ) ) 

       # divide by interval 
       /  h )
       )
}
