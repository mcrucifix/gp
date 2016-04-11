meanrm <- function (x) mean (x, na.rm=TRUE) 

pca <-
function(A)
{
  nlon = dim(A)[1]
  nlat = dim(A)[2]
  nr = dim(A)[3]

  A = matrix(A, nlon*nlat, nr)

  mean = apply(A,1,mean)
  Am = sweep( A,1, mean ) 

  o = which( ! is.na (mean ) )

  U = svd(Am[o,])

  O = A*NA 
  V = O 
  O[o,] = U$u

  O = array(O, c(nlon,nlat,nr))

  # reformatting 

  mean = matrix(mean, nlon, nlat)

  d     = U$d

  return(list(mean=mean, PCA=O, amps = U$v, d=d))

}
