meanrm <- function (x) mean (x, na.rm=TRUE) 

pca <-
function(A)
{
  nlon = dim(A)[1]
  nlat = dim(A)[2]
  nr = dim(A)[3]

  A = matrix(A, nlon*nlat, nr)

  means = apply(A,1,meanrm)
  Am = sweep( A,1, means) 

  o = which( ! is.na (means ) )

  U = svd(Am[o,])

  O = A*NA 
  O[o,] = U$u

  O = array(O, c(nlon,nlat,nr))

  # reformatting 

  means = matrix(means, nlon, nlat)

  d     = U$d

  return(list(mean=means, PCA=O, amps = U$v, d=d))

}
