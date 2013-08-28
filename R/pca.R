pca <-
function(A)
{
  nlon = dim(A)[1]
  nlat = dim(A)[2]
  nr = dim(A)[3]

  A = matrix(A, nlon*nlat, nr)

  mean = apply(A,1,mean)
  Am = sweep(A,1, apply(A,1,mean))

  o = which(! is.na(Am[,1]))

  U = svd(Am[o,])

  O = A; V=A
  O[o,] = U$u

  O = array(O, c(nlon,nlat,nr))
  mean = matrix(mean, nlon, nlat)

  d     = U$d

  return(list(mean=mean, PCA=O, amps = U$v, d=d))

}
