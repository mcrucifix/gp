# sparse matrix utility

i2a <- function(dims,m)
{
  # convert scalar indices into array indices
  m1 <- m-1
  n <- length(dims)
  c <- m1 %% dims[1]
  if (n > 2)
  {
    for (i in seq(2,(n-1)))
      c <- cbind(c, m1 %/% prod(dims[1:(i-1)]) %% (dims[i]))
  }
  c <- cbind(c,m1 %/% prod(dims[1:(n-1)]) ) + 1
  return(as.data.frame(c))
}


a2i <- function(X,n)
{
  if (is.null(n))
  { 
    return(X)
  } else 
  { 
  nc <- length(n)
  # convert array index (X) into integer index
  return(sum(sapply(seq(nc), function(i) enpower(n, i) * 
      (X[i] - 1))) + 1)
  }
}

a2s <- function(a)
{
  # convert dense array into sparse array
  a <- as.array(a)
  dims <- dim(a)
  m <- which(a != 0)
  indices <- i2a(dims, m)
  values  <- a[m]
  structure(list(indices=indices, values=values, dims=dim(a)), class="sparseMatrix")
}

s2a <- function(a)
{
  # convert sprase array into dense array
  dims <- a$dims
  M2 <- array(0, dim=a$dims)
  if (nrow(a$indices) >0) M2[convert(a$indices,dims)] = a$values
  M2
}

s_apply <- function(a, margins, fun,...)
{
  # equivalent of 'apply' function but for sparse arrrays
  # establish contigency table for all elements specified in margin
  # attention : func should return a scalar and this is not verified
  if (class(a) != "sparseMatrix") return (apply(a,margins,fun,...))
  dims          <- a$dims[margins]
  if (nrow(a$indices) == 0)
  {
    out = a; out$dims=dims; return(out)
  }
  itable1 <- a$indices[,margins,drop=FALSE]
  out_indices_c <- convert(itable1, dims)
  out_indices_u <- unique(out_indices_c)
  n <- length(out_indices_u)
  out <- list(indices = i2a(dims,out_indices_u), 
              values=rep(0,n), dims=a$dims[margins])
  for (j in seq(n))
  {
      values <- a$values[which(out_indices_c == out_indices_u[j])]
      out$values[j]  <- fun(values,...)
  }
  structure(out, class="sparseMatrix")
}
      
s_sum <- function(a)
{
  # sum for sparse matrices
  if (class(a) != "sparseMatrix") return (sum(a)) 
  sum(a$values)
}

s_sweep <- function(a, margins, stats, FUN='*',...)
{
  # equivalent of 'sweep' function but for sparse arrrays
  # establish contigency table for all elements specified in margin
  if (class(a) != "sparseMatrix") return (sweep(a,margins,stats,FUN,...))
  fun <- match.fun(FUN)
  # copy stats as a matrix, to have dims and all sorts
  stats <- as.matrix(stats)
  ds <- dim(stats)
  
  if ( any (ds != dim(a[margins]) )) 
    break ('stats or margins of wrong dimension ')
  itable1 <- a$indices[,margins,drop=FALSE]

  # look for relevant indices in stats matrix/array
  n <- nrow(itable1)
  out <- list(indices = a$indices, values=rep(0,length(a$indices)), dims=a$dims)
  if (n > 0)
  {
    stat_indices <- convert(itable1, dim(stats))
    out$values <- a$values * stats[stat_indices]
  }
  structure(out, class="sparseMatrix")
} 

s_aperm <- function(s, perm)
{
  # permute indices of sparse array
  # that's it : ultra fast. 
  n = ncol(s$indices)
  if (length(perm) != n) break("permute vector does not conform to matrix")
  if (length(perm) != length(unique(perm))) break("permute vector must include unique values")
  if (max(perm)    != n) break("permute vector excceds maximum allowed indice")
  out <- list(indices=s$indices[perm], values=s$values, dims=s$dims[perm])
  structure(out, class="sparseMatrix")
}

s_outer <- function(s1,s2, fun="*")
{
  # outer product between two sparse matrices
  # i dont think there is any significant bottleneck here
  s1size <- prod(s1$dims)
  ind_builder <- function(x,y) (x+s1size*(y-1) )
  dims    <- c(s1$dims,s2$dims)
  # if one of the two matrices is zero,
  # then the resulting outer matrix is also zero
  # with appropriate dimension
  if ((nrow (s1$indices) * nrow(s2$indices)) == 0)
  {
    return (list(indices=data.frame(NULL), values=NULL, 
                 dims=dims) )
  }
  # else, continue:
  indices <- i2a(dims, as.numeric(outer(convert(s1$indices,s1$dims), 
                             convert(s2$indices,s2$dims), ind_builder)))
  values  <- as.numeric(outer ( s1$values , s2$values, fun))
  out <- list(indices=indices, values=values, dims=dims)
  structure(out, class="sparseMatrix")
}

if (interactive())
{
  # create artificially sparse array

  M <- array(rnorm(120), dim=c(4,5,6))
  M[which(M[]>(-.4))] = 0. 

  M2 <- array(rnorm(40), dim=c(4,5,2))
  M2[which(M2[]>(-.4))] = 0. 

  stats=c(1,2,3,4)
  s <- a2s(M)

}


