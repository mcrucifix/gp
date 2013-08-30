
enpower <- function(n, i) { 
          dn = length(n) ; out = 1 ; 
          if (i > 1) for (i in seq(1, i-1)) 
          out=out*n[i] ; out
          }

convert <-
function(X, n) 
  { nc = ncol(X)
    print(nc)
    if (nrow(X) > 1)
     return (rowSums(sapply(seq(nc), function(i) enpower(n, i) * (X[[i]]-1) )) +1 )
    else
     return (sum(sapply(seq(nc), function(i) enpower(n, i) * (X[[i]]-1) )) +1 )
  }

