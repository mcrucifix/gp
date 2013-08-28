convert <-
function(X, n) 
  { nc = ncol(X)
    rowSums(sapply(seq(nc), function(i) n^(i-1) * (X[[i]]-1) )) +1 
  }
