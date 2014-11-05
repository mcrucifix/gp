isum <- function(x) 
{ 
  # introduced sparse matrix support
  # s_apply and s_sum will by themselves redirect
  # to apply and sum if array is not of class
  # sparse_array
  n = ifelse ( is.array(x), length(dim(x)), 0 ) 
  if (n  > 2) out <- apply(x, seq(n)[-c(1,2)], sum) else out <- sum(x) 
  return (out)
}


