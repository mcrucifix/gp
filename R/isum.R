isum <- function(x) 
{ 
  # introduced sparse matrix support
  # s_apply and s_sum will by themselves redirect
  # to apply and sum if array is not of class
  # sparse_array
  n <-  ifelse (class(x) == "sparseMatrix", length(x$dims), length(dim(x))) 
  if (n  > 2) out <- s_apply(x, seq(n)[-c(1,2)], sum) else out <- s_sum(x) 
  if (class(out) == "sparseMatrix") 
  {
    s2a(out)
  } else
  {
    return (out)
  }
}


