rsum <-
function(x) 
{ 
  n = length(dim(x)) 
  if (n  > 1) apply(x, 1, sum) else sum(x) 
}
