csum <-
function(x) { n = length(dim(x)) ; 
                      if (n  > 1) apply(x, seq(n)[-1], sum) else sum(x) }
