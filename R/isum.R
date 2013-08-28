isum <-
function(x) { n = length(dim(x)) ; 
                      if (n  > 2) apply(x, seq(n)[-c(1,2)], sum) else sum(x) }
