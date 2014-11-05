funcmu_quadratic <-
function(x) 
				{ 
          out = funcmu_linear(x)
          n = length(x)
          for (i in seq(n)) for (j in seq(i,n)) 
          out = c(out, x[i] * x[j])
					return(out)
				}
