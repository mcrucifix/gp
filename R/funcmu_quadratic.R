funcmu_quadratic <-
function(x) 
				{ k=(1+ncol(x))/2
					x = cbind(1, x , x[k-1]*x[k], x[k-1]*x[k+1], x[k+1]*x[k])
					return(x)
				}
