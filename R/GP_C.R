# overloads exp function for creation of sparse matrices
# commented out on
# Fri Oct  4 14:33:38 CEST 2013 as at the end
# we won't use sparse matrices
# exp <- function(x) Vectorize({ifelse(x<(-2), 0, .Primitive("exp")(x)) })

cov_mat	= function(lambda=lambda,X1=X1,X2=X2, covar=exp)
{
  # requires lambda to be a vector
  theta = lambda$theta
  nk = length(lambda$theta)
  nx = nrow(X1)
  ny = nrow(X2)
  RR = array(0, c(nx,ny,nk))
  for ( k in seq(nk) ) { RR[,,k] = outer(X1[,k], X2[,k], "-") / theta[k]  }
  R = covar(- ( apply(RR^2,c(1,2),sum) ))
}
			   
cov_mat_1 = function(lambda=lambda,X1=X1,X2=X2, covar=exp)
{
  theta = lambda$theta
  # here need to adda check of consistency between theta and ncol X1 and nelem of X2
  nk = length(theta)
  nx = nrow(X1)
  vx = matrix(0,nx,nk)
  for  ( k in seq(1,nk) )
  {
    vx[,k] =  ((X1 [,k] - X2[k])/theta[k])^2
  }
  R = covar ( - ( apply(vx,1,sum)) )
}



GP_C <-
function( X, Y ,lambda, regress='linear', covar=exp )
  # revision history
  # 4.10.2013 : added covar option + passed in output.
  #             Backward campatible.
{
  funcmu = get(sprintf('funcmu_%s',regress))
  if  ( ! is.function(funcmu)) stop ('invalid regression model')
  if  ( ! is.matrix(Y)) Y=t(t(Y))
  if  ( ! is.matrix(X)) X=t(t(X))

	muX  = t(apply(X, 1, funcmu))
  # note : if only a single row gets out of this this
  # probably means that muX is in fact a constant
  if ( nrow(muX) == 1) muX = (t(muX))

	n    <- nrow(X)      
	nn   <- ncol(X)
	nbr  <- n - ncol(X)
	nbrr <- n - ncol(X) - 2 
	
	R   <- cov_mat(lambda,X,X,covar) 
  # R1X <- solve(R,muX) # for P matrix, disgarding the nugget

  # apply nugget (a la Andrianakis et Challenor)
  Rt   <- R + diag(n) * lambda$nugget
  R1tX <- solve(Rt,muX) # for P matrix 

  dummy1 =  t(muX) %*% R1tX
  K      =  solve( dummy1, t(muX)) 
	betahat = K %*% solve(Rt, Y)


  dummy2  <- Y - muX %*% betahat 
  e       <- solve(Rt, dummy2) 
  # e is the notation of Oakley and OHagan 2004
  # equivalent to Nabila's formulation
  # See Bastos and O'Hagan 2009 
	sigma_hat_2 <-  t( dummy2 ) %*%  e / nbrr

  # mean error at design points (proportionnal to nugget)

  M      <- (lambda$nugget)^2  * t(dummy2) %*% solve ( t(Rt) %*% Rt , dummy2 ) / n

  # max error at design points (when nugget tends to infinity, then just regression)

  Minfty <-  t(dummy2) %*% dummy2  / n

  # e.g. Baston and OHagan, eq. 15
  # Note : Adrianakis and Challenor have the same
  # result but multiplied by nbrr
  # (their definition of sigma_hat has this additional
  # factor, which propagates down to the likelihood)

  log_REML = -  1/2*( (nbr)*log(diag(sigma_hat_2)) 
                + log(det(Rt)) + log( det(t(muX) %*% R1tX )) )

  # penalised log lik (eq. 17 of Ardianakis and Challenor ) 
  # epsilon defaults to one 


  if (is.null(lambda$epsilon)) lambda$epsilon = 1 
  log_pen_REML = log_REML - 2 * (M / Minfty / lambda$epsilon )


  EM_Cali = list(betahat=betahat, sigma_hat_2=sigma_hat_2, 
                 R=R, Rt = Rt,  muX = muX, X=X, Y=Y, lambda=lambda, e=e,
                 funcmu=funcmu, # R1X = R1X, 
                 R1tX=R1tX , log_REML = log_REML, 
                 log_pen_REML=log_pen_REML, covar = covar, nbrr=nbrr )

  attr(EM_Cali, "class") <-  "GP_Emul"
	return(EM_Cali)
}


BIC.GP_Emul <- function(E) 
  {
    k <- with(E, length(lambda$theta) + length(lambda$nugget))
    n <- with(E, length(Y) )
    bic <- with (E,  -2 * log_REML +  k * (log ( n ) - log ( 2 * pi ) )  )
    return(bic)
  }

logLik.GP_Emul <- function(E) E$log_REML
