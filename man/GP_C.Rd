\name{GP_C}
\alias{GP_C}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gaussian process calibration 
}
\description{
Provides a Gaussian process calibrated on experiment design X
with data Y and hyperparameters lambda. Different  base regression
functions are available.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GP_C(X, Y, lambda, regress = "linear")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  Matrix of \code{n} lines corresponding to experiment members, 
  and with \code{m} columns corresponding to different inputs. 
  Can be one column, but always need to have this matrix form 
  (use \code{as.matrix} if needed )}
  \item{Y}{
  Currently the Gaussian process is univariate. So, vector
  with \code{n} elements corresponding to outputs. 
}
  \item{covar}{
  Minus the covariance function. Defaults to \code{exp}. 
}
  \item{lambda}{
  List  made of (1) a vector code{theta}, with $m$ elements corresponding
  to roughness lengths associated with input variables,  and (2) 
%%     ~~Describe \code{lambda} here~~
}
  \item{regress}{One of \code{constant}, \code{linear} or \code{quadratic}}
}
\details{
}
\value{
List of 
  \item{betahat }{Linear regression coefficients (posterior mode)}
  \item{sigma_hat_2 }{Simulator variance (posterior mode)}
  \item{R }{Design Covariance matrice (nugget free)}
  \item{Rt }{Design Covariance matrix (with nugget acconted for)}
  \item{muX }{Input matrix (regression function applied)}
  \item{X,Y,lambda, funcmu }{same as inputs}
  \item{R1X, R1tX }{Output dummies used by \code{GP_P}}
  \item{log_REML, loc_pen_REML,nbrr}{(penalised) log-likelihood}
}
\references{
Jeremy Oakley and Anthony O\'Hagan, Bayesian Inference for the Uncertainty Distribution of Computer Model Outputs,  Biometrika, 89, 769--784  2002

Ioannis Andrianakis and Peter G. Challenor, The effect of the nugget on Gaussian process emulators of computer models,  Computational Statistics \& Data Analysis, 56, 4215--4228  2012
}
\author{
Michel Crucifix
}
\note{
%%  ~~further notes~~
}


\seealso{ 
\code{\link{L1O}}, \code{\link{GP_P}}
}

\examples{
 # univariate example
 X <- matrix(c(1,2,3,4,5,6,7), 7, 1)
 Y <- c(1.1, 2.1, 4.7, 1.3, 7.2, 8, 9)
 
 
 # will attempt to optimize lambda for different
 # models and give the resulting log-likelihood
 # assumes no nugget 

# optimises the log to guarantee positiveness

  loglik <- function(theta)
 {
   -GP_C(X, Y, lambda=list(theta=exp(theta), nugget=0.0), regress=regress)$log_REML
 }

 for (r in c('constant','linear','quadratic'))
 {
  regress=r
  o = optimize(loglik, c(-1,1))
  print (sprintf (" Model \%s, Lambda=\%f, LogLik = \%f \n ", r, exp(o$minimum), 
  -o$objective) ) 
 }

}
