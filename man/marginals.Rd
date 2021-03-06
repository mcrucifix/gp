\name{marginals}
\alias{marginals}
\title{Simple and double intergrals for calculation of main effects}
\description{
Estimate simple and double integrals for estimators of marginal
mens and variances calculations for Gaussian process  main effect and 
associated variance analysis}
\usage{
marginals(E, rhoarray, p, method='deterministic')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{
  A Gaussian process list as generated by \code{\link{GP_C}}}
%%     ~~Describe \code{E} here~~
  \item{rhoarray}{
  An array of \code{n1 x n2 x n3 ... } with as many dimensions as the
  number of design factors used to generate the Gaussian process.
  The array represent the discretized density function. The densitity
  must be associated with a list of levels correspnding to the mid-points
  of the distribution. See examples. }
  \item{p}{
  A vector of the indices over which the integral of squares are being estimated
  (see details) }
  \item{method}{
  The \code{determistic} method assuming a sum over all components of \code{rhoarray}
   is currently the only one supported. 
  }
}
\details{
 \code{varanal} will estimate the emulator variances (more precisely, integral
 of squared quantities !) associated with the variance of indices
 complementary to \code{p} (\code{-p} in the notation of Oakley and Ohagan, 2004, p. 761),
 and then take the expectation of this quantity over \code{p}.
 the less indices in "p", the more computationally intensive
 given that all (-p) x (-p) combinations have to be summed up.
 variances measures will then be obtained as follows,

 total sensitivty index associated with, e.g., vars 1 and 2 (assuming 5 factors in total:
    varanal(E, rho, c(1,2,3,4,5))  - varanal(E, rho, c(3,4,5))

 mean sensitivty index associated with the same variables
     varanal(E, rho, c(1,2) ) - varanal (E, rho, NULL )


 Incidentially, varanal(E, rho, c(1,2,3,4,5)) - varanal(E, rho, NULL) will return the total
 variance. 

}
\value{
A list with following elements:
 \item{simple_integrals}{
 Matrix with as many columns as levels supplied along \code{p} 
 with lines corresponding to simple integrals 
 Simple integrals over \code{-p}  (the complementary 
 indices to \code{p}), of  \code{y}: the GP mean, 
 and its regression and stochastic components, \code{y_m} and \code{y_g}, 
 respectively. }
 \item{simple_integrals}{
  Matrix with as many columns as levels supplied along \code{p} 
 with lines corresponding to
 double integrals over \code{-p x -p}, composed
 of vectors with as many elements as levels supplied along \code{p} 
 with elements corresponding to double integrals over \code{Sp}: the GP co-variance,
 \code{yyt}: the GP mean products, and its contributions associated 
 with the regression, co-variance between regression and stochastic components,
 and stocastic compenents respectively: \code{yyt_m},  \code{yyt_mg},  \code{yyt_g}. }
}
\references{
}
\author{
Michel Crucifix
}
\seealso{
\code{\link{varanal} } }
\examples{
# generate data
X = as.matrix( expand.grid(seq(3), seq(3), seq(3))) 
Y = apply (X, 1, function(x) { x[1] - 2*x[1]*x[2] + x[3] - 2} )
Y = Y + rnorm(length(Y), sd=0.1)

E = GP_C(X, Y, lambda=list(theta=c(1,1,1), nugget=0.1))

# define rho function : uniform distribution

n = 8

levels = list ( seq(0,3, l=n), seq(0,3, l=n), seq(0,3, l=n) )
rho    = array(1, c(n,n,n))
rho    = rho / sum(rho) # not necessary. 
attr(rho, "levels") <- levels  # this is necessary ! 

# mean sensitivity index associated with the index 2

E_2 = marginals(E, rho,  c(2))

plot(E_2$simple_integrals['y',], type='l', main='Main effect of factor 2')
lines(E_2$simple_integrals['y',] +sqrt ( E_2$double_integrals['Sp',]  ), 
   lty=2, col='blue')

# variances associated with Gaussian process variance

lines(E_2$simple_integrals['y',] -sqrt ( E_2$double_integrals['Sp',]  ), 
   lty=2, col='blue')


# variances associated with process mean : 

lines(E_2$simple_integrals['y',] +sqrt ( 
    E_2$simple_integrals['yyt',]  - E_2$simple_integrals['y',]^2 ), 
   lty=2, col='red')

lines(E_2$simple_integrals['y',] -sqrt ( 
    E_2$simple_integrals['yyt',]  - E_2$simple_integrals['y',]^2 ), 
   lty=2, col='red')



legend('bottomleft', c('sd associated with GP variance', 'sd associated with process mean'), 
                 lty=2, col=c('red','blue'))



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
