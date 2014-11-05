\name{L1O}
\alias{L1O}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Leave-one-out diagnostics for Gaussian process
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Provides leave-one-out diagnostics for Gaussian process 
given a design X with data Y and hyperparameters lambda. 
}
\usage{
L1O(X, Y, lambda, ...)
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
  \item{lambda}{
  List  made of (1) a vector code{theta}, with $m$ elements corresponding
  to roughness lengths associated with input variables,  and (2) 
%%     ~~Describe \code{lambda} here~~
}
  \item{...}{Other parameters passed to \code{\link{GP_C}}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Vector of 
  \item{Xi }{Input design }
  \item{Yi }{Supplied output}
  \item{mean }{GP output mean for each left-out experiment}
  \item{sd }{GP standard error for each left-out experiment}
  \item{mahalanobis }{GP mahalanobis for each left-out experiment}
 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Michel Crucifix
}
\note{
}

\seealso{ 
\code{\link{GP_C}}, \code{\link{GP_P}}
}
\examples{
 # univariate example
 X <- matrix(c(1,2,3,4,5,6,7), 7, 1)
 Y <- c(1.1, 2.1, 4.7, 1.3, 7.2, 8, 9)
 
 out <- L1O(X, Y, lambda=list(theta=1, nugget=0))
 
 require(Hmisc)

 with(out, errbar(Yi, mean, mean-sd, mean+sd))
 lines(c(-10,10),  c(-10,10), lty=2)
}


       
