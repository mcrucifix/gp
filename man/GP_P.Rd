\name{GP_P}
\alias{GP_P}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Provides output of calibrated gaussian process
}
\description{
Provides output of calibrated gaussian process
}
\usage{
GP_P(EM_Cali, x, calc_var = FALSE, extra_output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EM_Cali}{
   List output of \code{\link{GP_C}}
}
  \item{x}{
   Input 'design' for which the gaussian process will be estimated.
%%     ~~Describe \code{x} here~~
}

  \item{calc_var}{

  Logical. Set to \code{TRUE} means that the full covariance matrix
  associated with \code{x} will be estimated. More time consiming. 
  Defaults to \code{FALSE}
%
%%     ~~Describe \code{extra_output} here~~
}

  \item{extra_output}{
  Logical. \code{TRUE} if extra diagnostic output associated 
  with the different contributions to output are needed. 
%%     ~~Describe \code{calc_var} here~~
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List of
 \item{yp }{Output mean at points \code{x} }
 \item{Sp_diag }{Diagonal of co-variance matrix}
 \item{yp_mean, yp_gauss }{Mean and Gaussian_process contribution to output}
 \item{Sp}{If requested, full covariance matrix}
 \item{ht, cxx, cxx_star, hht, htt, ttt}{If requested (\code{extra_output = TRUE}), 
 additional outputs associated to GP estimation. See Oakley and OHagan for expl. }
 }
\references{
Jeremy Oakley and Anthony O\'Hagan, Bayesian Inference for the Uncertainty Distribution of Computer Model Outputs,  Biometrika, 89, 769--784  2002

Ioannis Andrianakis and Peter G. Challenor, The effect of the nugget on Gaussian process emulators of computer models,  Computational Statistics \& Data Analysis, 56, 4215--4228  2012
}
\author{
Michel Crucifix
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
 X <- matrix(c(1,2,3,4,5,6,7), 7, 1)
 Y <- c(1.1, 2.1, 4.7, 1.3, 7.2, 8, 6)
 x <- seq(0,9,0.01)
 x <- matrix(x, length(x), 1)
  
 # comparse constant and linear regression
 
 models = c(constant='constant',linear='linear')
 colors = c(constant='blue', linear='red')

 E <- lapply(models, function(m)  GP_C(X, Y, lambda=list(theta=1, nugget=0.1), regress=m) )
 O <- lapply(models, function(m) GP_P(E[[m]], x) )

 plot(X, Y, xlim=c(0,10), ylim=c(0,9))
 
 for (m in models)
 {
 lines(x, O[[m]]$yp, col=colors[m])
 lines(x, O[[m]]$yp + sqrt(O[[m]]$Sp), lty=2, col=colors[m])
 lines(x, O[[m]]$yp - sqrt(O[[m]]$Sp), lty=2, col=colors[m])
 }


 legend('topleft', models, col=colors, lty=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

