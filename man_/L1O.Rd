\name{L1O}
\alias{L1O}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
L1O(X, Y, lambda, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, Y, lambda, ...) 
{
    if (!is.matrix(Y)) 
        Y = t(t(Y))
    if (!is.matrix(X)) 
        X = t(t(X))
    n <- nrow(X)
    as.data.frame(t(sapply(seq(1:n), function(i) {
        Xi <- X[-i, ]
        Yi <- Y[-i, ]
        Ei <- GP_C(Xi, Yi, lambda, ...)
        x = X[i, , drop = FALSE]
        OUT = GP_P(Ei, x)
        o_mean <- OUT$yp
        o_sd <- as.numeric(sqrt(OUT$Sp_diag))
        MD <- ((o_mean - Y[i, ])/o_sd)^2
        c(Xi = X[i, ], Yi = Y[i, ], mean = o_mean, sd = o_sd, 
            mahalanobis = MD)
    })))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
