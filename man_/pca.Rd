\name{pca}
\alias{pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pca(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
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
function (A) 
{
    nlon = dim(A)[1]
    nlat = dim(A)[2]
    nr = dim(A)[3]
    A = matrix(A, nlon * nlat, nr)
    mean = apply(A, 1, mean)
    Am = sweep(A, 1, apply(A, 1, mean))
    o = which(!is.na(Am[, 1]))
    U = svd(Am[o, ])
    O = A
    V = A
    O[o, ] = U$u
    O = array(O, c(nlon, nlat, nr))
    mean = matrix(mean, nlon, nlat)
    d = U$d
    return(list(mean = mean, PCA = O, amps = U$v, d = d))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
