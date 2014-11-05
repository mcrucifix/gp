\name{pca}
\alias{pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standard PCA accounting for missing data
}
\description{
Standard (linear) PCA accounting for missing data, assuming a 2-d field
like structure for inputs.
}
\usage{
pca(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ Array, usually represinting a [longitude,latitudade, n] field, 
            where the two first dimensions represent the geographical distribution
            and the third dimension the sequence over which PCA is to be done 
            (time series or experiment ensemble ) 
           }
}
\details{
elementary SVD
}
\value{
List of
\item{mean}{2-d matrix of mean field over the third dimension}
\item{PCA}{Array \code{[lon, lat, n]}, where the third dimesion is the PC order, ranked
             in decreasing eigenvalue}
\item{amps}{Matrix \code{[n,n]} with each row corresponding to the field number, and
             columns to PCs}
\item{d}{Vector of eigenvalues}
}
\author{
Michel Crucifix
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
