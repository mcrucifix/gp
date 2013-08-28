\name{cintegral}
\alias{cintegral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cintegral(rhoarray, X, cgen, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rhoarray}{
%%     ~~Describe \code{rhoarray} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{cgen}{
%%     ~~Describe \code{cgen} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
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
function (rhoarray, X, cgen, p) 
{
    n <- dim(rhoarray)[1]
    m <- length(dim(rhoarray))
    lsum <- function(x) {
        inames <- names(x[[1]])
        shapes <- sapply(x[[1]], dim, simplify = FALSE, USE.NAMES = TRUE)
        if (is.null(inames)) {
            rsum(simplify(x))
        }
        else {
            sapply(inames, function(iname) {
                out <- rsum(sapply(x, function(x) {
                  x[[iname]]
                }))
                ishape <- shapes[[iname]]
                if (!is.null(ishape)) 
                  array(out, ishape)
                else out
            }, simplify = FALSE, USE.NAMES = TRUE)
        }
    }
    genmatrix <- function(ccmatrix, ip = 1) {
        rhorhomatrix <- outer(as.numeric(rho[[ip]]$mpp), as.numeric(rho[[ip]]$mpp))
        if (is.list(ccmatrix)) {
            sapply(ccmatrix, function(ccmatrix) {
                isum(sweep(ccmatrix, c(1, 2), rhorhomatrix, "*")) * 
                  rho[[ip]]$p
            }, USE.NAMES = TRUE, simplify = FALSE)
        }
        else {
            isum(sweep(ccmatrix, c(1, 2), rhorhomatrix, "*")) * 
                rho[[ip]]$p
        }
    }
    if (is.null(p)) {
        rho <- list(list(mpp = rhoarray/sum(rhoarray), p = 1))
        ccmatrix <- cgen(X)
        genmatrix(ccmatrix)
    }
    else {
        rho <- MarginalDens(rhoarray, p)
        Xgrid <- rho$grid
        rho <- rho$marginals
        lsum(lapply(seq(nrow(Xgrid$p)), function(ip) {
            marginal <- Xgrid$p[ip, ]
            names(marginal) <- colnames(Xgrid$p)
            grid <- convert(combinepmp(Xgrid$mp, Xgrid$p[ip, 
                , drop = FALSE]), n)
            ccmatrix <- cgen(X[grid, , drop = FALSE])
            genmatrix(ccmatrix, ip)
        }))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
