\name{MarginalDens}
\alias{MarginalDens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MarginalDens(rhoarray, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rhoarray}{
%%     ~~Describe \code{rhoarray} here~~
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
function (rhoarray, p) 
{
    n <- dim(rhoarray)[1]
    m <- length(dim(rhoarray))
    Xgrid <- gridpmp(m, p, n)
    print("Xgrid...")
    print(Xgrid)
    out <- lapply(seq(nrow(Xgrid$p)), function(i) normalise(rhoarray[convert(combinepmp(Xgrid$mp, 
        Xgrid$p[i, , drop = FALSE]), n)]))
    return(list(grid = Xgrid, marginals = out))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
