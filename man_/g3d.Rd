\name{g3d}
\alias{g3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
g3d(x, y, h, omega)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{omega}{
%%     ~~Describe \code{omega} here~~
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
function (x, y, h, omega) 
{
    r <- lapply(dim(y), function(i) {
        seq(1, i - 1)
    })
    names(r) <- c("x", "y", "z")
    midvalue <- function(x) {
        (x[-length(x)] + x[-1])/2
    }
    newgrid <- as.data.frame(apply(x, 2, midvalue))
    print("newgrid")
    print(newgrid)
    print("dimy")
    print(dim(y))
    new_egrid <- lapply(expand.grid(newgrid), function(i) array(i, 
        dim(y) - 1))
    return(list(grid = newgrid, var = (omega$pre * abs((new_egrid$ecosw * 
        (y[r$x + 1, r$y, r$z] - y[r$x, r$y, r$z])) - (new_egrid$esinw * 
        (y[r$x, r$y + 1, r$z] - y[r$x, r$y, r$z]))) + (omega$eps * 
        abs(sqrt(1 - new_egrid$eps^2) * (y[r$x, r$y, r$z + 1] - 
            y[r$x, r$y, r$z]))))/h))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
