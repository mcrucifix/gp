\name{verb}
\alias{verb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
verb(fname, expr = expr, verbose = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fname}{
%%     ~~Describe \code{fname} here~~
}
  \item{expr}{
%%     ~~Describe \code{expr} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (fname, expr = expr, verbose = NULL) 
{
    if (is.null(verbose)) {
        if (exists("verbose", parent.env(environment()))) {
            verbose = get("verbose", parent.env(environment()))
        }
        else if (exists("verbose", .GlobalEnv)) {
            verbose = get("verbose", .GlobalEnv)
        }
        else break("verbose must be defined either in the GlobalEnv, in the\n                environment of the calling function, or as a parameter ")
    }
    if (verbose) {
        cat(sprintf("calling function \%s ...", fname))
        eval(expr)
        cat("... done \n")
    }
    else eval(expr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
