\name{define_mask}
\alias{define_mask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
define_mask()
}
%- maybe also 'usage' for other objects documented here.
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
function () 
{
    blank <- function() matrix(NA, 64, 32)
    tropforest <- blank()
    tropforest[53:58, 13:16] = 1
    tropforest[4:7, 13:16] = 1
    ssahel <- blank()
    ssahel[(62:70)\%\%64 + 1, 19:21] = 1
    arctic <- blank()
    arctic = matrix(0, 64, 32)
    arctic[-(50:60), 31:32] = 1
    arctic[, 32] = 1
    arctic[28:42, 30:31] = 1
    arctic[1:10, 30:31] = 1
    arctic[62:64, 30:31] = 1
    necan <- blank()
    necan[40:47, 27:29] = 1
    siberia <- blank()
    siberia[12:22, 27:29] = 1
    siberia[12:25, 27:29] = 1
    siberia[26:31, 28:29] = 1
    rm("blank")
    return(as.list(environment()))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
