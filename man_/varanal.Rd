\name{varanal}
\alias{varanal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
varanal(E, rhoarray, x, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{
%%     ~~Describe \code{E} here~~
}
  \item{rhoarray}{
%%     ~~Describe \code{rhoarray} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (E, rhoarray, x, p) 
{
    test <- FALSE
    cgen <- function(x) {
        OUT = GP_P(E, x, calc_var = TRUE, extra_output = TRUE)
        list(cxx = OUT$cxx, hht = OUT$hht, htt = OUT$htt, ttt = OUT$ttt)
    }
    sgen <- function(x) {
        OUT = GP_P(E, x, calc_var = FALSE, extra_output = TRUE)
        list(ht = OUT$ht, tt = OUT$tt)
    }
    verb("sintegral", {
        simple_integrals <- sintegral(rhoarray, x, sgen)
    })
    verb("cintegral", {
        triple_integrals <- cintegral(rhoarray, x, cgen, p)
    })
    Up <- triple_integrals[["cxx"]]
    Pp <- triple_integrals[["ttt"]]
    Qp <- triple_integrals[["hht"]]
    Sp <- triple_integrals[["htt"]]
    Rp <- simple_integrals[["ht"]]
    Tp <- simple_integrals[["tt"]]
    A <- E$Rt
    H <- E$muX
    Winv <- t(H) \%*\% solve(A, H)
    e <- E$e
    tr <- function(x) sum(diag(x))
    emulator_variance <- E$sigma_hat_2 * (Up - tr(solve(A, Pp)) + 
        tr(solve(Winv, (Qp - Sp \%*\% solve(A, H) + t(H) \%*\% solve(A, 
            t(Sp)) + t(H) \%*\% solve(A, Pp) \%*\% solve(A, H)))))
    betahat <- E$betahat
    simulator_variance <- tr(t(e) \%*\% Pp \%*\% e) + 2 * tr(t(betahat) \%*\% 
        Sp \%*\% e) + tr(t(betahat) \%*\% Qp \%*\% betahat)
    squared_mean <- (Rp \%*\% betahat + Tp \%*\% e)^2
    OUT <- list(ev = emulator_variance, sv = simulator_variance, 
        sm = squared_mean)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
