\name{GP_C}
\alias{GP_C}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GP_C(X, Y, lambda, regress = "linear")
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
  \item{regress}{
%%     ~~Describe \code{regress} here~~
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
function (X, Y, lambda, regress = "linear") 
{
    funcmu = get(sprintf("funcmu_\%s", regress))
    if (!is.function(funcmu)) 
        stop("invalid regression model")
    if (!is.matrix(Y)) 
        Y = t(t(Y))
    if (!is.matrix(X)) 
        X = t(t(X))
    muX = t(apply(X, 1, funcmu))
    if (nrow(muX) == 1) 
        muX = (t(muX))
    n <- nrow(X)
    nn <- ncol(X)
    nbr <- n - ncol(X)
    nbrr <- n - ncol(X) - 2
    R <- cov_mat(lambda, X, X)
    R1X <- solve(R, muX)
    Rt <- R + diag(n) * lambda$nugget
    R1tX <- solve(Rt, muX)
    dummy1 = t(muX) \%*\% R1tX
    K = solve(dummy1, t(muX))
    betahat = K \%*\% solve(Rt, Y)
    dummy2 <- Y - muX \%*\% betahat
    e <- solve(Rt, dummy2)
    sigma_hat_2 <- t(dummy2) \%*\% e/nbrr
    M <- (lambda$nugget)^2 * t(dummy2) \%*\% solve(t(Rt) \%*\% Rt, 
        dummy2)/n
    Minfty <- t(dummy2) \%*\% dummy2/n
    log_REML = -1/2 * ((nbr) * log(diag(sigma_hat_2)) + log(det(Rt)) + 
        log(det(t(muX) \%*\% R1tX)))
    if (is.null(lambda$epsilon)) 
        lambda$epsilon = 1
    log_pen_REML = log_REML - 2 * (M/Minfty/lambda$epsilon)
    EM_Cali = list(betahat = betahat, sigma_hat_2 = sigma_hat_2, 
        R = R, Rt = Rt, muX = muX, X = X, Y = Y, lambda = lambda, 
        e = e, funcmu = funcmu, R1X = R1X, R1tX = R1tX, log_REML = log_REML, 
        log_pen_REML = log_pen_REML)
    return(EM_Cali)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
