\name{GP_P}
\alias{GP_P}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GP_P(EM_Cali, x, calc_var = FALSE, extra_output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EM_Cali}{
%%     ~~Describe \code{EM_Cali} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{calc_var}{
%%     ~~Describe \code{calc_var} here~~
}
  \item{extra_output}{
%%     ~~Describe \code{extra_output} here~~
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
function (EM_Cali, x, calc_var = FALSE, extra_output = FALSE) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    X <- EM_Cali$X
    Y <- EM_Cali$Y
    R <- EM_Cali$R
    Rt <- EM_Cali$Rt
    R1X <- EM_Cali$R1X
    R1tX <- EM_Cali$R1tX
    lambda <- EM_Cali$lambda
    betahat <- EM_Cali$betahat
    sigma_hat_2 <- EM_Cali$sigma_hat_2
    funcmu <- EM_Cali$funcmu
    muX <- EM_Cali$muX
    mux <- t(apply(x, 1, funcmu))
    if (length(funcmu(1)) == 1) {
        mux = t(mux)
    }
    nx <- nrow(x)
    n <- nrow(X)
    nq <- length(betahat)
    nbrr <- n - nq - 2
    r <- if (nx == 1) {
        matrix(cov_mat_1(lambda, X, x), n, 1)
    }
    else {
        cov_mat(lambda, X, x)
    }
    yp = mux \%*\% betahat + t(r) \%*\% solve(Rt, (Y - muX \%*\% betahat))
    if (calc_var) {
        rr <- if ((nx) == 1) {
            1
        }
        else {
            cov_mat(lambda, x, x)
        }
        rr <- rr + diag(nx) * lambda$nugget
        P <- (mux - t(r) \%*\% R1X)
        cxx = rr - t(r) \%*\% solve(Rt, r) + (P \%*\% solve((t(muX) \%*\% 
            R1tX), t(P)))
        Sp = kronecker(cxx, sigma_hat_2)
        Sp_diag = diag(Sp)
        if (extra_output) 
            Emul_pred = list(yp = yp, S = cxx, Sp = Sp, Sp_diag = Sp_diag, 
                cxx = cxx, r = r, ht = mux, hht = aperm(t(mux) \%o\% 
                  mux, c(2, 3, 1, 4)), htt = aperm(t(mux) \%o\% 
                  t(r), c(2, 3, 1, 4)), ttt = aperm(r \%o\% t(r), 
                  c(2, 3, 1, 4)))
        else Emul_pred = list(yp = yp, S = cxx, Sp = Sp, Sp_diag = Sp_diag)
    }
    else {
        rr = 1 + lambda$nugget
        dummy1 = (t(muX) \%*\% R1tX)
        Sp_diag = rep(0, nx)
        for (j in seq(1, nx)) {
            rj = r[, j, drop = FALSE]
            P <- (mux[j, ] - t(rj) \%*\% R1X)
            cxx = rr - (t(rj) \%*\% solve(Rt, rj)) + P \%*\% solve(dummy1, 
                t(P))
            Sp_diag[j] = sigma_hat_2 * cxx
        }
        if (extra_output) 
            Emul_pred = list(yp = yp, Sp_diag = Sp_diag, ht = mux, 
                tt = t(r))
        else Emul_pred = list(yp = yp, Sp_diag = Sp_diag)
    }
    return(Emul_pred)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
