dintegral <- function (rhoarray, X, cgen, p) 
{
    n <- dim(rhoarray)
    m <- length(dim(rhoarray))

    genmatrix <- function(ccmatrix, ip = 1) {
        rhorhomatrix <- outer(as.numeric(rho[[ip]]$mpp), as.numeric(rho[[ip]]$mpp))
        if (class(ccmatrix) == "list") {
            sapply(ccmatrix, function(ccmatrix) {
                isum(sweep(ccmatrix, c(1, 2), rhorhomatrix, "*")) 
            }, USE.NAMES = TRUE)
        }
        else {
            isum(sweep(ccmatrix, c(1, 2), rhorhomatrix, "*")) 
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
        print('ici')
        print(rho)
        print('-- ici --')
        sapply(seq(nrow(Xgrid$p)), function(ip) {
            marginal <- Xgrid$p[ip, ]
            names(marginal) <- colnames(Xgrid$p)
            grid <- convert(combinepmp(Xgrid$mp, Xgrid$p[ip, 
                , drop = FALSE]), n)
            ccmatrix <- cgen(X[na.fail(grid), , drop = FALSE])
            genmatrix(ccmatrix, ip)
        })
    }
}
