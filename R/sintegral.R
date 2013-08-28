sintegral <-
function(rhoarray, X, cgen, p)
{
 f <- function(hmatrix) csum(sweep(hmatrix, 1, as.numeric(rhoarray),'*'))
 hmatrix <- cgen(X)
 if (is.list(hmatrix)) sapply(hmatrix, f, simplify=FALSE) else f(hmatrix)
}
