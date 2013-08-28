norm <-
function (x, type = c("O", "I", "F", "M", "2")) 
{
    if (identical("2", type)) {
        svd(x, nu = 0L, nv = 0L)$d[1L]
    }
    else .Call("La_dlange", x, type, PACKAGE = "base")
}
