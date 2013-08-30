gridpmp <-
function(m,pin,n) 
 {
  # build two grids (e.g. for intergration) over subspaces p and mp
  # p can be supplies as a vector of indices
  # e.g.: c(1,3). M is the total number of colons
  # n is a tuple with number of rows

  cX <-  seq(m)
  p <- list(p=pin, mp=cX[-pin])
  # 
  p <- lapply(p, function(p) { p; names(p)=p; p} )
  l <- lapply(p,length)

  # grids
  mpgrid = lapply(p, function(p) 
                 expand.grid(sapply(p, function(x) seq(n[x]), USE.NAMES=TRUE, simplify=FALSE)))
  mpgrid
  }
