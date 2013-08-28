apply_mask <-
function(A, mask)
{
  # not weighted according to latitude... 
  OUT = sweep(A, c(1,2), mask, "*")
  apply(OUT, 3, mean, na.rm=TRUE)
}
