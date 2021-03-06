verb <-
function (fname, expr=expr, verbose=NULL)
{
 if ( is.null (verbose))
  # try to catch  verbose from the parent environment
  {
   if ( exists('verbose', parent.env(environment())) )
    { verbose  = get('verbose', parent.env(environment())) }
   else if ( exists('verbose', .GlobalEnv ))
   # try to catch it from the global environment
    { verbose = get('verbose', .GlobalEnv) }
   else verbose = FALSE
  } 
  if (verbose) 
        {
         cat (sprintf('calling function %s ...', fname))
         eval(expr)
         cat ('... done \n')
        } else eval (expr) 
}
