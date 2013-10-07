varanal <- function( E, rhoarray, x, p)
{
  # fast version of varanal (that shoud eventually superseed varanal
  # this fast version does not use the 'extra_output' facility.
  # It rather uses the full co-variance matrix and meanstates
  # MCR Fri Oct  4 13:44:19 CEST 2013

  cgen <- function(x) 
  { 
    OUT = GP_P(E,x, calc_var=TRUE, extra_output=FALSE)
    list(   Sp = OUT$Sp, 
           yyt = drop(outer(OUT$yp,t(OUT$yp))) ,
         yyt_m  =   drop(outer(OUT$yp_mean,t(OUT$yp_mean))) ,
         yyt_mg = 2*drop(outer(OUT$yp_mean,t(OUT$yp_gaus))) ,
         yyt_g  =   drop(outer(OUT$yp_gaus,t(OUT$yp_gaus))) )
   }

  sgen <- function(x) 
  { 
     OUT = GP_P(E,x, calc_var=FALSE, extra_output=FALSE)
           list(  y = OUT$yp, 
                y_m = OUT$yp_mean, 
                y_g = OUT$yp_gaus 
                ) 
  }


  verb ( 'sintegral', { simple_integrals <- sintegral(rhoarray, x , sgen )  } ) 
  verb ( 'cintegral', { triple_integrals <- cintegral(rhoarray, x , cgen , p) }) 

  # original notation
  # split up emul variance into mean and gaussian process contributions 
  # added Fri Oct  4 14:14:50 CEST 2013
  Sp <-   triple_integrals[['Sp']]
  YY <-   triple_integrals[['yyt']]
  YY_m  <- triple_integrals[['yyt_m']]
  YY_mg <- triple_integrals[['yyt_mg']]
  YY_g  <- triple_integrals[['yyt_g']]

  Y    <- simple_integrals[['y']]
  Y_m  <- simple_integrals[['y_m']]
  Y_g  <- simple_integrals[['y_g']]

  OUT <- list (
    ev    = Sp,
    sv    = YY - Y*Y,
    sv_m  = YY_m - Y_m*Y_m,
    sv_mg = YY_mg - 2*Y_m*Y_g,
    sv_g  = YY_g - Y_g*Y_g 
  )
  # print output
  OUT 
}


