marginals <- function (E, rhoarray, p, method='deterministic') {

  if (method != 'deterministic') 
  {
    message ('Only the deterministic method is supported in this package')
  }

    cgen <- function(x) {
        OUT = GP_P(E, x, calc_var = TRUE, extra_output = FALSE)
        list(Sp = OUT$Sp, yyt = drop(outer(OUT$yp, t(OUT$yp))))
    }
    sgen <- function(x) {
        OUT = GP_P(E, x, calc_var = FALSE, extra_output = FALSE)
        list(y = OUT$yp, y_m = OUT$yp_mean, y_g = OUT$yp_gaus, 
             yyt = OUT$yp^2, 
             yyt_m = OUT$yp_mean * OUT$yp_mean,
            yyt_mg = OUT$yp_mean * OUT$yp_gaus, 
             yyt_g = OUT$yp_gaus * OUT$yp_gaus)
    }
   
    x   <- expand.grid(attr(rhoarray,'levels'))
    double_integrals <- dintegral(rhoarray, x, cgen, p)
    simple_integrals <- sintegral_p ( rhoarray, x, sgen, p )



    OUT <- list(simple_integrals = simple_integrals, 
                double_integrals = double_integrals ) 
    OUT
}
