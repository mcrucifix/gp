varanal <-
function( E, rhoarray, x, p)
{
  cgen <- function(x) { 
             OUT = GP_P(E,x, calc_var=TRUE, extra_output=TRUE)
             list(cxx = OUT$cxx, hht = OUT$hht,  htt = OUT$htt,  ttt = OUT$ttt ) }

  sgen <- function(x) { 
             OUT = GP_P(E,x, calc_var=FALSE, extra_output=TRUE)

             list(ht = OUT$ht, tt = OUT$tt) }


  verb ( 'sintegral', { simple_integrals <- sintegral(rhoarray, x , sgen )  } ) 
  verb ( 'cintegral', {triple_integrals <- cintegral(rhoarray, x , cgen , p) }) 

  # refer to the notation of Oackley and OHagan,  p. 762
  Up <- triple_integrals[['cxx']]
  Pp <- triple_integrals[['ttt']]
  Qp <- triple_integrals[['hht']]
  Sp <- triple_integrals[['htt']]
  Rp <- simple_integrals[['ht']]
  Tp <- simple_integrals[['tt']]
  A  <- E$Rt
  H  <- E$muX
  Winv <- t(H) %*%  solve(A,H)
  e <- E$e

  tr <- function(x) sum(diag(x))

  emulator_variance <- E$sigma_hat_2 * ( Up - tr(solve(A, Pp)) +
                       tr( solve( Winv,  (Qp - Sp %*% solve(A,H) - t(H) %*% solve(A, t(Sp)) + 
                               t(H) %*% solve (A, Pp)%*% solve ( A, H))) ))

  betahat <- E$betahat

  simulator_variance <- tr(t(e) %*% Pp %*% e) + 2*tr(t(betahat) %*% Sp %*% e) +
                        tr(t(betahat) %*% Qp %*% betahat)



  squared_mean <- (Rp %*% betahat + Tp %*% e ) ^2
  # the above is not entirely correct. This is  ( E*(E(Y)) )^2 whil we need
  # E* ( E^2(Y) ). Need to correct with variances as given in the top formula
  # of p. 761, in fact naturally obtained when the users choses 'NULL' as the 'p'
  # space. 

  OUT <- list ( ev = emulator_variance, sv = simulator_variance, sm = squared_mean ) 
  OUT
}
