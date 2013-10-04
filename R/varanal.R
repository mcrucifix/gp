varanal_deprecated <-
function( E, rhoarray, x, p)
  # as it is: high inefficiency due to the computation of very big
  # arrays which are, in fact, sparse (or at least dominated by a 
  # few elements (arrays computed within  GP_P)
  # NOTA : deprecated from version 0.1.2
  # as too inefficient. Replaced by varanal_fast
 
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
  # this is the squared mean associated with the simulator variance
  # i.e.  ( E*(E(Y)) )^2. 
  # However we also need 
  # E* ( E^2(Y) ),  to be subtracted from 'ev'. 
  # This quantity involves an integral of the covariance
  # structure over the whole domain, and  is
  # naturally obtained when the users choses 'NULL' as the 'p'
  # space. This is, though, _extremely_ computing expensive. 
  # so we don't want the user to compute this for every combinations
  # of the 'p'. 
  # this seems  however to be the only solution. 

  OUT <- list ( ev = emulator_variance, sv = simulator_variance -  squared_mean ) 
  OUT
  }
