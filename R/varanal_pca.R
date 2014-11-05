varanal_pca <-
function( EList, rhoarray, x, p)
  # original varanal routine recycled
  # for pca analysis
  # E is a list of emulators obtained for the different
  # eigenvalues (or, more generally, different components
  # of a field to be reconstructed by summing
  # ATTENTION : it is supposed that all emulators
  # have the same lambda (theta and nugget)
 
  {

  x   <- expand.grid(attr(rhoarray,'levels'))
  E <- EList[[1]]
  cgen <- function(x) { 
             OUT = GP_P(E,x, calc_var=TRUE, extra_output=TRUE)
             list(cxx = OUT$cxx, hht = OUT$hht,  htt = OUT$htt,  ttt = OUT$ttt ) }

  sgen <- function(x) { 
             OUT = GP_P(E,x, calc_var=FALSE, extra_output=TRUE)
             list(ht = OUT$ht, tt = OUT$tt) }


  verb ( 'sintegral', { simple_integrals <- sintegral(rhoarray, x , sgen )  } ) 
  verb ( 'cintegral', {triple_integrals <- cintegral(rhoarray, x , cgen , p) }) 

  # refer to the notation of Oackley and OHagan,  p. 762
  # all the integrals depend only on the emulator lambda, the training set 
  # and the design, hence only need to be computed once for all components
  Up <- triple_integrals[['cxx']]
  Pp <- triple_integrals[['ttt']]
  Qp <- triple_integrals[['hht']]
  Sp <- triple_integrals[['htt']]
  Rp <- simple_integrals[['ht']]
  Tp <- simple_integrals[['tt']]
  # now begins 
  A  <- E$Rt
  H  <- E$muX
  Winv <- t(H) %*%  solve(A,H)
  e <- E$e

  tr <- function(x) sum(diag(x))

  emulator_main_variance <- ( Up - tr(solve(A, Pp)) +
           tr( solve( Winv,  (Qp - Sp %*% solve(A,H) - t(H) %*% solve(A, t(Sp)) + 
               t(H) %*% solve (A, Pp)%*% solve ( A, H))) ))

  n <- length(EList)
  ev <- matrix(0,n,n)
  sv <- matrix(0,n,n); svm <- matrix(0,n,n); svg <- matrix(0,n,n)
  svmg <- matrix(0,n,n);

  for (i in seq(1,n))
    for (j in seq(1,n))
    {

        ev[i,j] = emulator_main_variance *
         t(EList[[i]]$Y - EList[[i]]$muX %*% EList[[i]]$betahat) %*% solve(A,
         (EList[[j]]$Y - EList[[j]]$muX %*% EList[[j]]$betahat))  / E$nbrr
  
        svg[i,j] <- tr(t(EList[[i]]$e) %*% Pp %*% EList[[j]]$e) -
                     (Tp %*% EList[[i]]$e) * (Tp %*% EList[[j]]$e)
        svmg[i,j] <- 2*tr(t(EList[[i]]$betahat) %*% Sp %*% EList[[j]]$e) -
                     2* ( Rp %*% EList[[i]]$betahat )*(Tp %*% EList[[j]]$e )
        svm[i,j] <-  tr(t(EList[[i]]$betahat) %*% Qp %*% EList[[j]]$betahat) -
                     (Rp %*% EList[[i]]$betahat) * (Rp %*% EList[[j]]$betahat)
     }
  OUT <- c ( ev = ev, sv=svg+svmg+svm, svm=svm, svmg=svmg, svg=svg)
  OUT
  }
