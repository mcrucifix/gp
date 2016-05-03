GP_P <-
function( EM_Cali, x, calc_var=FALSE, extra_output=FALSE) 
{
  # revision history
  # 4.10.03 : added covar
  if (! is.matrix(x)) x=as.matrix(x)
  X <- EM_Cali$X
  Y <- EM_Cali$Y
  # the Rt includes the nugget 
  R <- EM_Cali$R
  e <- EM_Cali$e
  Rt <- EM_Cali$Rt
  #R1X  <- EM_Cali$R1X
  R1tX <- EM_Cali$R1tX
  covar  <- EM_Cali$covar

  lambda <- EM_Cali$lambda
  betahat   <- EM_Cali$betahat
  sigma_hat_2 <- EM_Cali$sigma_hat_2
  funcmu       <- EM_Cali$funcmu
  muX  <- EM_Cali$muX
  mux  <- t(apply(x, 1, funcmu)) 
  
  # if only a constant in returned need to transpose output
  if (length(funcmu(1)) == 1) {mux = t(mux)}

  nx   <- nrow(x)
  n    <- nrow(X)    
  # check if Nabila's code consistent with this. Was previousy
  # (in error) ncol(X)
  nq   <- length(betahat)
	nbrr <- n - nq - 2 

  # r is T(x) in the challenor paper. No nugget ! 
  r  <- if (nx == 1) { matrix(cov_mat_1(lambda,X,x,covar),n,1)
                     } else {cov_mat(lambda,X,x,covar)}

  yp_mean = mux %*% betahat 
  yp_gaus = t(r) %*% e

  yp = yp_mean + yp_gaus 

  if (calc_var) # do we need the full output covariance matrix ? 
  {
    rr <-  if ((nx)==1) {1}   else {cov_mat(lambda,x ,x,covar)}

    # save a copy of nugget-free covariance
    rr_nuggetfree <- rr

    rr <- rr + diag(nx) * lambda$nugget
    # correction : the nugget must not only be added to the diag, but in fact
    # to any couple of inputs that would be identical 

    Ind = cbind(which(duplicated(x)),nrow(x)+1-which(duplicated(apply(x, 2, rev) )))
    Ind = rbind(Ind, Ind[, c(2,1)])
    rr[Ind] = rr[Ind] + lambda$nugget 



    # see email to Adrianikis and Challenor 5.08.2013 without response
    # whether one hould use R1X (i. e. : A-1) or R1tx (i. e. Atilde - 1)
    # new email on 22 august
    # ok. response on 9 septembre : Must use Atilde. 
    # equation modified accordingly
    P  <- ( mux - t(r) %*% R1tX )

    ## in the notation Oakley OHagan : 

   
    cxx0 <- -  t(r) %*% solve ( Rt, r ) + ( P  %*% solve ( (t(muX) %*% R1tX )  ,  t(P) ) )

    cxx_star = rr  + cxx0
    cxx_star_nuggetfree = rr_nuggetfree  + cxx0

    Sp = cxx_star * as.numeric(sigma_hat_2)
    Sp_nuggetfree = cxx_star_nuggetfree * as.numeric(sigma_hat_2)

    Sp_diag = diag(Sp)
    Sp_diag_nuggetfree = diag(Sp_nuggetfree)

    if (extra_output)
    {
      # extra output useful to compute the triple integrals fro
      # general sensitivity analysis after Oakley and Ohagan
      # attention : our ht is the transpose of their 'h', but
      # easier for the calculations that will follow

      # output sparse matrices for cxx, htt, ttt
      # to do : this needs to be an option
      # i believe that for hht this may be counterproductive
      # or also if lambda is large so that in practice there are
      # little excluded points, or at last if the 
      # user does not cut exponential tails
      # to do: select sparse matrix or dense matrix
      # depending on th proportion of zeros in r. 
      
      # extra note on Fri Oct  4 14:04:32 CEST 2013
      # it seems that all of this is actually a waste of 
      # computing time : theses output are not used for
      # varanal_fast (which was verified to produce the same
      # output as varanal). 
 
      mux = as.matrix(mux)
      r = as.matrix(r)
      Emul_pred = list(yp=yp, yp_mean=yp_mean, yp_gaus=yp_gaus, Sp=Sp, Sp_nuggetfree = Sp_nuggetfree, Sp_diag = Sp_diag, Sp_diag_nuggetfree = Sp_diag_nuggetfree, r=r, ht=mux, 
                       cxx = rr, 
                       cxx_star = cxx_star,
                       hht = aperm ( outer(t(mux), mux) , c(2,3,1,4) ),
                       htt = aperm(outer(t(mux), t(r) ), c(2,3,1,4)) , 
                       ttt = aperm(outer(r, t(r) )  , c(2,3,1,4))   
                       )
    }
    else 
    Emul_pred = list(yp=yp, Sp=Sp, Sp_diag = Sp_diag, 
     Sp_diag_nuggetfree, yp_mean=yp_mean, yp_gaus=yp_gaus)

  } else
  {
    rr_nuggetfree = 1.
    rr =  rr_nuggetfree + lambda$nugget
    dummy1 = (t(muX) %*% R1tX ) 
    Sp_diag=rep(0,nx)
    Sp_diag_nuggetfree=rep(0,nx)
    for (j in seq(1,nx))
      {
        rj = r[,j, drop=FALSE]
        # bug corrected on Fri Nov 14 22:48:46 CET 2014
        # R1X below -> R1tX
        P  <- ( mux[j,] - t(rj) %*% R1tX )
        cxx0 = - (t(rj) %*% solve(Rt, rj)) + P %*%  solve ( dummy1 , t(P) )  
        cxx = rr + cxx0
        cxx_nuggetfree = rr_nuggetfree + cxx0
        Sp_diag[j] = sigma_hat_2 * cxx 
        Sp_diag_nuggetfree[j] = sigma_hat_2 * cxx_nuggetfree
      }
    if (extra_output) 
      Emul_pred = list(yp=yp, Sp_diag=Sp_diag, Sp_diag_nuggetfree = Sp_diag_nuggetfree, ht=mux, tt=t(r))
    else 
      Emul_pred = list(yp=yp, Sp_diag = Sp_diag,  Sp_diag_nuggetfree = Sp_diag_nuggetfree, yp_mean=yp_mean, yp_gaus=yp_gaus)
  }

}



