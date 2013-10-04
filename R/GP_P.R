GP_P <-
function( EM_Cali, x, calc_var=FALSE, extra_output=FALSE) 
{
  if (! is.matrix(x)) x=as.matrix(x)
  X <- EM_Cali$X
  Y <- EM_Cali$Y
  R <- EM_Cali$R
  Rt <- EM_Cali$Rt
  R1X  <- EM_Cali$R1X
  R1tX <- EM_Cali$R1tX

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
  r  <- if (nx == 1) { matrix(cov_mat_1(lambda,X,x),n,1)
                     } else {cov_mat(lambda,X,x)}
  yp = mux %*% betahat + t(r) %*% solve ( Rt,  (Y- muX %*% betahat) )

  if (calc_var) # do we need the full output covariance matrix ? 
  {
    rr <-  if ((nx)==1) {1}   else {cov_mat(lambda,x ,x)}
    rr <- rr + diag(nx) * lambda$nugget

    # see email to Adrianikis and Challenor 5.08.2013 without response
    # whether one hould use R1X (i. e. : A-1) or R1tx (i. e. Atilde - 1)
    # new email on 22 august
    # ok. response on 9 septembre : Must use Atilde. 
    # equation modified accordingly
    P  <- ( mux - t(r) %*% R1tX )

    ## in the notation Oakley OHagan : 

   
    cxx_star = rr -  t(r) %*% solve ( Rt, r ) + 
             ( P  %*% solve ( (t(muX) %*% R1tX )  ,  t(P) ) )

    Sp = kronecker(cxx_star,sigma_hat_2) 
    Sp_diag = diag(Sp)
    
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
      
      tr <- a2s(t(r))
      tmux <- a2s(t(mux))
      Emul_pred = list(yp=yp, Sp=Sp, Sp_diag = Sp_diag, r=r, ht=mux, 
                       cxx = rr, 
                       hht = aperm ( outer(t(mux), mux) , c(2,3,1,4) ),
                       htt = s_aperm(s_outer(tmux, tr ), c(2,3,1,4)) , 
                       ttt = s_aperm(s_outer(a2s(r), tr )  , c(2,3,1,4))   
                       )
    }
    else
    Emul_pred = list(yp=yp, Sp=Sp, Sp_diag = Sp_diag)

  } else
  {
    rr =  1 + lambda$nugget
    dummy1 = (t(muX) %*% R1tX ) 
    Sp_diag=rep(0,nx)
    for (j in seq(1,nx))
      {
        rj = r[,j, drop=FALSE]
        P  <- ( mux[j,] - t(rj) %*% R1X )
        cxx = rr - (t(rj) %*% solve(Rt, rj)) + P %*%  solve ( dummy1 , t(P) )  
        Sp_diag[j] = sigma_hat_2 * cxx 
      }
    if (extra_output)
     Emul_pred = list(yp=yp, Sp_diag=Sp_diag, ht=mux, tt=t(r))
    else
     Emul_pred = list(yp=yp, Sp_diag=Sp_diag)
  }

  return(Emul_pred)

}
