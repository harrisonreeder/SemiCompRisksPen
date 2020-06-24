nll_pen_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                         hazard, frailty, model,
                         basis1, basis2, basis3, basis3_y1,
                         dbasis1, dbasis2, dbasis3,
                         penalty, lambda, a, penweights_list,
                         pen_mat_w_lambda, mu_smooth_fused){

  #general regularized ll function
  #this function should be able to return two possible things:
  #1. the standard penalized negative log likelihood, including parameter-wise penalties and fused penalties
  #2. the negative log likelihood, plus the original parameter-wise penalty, plus the nesterov-smoothed fused penalty
  #if the convex parameterwise lasso penalty is also being smoothed, then this function does not need to exist, because the 'smoothed_obj_func' above
  #is all you need.


  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  nll <- nll_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                  frailty=frailty, hazard=hazard, model=model,
                  basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)

  #here, we're not passing in the fused penalty, because whether mu_smooth_fused is 0 or nonzero, it is correctly estimated below
  pen <- pen_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,
                  penalty=penalty,lambda=lambda, a=a,
                  penweights_list=penweights_list,
                  pen_mat_w_lambda = NULL)

  out <- nll + (n * pen)

  #include nesterov smoothed result here instead of in above pen_func! This also returns the original fused penalty if mu_smooth_fused=0
  out <- out + n * smoothed_fused_lasso_func(para = para,
                                             pen_mat_w_lambda = pen_mat_w_lambda,
                                             mu_smooth_fused = mu_smooth_fused)

  return(out)
}
