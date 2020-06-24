smooth_obj_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                            hazard, frailty, model,
                            basis1, basis2, basis3, basis3_y1,
                            dbasis1, dbasis2, dbasis3,
                            penalty, lambda, a, penweights_list,
                            mu_smooth=0, #only set this to be non-zero if you are trying to work with FULLY SMOOTHED objective function
                            pen_mat_w_lambda, mu_smooth_fused){

  #function that combines the nll, the smooth concave part of any parameter-wise penalties, and optionally the nesterov-smoothed penalty
  #this function should be able to return two possible things:
  #1. The negative log-likelihood, plus the smooth nonconvex part of any penalties
  #this possibility is for the proximal gradient method that solves the fused lasso proximal operator,
  #so we just want the smooth part of the objective function
  #2. The above smooth function, plus the nesterov-smoothed lasso and fused lasso components


  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  out <- nll_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                  frailty=frailty, hazard=hazard, model=model,
                  basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)

  #note that this only accounts for the elementwise penalties, because the fused penalty is either convex, or smoothed in the next step
  if(penalty %in% c("scad","mcp")){
    componentwise_pen <- pen_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,
                                  penalty=penalty,lambda=lambda, a=a,
                                  penweights_list=penweights_list,
                                  pen_mat_w_lambda = NULL)

    neg_componentwise_pen <- pen_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,
                                      penalty="lasso",lambda=lambda, a=a,
                                      penweights_list=penweights_list,
                                      pen_mat_w_lambda = NULL)
    out <- out  + (n * componentwise_pen) - (n * neg_componentwise_pen)
  }

  #if we are also nesterov smoothing the remaining convex lasso penalty, include that!
  if(!is.null(mu_smooth) && mu_smooth > 0){
    out <- out + n * smoothed_lasso_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,lambda=lambda,
                                         mu_smooth=mu_smooth,penweights_list=penweights_list)
  }

  #if we are also nesterov smoothing the fused penalty, include that!
  if(!is.null(mu_smooth_fused) && mu_smooth_fused > 0){
    out <- out + n * smoothed_fused_lasso_func(para = para, pen_mat_w_lambda = pen_mat_w_lambda, mu_smooth_fused = mu_smooth_fused)
  }

  return(out)
}




smooth_obj_grad_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                                 hazard, frailty, model,
                                 basis1, basis2, basis3, basis3_y1,
                                 dbasis1, dbasis2, dbasis3,
                                 penalty, lambda, a,
                                 penweights_list, mu_smooth,
                                 pen_mat_w_lambda, mu_smooth_fused){
  #general gradient function of smooth part of penalized negative loglikelihood
  #corresponds to the smooth_obj_func from above

  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  out <- ngrad_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2, Xmat1=Xmat1, Xmat2=Xmat2,
                    Xmat3=Xmat3, frailty=frailty, hazard=hazard, model=model,
                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)

  if(tolower(penalty) %in% c("scad","mcp")){
    out <- out + n * pen_concave_part_prime_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,
                                                 penalty=penalty,lambda=lambda, a=a,
                                                 penweights_list=penweights_list)
  }

  #if we are also nesterov smoothing the remaining convex lasso penalty, include that!
  if(!is.null(mu_smooth) && mu_smooth > 0){
    out <- out + n * smoothed_lasso_prime_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,lambda=lambda,
                                               mu_smooth=mu_smooth, penweights_list=penweights_list)
  }

  if(!is.null(mu_smooth_fused) && mu_smooth_fused > 0){
    out <- out + n * smoothed_fused_lasso_prime_func(para=para, pen_mat_w_lambda=pen_mat_w_lambda, mu_smooth_fused=mu_smooth_fused)
  }

  return( out )

}





smooth_obj_lqa_pen_func <- function(para, prev_para,
                                    y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                                    hazard, frailty, model,
                                    basis1, basis2, basis3, basis3_y1,
                                    dbasis1, dbasis2, dbasis3,
                                    penalty, lambda, a, penweights_list,
                                    mu_smooth=0,  #only set this to be non-zero if you are trying to work with FULLY SMOOTHED objective function
                                    pen_mat_w_lambda, mu_smooth_fused,step_size){

  #function for the lasso-penalized local quadratic approximation of the smooth objective function
  #a quadratic expansion of the nll around prev_para, evaluated at para
  #following Wang (2014), equation 3.7


  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  #the original smooth objective at prev_para
  obj_val <- smooth_obj_func(para=prev_para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                             Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                             frailty=frailty, hazard=hazard, model=model,
                             basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                             dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                             penalty=penalty, lambda=lambda, a=a,
                             penweights_list=penweights_list, mu_smooth=mu_smooth,
                             pen_mat_w_lambda=pen_mat_w_lambda, mu_smooth_fused=mu_smooth_fused)

  #the corresponding grad of the smoothed objective function, evaluated at prev_para
  #oops this function is technically not defined until below but go look for it it's there!
  smooth_ngrad_temp <- smooth_obj_grad_func(para=prev_para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                            Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                            frailty=frailty, hazard=hazard, model=model,
                                            basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                            dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                            penalty=penalty, lambda=lambda, a=a,
                                            penweights_list=penweights_list, mu_smooth=mu_smooth,
                                            pen_mat_w_lambda=pen_mat_w_lambda, mu_smooth_fused=mu_smooth_fused)

  #the linear term of the taylor expansion
  out <- obj_val + t(smooth_ngrad_temp) %*% (para - prev_para)

  #the quadratic term of the expansion (notice I had to multiply by n, originally forgot that)
  out <- out + n/(2*step_size) * sum((para - prev_para)^2)

  #As a last step, add in any non-smooth penalties that remain after we've accounted for the smoothness above
  if(is.null(mu_smooth) || mu_smooth==0){
    pen <- pen_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,
                    penalty="lasso",lambda=lambda, a=a,
                    penweights_list=penweights_list,
                    pen_mat_w_lambda = NULL) #notice, we are separating out fused lasso piece, which is separately added below.
    out <- out + (n * pen)
  }

  #Separately, we add any non-smooth fused lasso that remains.
  if(is.null(mu_smooth_fused) || mu_smooth_fused==0){
    out <- out + n * smoothed_fused_lasso_func(para = para,
                                               pen_mat_w_lambda = pen_mat_w_lambda,
                                               mu_smooth_fused = mu_smooth_fused)
  }

  return(out)
}
