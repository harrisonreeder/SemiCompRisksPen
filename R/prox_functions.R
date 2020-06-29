##***********************##
##Helper Proximal operator functions##

lasso_prox_internal <- function(beta,lambda,step_size,penweights){

  #function to soft-threshold a vector at lambda*penweights*step_size
  #step_size is step size

  if(!is.null(penweights) && length(penweights)==length(beta)){
    return( as.matrix(  sign(beta)*pmax(0, abs(beta)-(step_size*lambda*penweights)) ) )
  } else{
    # warning("weights supplied to lasso function are NULL, or of incorrect length. ignoring weights.")
    return( as.matrix(  sign(beta)*pmax(0, abs(beta)-(step_size*lambda)) ) )
  }
}

prox_func <- function(para, prev_para, nP1, nP2, nP3, step_size,
                      penalty, lambda, penweights_list,
                      pen_mat_w,pen_mat_w_eig=NULL,lambda_f_vec,
                      mu_smooth_fused, ball_R=Inf){

  #perform proximal operator based on convex part of penalty term (following Yao (2018))

  # a actually never gets used here...
  # check_pen_params(penalty,penalty_fusedcoef,penalty_fusedbaseline,
  #                  lambda,lambda_fusedcoef,lambda_fusedbaseline,a)

  #redefine lambda and lambda_fusedcoef depending on a single value or three separate values are given
  if(length(lambda)==1){
    lambda1 <- lambda2 <- lambda3 <- lambda
  } else if(length(lambda)==3){
    lambda1 <- lambda[1]; lambda2 <- lambda[2]; lambda3 <- lambda[3]
  } else{ stop("lambda is neither a single value or a 3-vector!!") }

  nPtot <- length(para)
  nP0 <- length(para) - nP1 - nP2 - nP3

  ###COMPUTE THE FUSED LASSO PROXIMAL STEP USING ADMM
  eps_num <- min(sqrt(.Machine$double.eps), 100 * .Machine$double.eps) #definition from Smurf code
  if(!is.null(pen_mat_w) && mu_smooth_fused == 0){
    para_fl <- admm_po_cpp(beta_tilde = para,
                           slambda = lambda_f_vec * step_size,
                           penmat = pen_mat_w,
                           Q = if(!is.null(pen_mat_w_eig)) pen_mat_w_eig$Q else as.matrix(0),
                           eigval =  if(!is.null(pen_mat_w_eig)) pen_mat_w_eig$eigval else 0,
                           fast = if(!is.null(pen_mat_w_eig)) all(abs(pen_mat_w_eig$eigval) >= eps_num) else FALSE,
                           maxiter = 1e4, rho = 1,
                           beta_old = prev_para)
  } else{
    para_fl <- para
  }


  prox_out <- para_fl #i think this could equivalently be para but whatever

  #break out the beta vectors from the larger parameter vector, using nP0 to correctly pad out the baseline hazard and theta parameters
  if(nP1 != 0){
    beta1 <- para_fl[(1+nP0):(nP0+nP1)]
    prox_out[(1+nP0):(nP0+nP1)] <- lasso_prox_internal(beta=beta1,lambda=lambda1,step_size=step_size,penweights=penweights_list[["coef1"]])
  }
  if(nP2 != 0){
    beta2 <- para_fl[(1+nP0+nP1):(nP0+nP1+nP2)]
    prox_out[(1+nP0+nP1):(nP0+nP1+nP2)] <- lasso_prox_internal(beta=beta2,lambda=lambda2,step_size=step_size,penweights=penweights_list[["coef2"]])
  }
  if(nP3 != 0){
    beta3 <- para_fl[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    prox_out[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)] <- lasso_prox_internal(beta=beta3,lambda=lambda3,step_size=step_size,penweights=penweights_list[["coef3"]])
  }

  #now, add projection of the covariates onto the ball of radius R to potentially accomodate the constraints of Wang (2014)
  #this is equivalent to ridge regression
  if(nP1 + nP2 + nP3 > 0 && !is.null(ball_R) && !is.infinite(ball_R)){
    temp_norm2 <- sum(prox_out[(1+nP0):(nP0+nP1+nP2+nP3)]^2)
    if(temp_norm2 > ball_R){
      prox_out[(1+nP0):(nP0+nP1+nP2+nP3)] <- prox_out[(1+nP0):(nP0+nP1+nP2+nP3)] * ball_R / temp_norm2
    }
  }

  # names(prox_out) <- names(para)
  return(prox_out)
}
