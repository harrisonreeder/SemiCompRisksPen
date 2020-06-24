
##Nesterov-smoothed lasso and fused lasso functions##
##*************************************************##

#helper function for the nesterov smoothed approximation of the lasso
#Defined in Chen (2012) equation 3.5, with the solution alphastar defined
#in proposition 2 (no extension to general penalties, for which we specifically use the fused version below)
smoothed_lasso_internal <- function(para, lambda, mu_smooth, penweights){

  if(is.null(mu_smooth)){
    return(0)
  }

  if(!is.null(penweights) && length(penweights)==length(para)){
    temp_vec <- abs(para) * lambda * penweights
  } else{
    temp_vec <- abs(para) * lambda
  }
  if(mu_smooth==0){
    return(sum(temp_vec))
  } else{
    alphastar <-  sign(temp_vec/mu_smooth)*pmin(1, abs(temp_vec/mu_smooth))
    return( as.vector( t(alphastar) %*% temp_vec - mu_smooth*sum(alphastar^2)/2 ) )
  }
}

#helper function for the gradient of the nesterov smoothed approximation of the lasso
# (no extension to general penalties, for which we specifically use the fused version below)
smoothed_lasso_prime_internal <- function(para, lambda, mu_smooth, penweights){

  if(is.null(mu_smooth)){
    return(numeric(length(para)))
  }
  if(mu_smooth==0){
    return(numeric(length(para)))
  }

  if(!is.null(penweights) && length(penweights)==length(para)){
    temp_vec <- abs(para) * lambda * penweights / mu_smooth
  } else{
    temp_vec <- abs(para) * lambda / mu_smooth
  }
  alphastar <-  sign(temp_vec)*pmin(1, abs(temp_vec))

  return( as.vector( alphastar ) )
}

# function to compute nesterov smoothed parameterwise lasso term for penalized (negative) log likelihood
# if mu_smooth = 0, this should return the true lasso penalty (with any associated w)
# this computes it on the 'mean' scale, i.e., the penalty is not scaled by the number of observations, which can happen in the calling function
smoothed_lasso_func <- function(para,nP1,nP2,nP3, lambda, mu_smooth, penweights_list){
  #
  # if(is.null(mu_smooth)){
  #   return(0)
  # }
  # if(mu_smooth==0){
  #   return(0)
  # }

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

  ##COMPUTE PENALTY##
  ##***************##
  pen <- 0
  #break out the beta vectors from the larger parameter vector, using nP0 to correctly pad out the baseline hazard and theta parameters
  if(nP1 != 0){
    # beta1 <- para[(1+nP0):(nP0+nP1)]
    pen <- pen + sum(smoothed_lasso_internal(para=para[(1+nP0):(nP0+nP1)], lambda=lambda1,
                                             mu_smooth=mu_smooth, penweights=penweights_list[["coef1"]]))
  }
  if(nP2 != 0){
    # beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    pen <- pen + sum(smoothed_lasso_internal(para=para[(1+nP0+nP1):(nP0+nP1+nP2)], lambda=lambda2,
                                             mu_smooth=mu_smooth, penweights=penweights_list[["coef2"]]))
  }
  if(nP3 != 0){
    # beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    pen <- pen + sum(smoothed_lasso_internal(para=para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)], lambda=lambda3,
                                             mu_smooth=mu_smooth, penweights=penweights_list[["coef3"]]))
  }
  return(pen)
}



# function to compute nesterov smoothed parameterwise lasso term for penalized (negative) gradient
# if mu_smooth = 0, this should return 0
# this computes it on the 'mean' scale, i.e., the penalty is not scaled by the number of observations, which can happen in the calling function
smoothed_lasso_prime_func <- function(para,nP1,nP2,nP3,lambda,mu_smooth,penweights_list){

  #redefine lambda and lambda_fusedcoef depending on a single value or three separate values are given
  if(length(lambda)==1){
    lambda1 <- lambda2 <- lambda3 <- lambda
  } else if(length(lambda)==3){
    lambda1 <- lambda[1]; lambda2 <- lambda[2]; lambda3 <- lambda[3]
  } else{ stop("lambda is neither a single value or a 3-vector!!") }

  nPtot <- length(para)
  nP0 <- length(para) - nP1 - nP2 - nP3

  grad_part <- numeric(nPtot)

  #break out the beta vectors from the larger parameter vector, using nP0 to correctly pad out the baseline hazard and theta parameters
  if(nP1 != 0){
    beta1 <- para[(1+nP0):(nP0+nP1)]
    grad_part[(1+nP0):(nP0+nP1)] <- smoothed_lasso_prime_internal(para=beta1, lambda=lambda1,
                                                                  mu_smooth=mu_smooth,penweights=penweights_list[["coef1"]])
  }
  if(nP2 != 0){
    beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    grad_part[(1+nP0+nP1):(nP0+nP1+nP2)] <- smoothed_lasso_prime_internal(para=beta2, lambda=lambda2,
                                                                          mu_smooth=mu_smooth,penweights=penweights_list[["coef2"]])
  }
  if(nP3 != 0){
    beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    grad_part[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)] <- smoothed_lasso_prime_internal(para=beta3, lambda=lambda3,
                                                                                  mu_smooth=mu_smooth,penweights=penweights_list[["coef3"]])
  }

  # names(grad_part) <- names(para)
  return(grad_part)
}



#function for the nesterov smoothed approximation of the fused lasso
#Defined in Chen (2012) equation 3.5, with the solution alphastar defined
#in proposition 2
smoothed_fused_lasso_func <- function(para,pen_mat_w_lambda, mu_smooth_fused){
  if(is.null(pen_mat_w_lambda) | is.null(mu_smooth_fused)){
    return(0)
  }
  if(mu_smooth_fused==0){
    return(norm(pen_mat_w_lambda %*% para, type = "1"))
  }
  temp_vec <- pen_mat_w_lambda %*% para
  alphastar <-  sign(temp_vec/mu_smooth_fused)*pmin(1, abs(temp_vec/mu_smooth_fused))

  return( as.vector( t(alphastar) %*% temp_vec - mu_smooth_fused*sum(alphastar^2)/2 ) )
}


#function for the gradient of the nesterov smoothed approximation of the fused lasso
smoothed_fused_lasso_prime_func <- function(para,pen_mat_w_lambda, mu_smooth_fused){
  if(is.null(pen_mat_w_lambda) | is.null(mu_smooth_fused)){
    return(0)
  }
  if(mu_smooth_fused==0){
    return(0)
  }
  temp_vec <- pen_mat_w_lambda %*% para /mu_smooth_fused
  alphastar <-  sign(temp_vec)*pmin(1, abs(temp_vec))

  return( as.vector( t(pen_mat_w_lambda) %*% alphastar ) )
}

