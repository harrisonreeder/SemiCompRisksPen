##******************************************##
####Define Regularization Helper Functions####
##******************************************##


penweights_internal <- function(parahat,D,penweight_type="adaptive",addl_penweight){
  ####function to compute adaptive weights based on input
  ####written to be general enough to apply to either univariate or fused penalty
  if(penweight_type == "adaptive"){
    if(!is.null(D)){
      if(!is.null(parahat)){
        ada_weight <- abs(D %*% parahat)#^(-a) #originally, the adaptive lasso allows the weight to be scaled by an exponent, but that causes confusion I feel.
        if(any(!is.finite(ada_weight))){
          warning("at least one adaptive fusion weight is not finite, setting to 0") #wait, it may not make sense to set to 0
          ada_weight[!is.finite(ada_weight)] <- 0
        }
      } else {
        ada_weight <- rep(1,nrow(D))
      }
    } else if(!is.null(parahat)){
      ada_weight <- abs(parahat)#^(-a) #originally, the adaptive lasso allows the weight to be scaled by an exponent, but that causes confusion I feel.
    } else {
      ada_weight <- rep(1,length(parahat))
    }

    if(!is.null(addl_penweight) ){
      if(length(ada_weight)==length(addl_penweight)){
        ada_weight <- ada_weight * addl_penweight
      } else{
        stop("supplied additional weight vector is incorrect dimension")
      }
    }
    return(ada_weight)

  } else{ #if not adaptive weights, then just return whatever additional weights might exist
    return(addl_penweight)
  }
}

pen_internal <- function(para, penalty, lambda, a, D, penweights){
  #unified function that returns penalties--  ASSUMES CHECKS HAVE BEEN DONE
  #D is a matrix of contrasts, with number of columns equal to the total number of parameters, rows equal to the total number of differences being `fused'
  #e.g., to fuse first two transitions' betas, set D <- cbind(matrix(data=0,nrow=nP1,ncol=7), diag(nP1), -diag(nP2), matrix(data=0,nrow=nP1,ncol=nP1))

  if(is.null(D)){
    beta <- abs(para)
  } else{
    beta <- as.vector(abs(D %*% para))
  }

  if(penalty == "lasso"){
    out <- abs(beta) * lambda
  } else if(penalty == "scad"){
    beta <- abs(beta)
    ind1 <- ifelse(beta <= lambda,1,0)
    ind2 <- ifelse(beta > lambda & beta <= a*lambda ,1,0)
    ind3 <- ifelse(beta > a*lambda ,1,0)
    out <- ind1 * lambda*beta +
      ind2 * (2*a*lambda*beta - beta^2 - lambda^2)/(2*(a-1)) +
      ind3 * lambda^2*(a+1)/2
  } else if(penalty == "mcp"){
    beta <- abs(beta)
    ind1 <- rep(0, length(beta))
    ind1[beta <= a*lambda] <- 1
    out <- ind1 * (lambda*beta - beta^2/(2*a)) +
      (1-ind1)* (a*lambda^2)/2
  } else{
    out <- rep(0, length(beta)) #vector of 0's
  }

  #in principle, any penalty can have weights added, e.g., to set some terms to 0
  if(!is.null(penweights) && length(penweights)==length(beta)){
    return(penweights * out)
  } else{
    # warning("weights supplied to lasso function are NULL, or of incorrect length. ignoring weights.")
    return(out)
  }
}

pen_prime_internal <- function(beta, penalty, lambda, a, penweights){

  ##Function implementing the 'derivative' of various penalty terms (defined everywhere except where beta=0)
  #beta and lambda are in all of them, a is in SCAD, MCP and adalasso, and then betahat is in adalasso

  if(penalty == "lasso"){
    out <- lambda
  } else if(penalty == "scad"){
    beta <- abs(beta)
    ind2 <- ind1 <- rep(0, length(beta))
    ind1[beta > lambda] <- 1
    ind2[beta <= (lambda * a)] <- 1
    out <- lambda * (1 - ind1) +
      ((lambda * a) - beta) * ind2/(a - 1) * ind1
  } else if(penalty == "mcp"){ #careful of the sign here! This is from breheny 2-29 slide 15
    ind1 <- rep(0, length(beta))
    ind1[abs(beta) <= a*lambda] <- 1
    # out <- ind1*(lambda - abs(beta)/a)*sign(beta)
    out <- ind1*(lambda - abs(beta)/a)
  } else{
    out <- rep(0, length(beta)) #vector of 0's
  }

  #in principle, any penalty can have weights added, e.g., to set some terms to 0
  if(!is.null(penweights) && length(penweights)==length(beta)){
    return(penweights * out)
  } else{
    # warning("weights supplied to lasso function are NULL, or of incorrect length. ignoring weights.")
    return(out)
  }
}

##***********************************##
####Define Regularization Functions####
##***********************************##

pen_func <- function(para,nP1,nP2,nP3,
                     penalty, lambda, a,
                     penweights_list,
                     pen_mat_w_lambda){
  #function to compute complete penalty term for penalized (negative) log likelihood
  #this computes it on the 'mean' scale, i.e., the penalty is not scaled by the number of observations, which can happen in the calling function

  #checks on 'a' to ensure no errors
  if(is.null(a)){
    a <- switch(tolower(penalty),"scad"=3.7, "mcp"=3, "adalasso"=1, "lasso"=1)
  }

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
    pen <- pen + sum(pen_internal(para=para[(1+nP0):(nP0+nP1)],
                                  penalty=penalty,lambda=lambda1,a=a,D=NULL,penweights=penweights_list[["coef1"]]))
  }
  if(nP2 != 0){
    # beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    pen <- pen + sum(pen_internal(para=para[(1+nP0+nP1):(nP0+nP1+nP2)],
                                  penalty=penalty,lambda=lambda2,a=a,D=NULL,penweights=penweights_list[["coef2"]]))
  }
  if(nP3 != 0){
    # beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    pen <- pen + sum(pen_internal(para=para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)],
                                  penalty=penalty,lambda=lambda3,a=a,D=NULL,penweights=penweights_list[["coef3"]]))
  }
  ##If fused lasso is being used, then the following comes into play##
  ##THIS MATCHES THE \|C\bbeta\| formulation of the fused penalty seen in Chen (2012)##
  if(!is.null(pen_mat_w_lambda)){
    pen <- pen + norm(pen_mat_w_lambda %*% para, type = "1")
  }

  return(pen)
}


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
