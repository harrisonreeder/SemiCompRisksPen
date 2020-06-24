####MM APPROACH####

##Function implementing majorization of various penalty terms (based on a perturbed quadratic approximation around a particular beta0)
#these formulae follow from hunter and li p. 1624, and are useful for plotting
pen_maj_internal <- function(beta,beta0,penalty,lambda,a,penweights,mm_epsilon){

  beta0_pen <- pen_internal(para=beta0,penalty=penalty,lambda=lambda,a=a,D=NULL,penweights=penweights)
  beta0_pen_prime <- pen_prime_internal(para=beta0,penalty=penalty,lambda=lambda,a=a,D=NULL,penweights=penweights)

  out <- beta0_pen +
    (beta^2 - beta0^2) * beta0_pen_prime /
    #              (2*abs(mm_epsilon+abs(beta0))) #Hunter & Li's Formulation of the approximation
    (2*sqrt(mm_epsilon+beta0^2)) #Oelker & Tutz's Formulation of the approximation

  #in principle, any penalty can have weights added, e.g., to set some terms to 0
  if(!is.null(penweights) && length(penweights)==length(beta)){
    return(penweights * out)
  } else{
    # warning("weights supplied to lasso function are NULL, or of incorrect length. ignoring weights.")
    return(out)
  }
}

##Helper function for majorized gradients and hessians in fused lasso
##this specification follows from Sennhenn-Reuler & Kneib (2015) p.4645 (equation for P_lambda)
##as well as Oelker & Tutz (2017) p. 105 (equation for A_lambda)
##In short, this function creates the fused lasso-related components of E_k (Hunter & Li), P_lambda (Sennhenn-Reulen), or A_lambda (Oelker & Tutz)
lasso_prime_mat_internal <- function(para,lambda,D,penweights,mm_epsilon){
  #D is a matrix of contrasts, with number of columns equal to the total number of parameters, rows equal to the total number of differences being `fused'
  #e.g., to fuse first two transitions' betas, set D <- cbind(matrix(data=0,nrow=nP1,ncol=7), diag(nP1), -diag(nP2), matrix(data=0,nrow=nP1,ncol=nP1))

  # temp_weight <- lambda/(abs(D %*% para) + mm_epsilon) #using Hunter & Li's formulation, rather than Sennhenn-Reulen's
  temp_weights <- lambda/sqrt((D %*% para)^2 + mm_epsilon) #This would be Sennhenn-Reulen's / Oelker & Tutz's formulation.

  #in principle, any penalty can have weights added, e.g., to set some terms to 0
  if(!is.null(penweights) && length(penweights)==length(temp_weights)){
    temp_weights <- penweights * temp_weights
  }
  return( t(D) %*% diag(as.vector(temp_weights)) %*% D )
}

#function to compute the majorized penalty matrix used to compute majorized gradient and hessian
#this is also on the 'mean' rather than the 'sum' scale, so must be scaled appropriately if grad/hess depend on sample size
pen_maj_mat_func <- function(para,nP1,nP2,nP3,
                             penalty, lambda, a,
                             penalty_fusedcoef, lambda_fusedcoef,
                             penalty_fusedbaseline, lambda_fusedbaseline,
                             penweights_list, mm_epsilon, hazard){

  #checks on 'a' to ensure no errors
  if(is.null(a)){
    a <- switch(tolower(penalty),"scad"=3.7, "mcp"=3, "lasso"=1)
  }

  check_pen_params(penalty,penalty_fusedcoef,penalty_fusedbaseline,
                   lambda,lambda_fusedcoef,lambda_fusedbaseline,a)

  #redefine lambda and lambda_fusedcoef depending on a single value or three separate values are given
  if(length(lambda)==1){
    lambda1 <- lambda2 <- lambda3 <- lambda
  } else if(length(lambda)==3){
    lambda1 <- lambda[1]; lambda2 <- lambda[2]; lambda3 <- lambda[3]
  } else{ stop("lambda is neither a single value or a 3-vector!!") }

  nPtot <- length(para)
  nP0 <- length(para) - nP1 - nP2 - nP3

  Ek_vec <- numeric(nPtot) #start with a vector of 0's

  #break out the beta vectors from the larger parameter vector, using nP0 to correctly pad out the baseline hazard and theta parameters
  if(nP1 != 0){
    beta1 <- para[(1+nP0):(nP0+nP1)]
    Ek_vec[(1+nP0):(nP0+nP1)] <- pen_prime_internal(beta=beta1,lambda=lambda1,penalty=penalty,a=a,penweights=penweights_list[["coef1"]]) /
      (sqrt(mm_epsilon+beta1^2)) #Oelker & Tutz's Formulation of the approximation
    #               (mm_epsilon+abs(beta1)), #Hunter & Li's Formulation of the approximation

  }
  if(nP2 != 0){
    beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    Ek_vec[(1+nP0+nP1):(nP0+nP1+nP2)] <- pen_prime_internal(beta=beta2,lambda=lambda2,penalty=penalty,a=a,penweights=penweights_list[["coef2"]]) /
      (sqrt(mm_epsilon+beta2^2)) #Oelker & Tutz's Formulation of the approximation
    #               (mm_epsilon+abs(beta2)), #Hunter & Li's Formulation of the approximation
  }
  if(nP3 != 0){
    beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    Ek_vec[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)] <- pen_prime_internal(beta=beta3,lambda=lambda3,penalty=penalty,a=a,penweights=penweights_list[["coef3"]]) /
      (sqrt(mm_epsilon+beta3^2)) #Oelker & Tutz's Formulation of the approximation
    #               (mm_epsilon+abs(beta3)) #Hunter & Li's Formulation of the approximation
  }

  Ek <- diag(Ek_vec)
  # print(Ek)

  ##Then, if fused lasso is being used, then the following comes into play##
  if(tolower(penalty_fusedcoef) %in% c("fusedlasso","adafusedlasso")){

    if(length(lambda_fusedcoef)==1){
      lambda_fusedcoef12 <- lambda_fusedcoef13 <- lambda_fusedcoef23 <- lambda_fusedcoef
    } else if(length(lambda_fusedcoef)==3){
      lambda_fusedcoef12 <- lambda_fusedcoef[1]; lambda_fusedcoef13 <- lambda_fusedcoef[2];lambda_fusedcoef23 <- lambda_fusedcoef[3]
    } else{ stop("lambda_fusedcoef is neither a single value or a 3-vector!!") }

    #create a difference matrix connecting first and second sets of betas, left padded by 7 columns of zeros for the hazards
    if(lambda_fusedcoef12 != 0){
      stopifnot(nP1==nP2)
      D_temp <- cbind(matrix(data=0,nrow=nP1,ncol=nP0), diag(nP1), -diag(nP1), matrix(data=0,nrow=nP1,ncol=nP1))
      Ek <- Ek + lasso_prime_mat_internal(para=para,lambda=lambda_fusedcoef12,D=D_temp,penweights=penweights_list[["fusedcoef12"]],mm_epsilon=mm_epsilon)
    }
    #create a difference matrix connecting first and third sets of betas, left padded by 7 columns of zeros for the hazards
    if(lambda_fusedcoef13 != 0){
      stopifnot(nP1==nP3)
      D_temp <- cbind(matrix(data=0,nrow=nP1,ncol=nP0), diag(nP1),matrix(data=0,nrow=nP1,ncol=nP1),-diag(nP1))
      Ek <- Ek + lasso_prime_mat_internal(para=para,lambda=lambda_fusedcoef13,D=D_temp,penweights=penweights_list[["fusedcoef13"]],mm_epsilon=mm_epsilon)
    }
    #create a difference matrix connecting second and third sets of betas, left padded by 7 columns of zeros for the hazards
    if(lambda_fusedcoef23 != 0){
      stopifnot(nP2==nP3)
      D_temp <- cbind(matrix(data=0,nrow=nP2,ncol=nP0), matrix(data=0,nrow=nP2,ncol=nP2),diag(nP2),-diag(nP2))
      Ek <- Ek + lasso_prime_mat_internal(para=para,lambda=lambda_fusedcoef23,D=D_temp,penweights=penweights_list[["fusedcoef23"]],mm_epsilon=mm_epsilon)
    }
  }

  #Finally, add penalty of fusion for baseline parameters
  if(tolower(penalty_fusedbaseline) != "none" & !(is.null(hazard))){
    if(hazard != "weibull"){
      stop("non-Weibull not yet implemented")
    }

    if(length(lambda_fusedbaseline)==1){
      lambda_fusedbaseline12 <- lambda_fusedbaseline13 <- lambda_fusedbaseline23 <- lambda_fusedbaseline
    } else if(length(lambda_fusedbaseline)==3){
      lambda_fusedbaseline12 <- lambda_fusedbaseline[1]
      lambda_fusedbaseline13 <- lambda_fusedbaseline[2]
      lambda_fusedbaseline23 <- lambda_fusedbaseline[3]
    } else{ stop("lambda_fusedbaseline is neither a single value or a 3-vector!!") }

    #create a difference matrix connecting first and second sets of baseline parameters
    if(lambda_fusedbaseline12 != 0){
      D_temp <- matrix(data=0,nrow=2,ncol=nPtot)
      D_temp[1,1] <- D_temp[2,2] <- 1
      D_temp[1,3] <- D_temp[2,4] <- -1
      Ek <- Ek + lasso_prime_mat_internal(para=para,lambda=lambda_fusedbaseline12,D=D_temp,penweights=penweights_list[["fusedbaseline12"]],mm_epsilon=mm_epsilon)
    }
    #create a difference matrix connecting first and third sets of baseline parameters
    if(lambda_fusedbaseline13 != 0){
      D_temp <- matrix(data=0,nrow=2,ncol=nPtot)
      D_temp[1,1] <- D_temp[2,2] <- 1
      D_temp[1,5] <- D_temp[2,6] <- -1
      Ek <- Ek + lasso_prime_mat_internal(para=para,lambda=lambda_fusedbaseline13,D=D_temp,penweights=penweights_list[["fusedbaseline13"]],mm_epsilon=mm_epsilon)
    }
    #create a difference matrix connecting second and third sets of baseline parameters
    if(lambda_fusedbaseline23 != 0){
      D_temp <- matrix(data=0,nrow=2,ncol=nPtot)
      D_temp[1,3] <- D_temp[2,4] <- 1
      D_temp[1,5] <- D_temp[2,6] <- -1
      Ek <- Ek + lasso_prime_mat_internal(para=para,lambda=lambda_fusedbaseline23,D=D_temp,penweights=penweights_list[["fusedbaseline23"]],mm_epsilon=mm_epsilon)
    }
  }

  return(Ek)
}

#function to compute the majorized penalty where para0 is the point of expansion, and para is the current point
#this computes it on the 'mean' scale, i.e., the penalty is not scaled by the number of observations, which can happen in the calling function
maj_func <- function(para, para0, nP1, nP2, nP3,
                     penalty, lambda, a,
                     penalty_fusedcoef, lambda_fusedcoef,
                     penalty_fusedbaseline, lambda_fusedbaseline, pen_mat_w_lambda,
                     penweights_list, mm_epsilon, hazard){
  # browser()
  #penalty function assessed at the null point
  pen0 <- pen_func(para=para0,nP1=nP1,nP2=nP2,nP3=nP3,
                   penalty=penalty,lambda=lambda, a=a,
                   penweights_list=penweights_list, pen_mat_w_lambda = pen_mat_w_lambda)

  #majorized penalty matrix assessed at the null point
  Ek0 <- pen_maj_mat_func(para=para0,nP1=nP1,nP2=nP2,nP3=nP3,
                          penalty=penalty,lambda=lambda, a=a,mm_epsilon=mm_epsilon,
                          penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                          penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                          penweights_list=penweights_list, hazard=hazard)

  #formula from Oelker & Tutz (2017) p. 105 formula for A_lambda, except I am fixing what I presume is a typo
  #where there should be a subtraction instead of an addition between the two quadratic forms
  maj <- pen0 + ( t(para) %*% Ek0 %*% para - t(para0) %*% Ek0 %*% para0 ) / 2

  return(as.numeric(maj))
}


