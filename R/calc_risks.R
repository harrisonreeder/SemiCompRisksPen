calc_risk_WB <- function(para, X1, X2, X3,
                         t_cutoff, tol=1e-3, frailty=TRUE,
                         model = "semi-markov"){
  #notice reduced default tolerance
  # browser()
  ##TO START, EXTRACT POINT ESTIMATES OF ALL PARAMETERS FROM MODEL OBJECT##
  ##*********************************************************************##

  nP0 <- if(frailty) 7 else 6

  if(!is.null(X1) && !(ncol(X1)==0)){
    nP1 <- ncol(X1)
    beta1 <- para[(1+nP0):(nP0+nP1)]
    eta1 <- as.vector(X1 %*% beta1)
  } else{
    nP1 <- 0
    eta1 <- 0
  }
  if(!is.null(X2) && !(ncol(X2)==0)){
    nP2 <- ncol(X2)
    beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    eta2 <- as.vector(X2 %*% beta2)
  } else{
    nP2 <- 0
    eta2 <- 0
  }
  if(!is.null(X3) && !(ncol(X3)==0)){
    nP3 <- ncol(X3)
    beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    eta3 <- as.vector(X3 %*% beta3)
  } else{
    nP3 <- 0
    eta3 <- 0
  }

  n <- max(1,nrow(X1),nrow(X2),nrow(X3))

  stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
  alpha1 = exp(para[2])
  alpha2 = exp(para[4])
  alpha3 = exp(para[6])
  kappa1 = exp(para[1])
  kappa2 = exp(para[3])
  kappa3 = exp(para[5])

  ##NOW, USE THE ESTIMATES TO COMPUTE THE PREDICTIONS##
  ##*************************************************##

  #first, compute some helper quantities
  #note this might crap out if there are no covariates in an arm?
    h1_const = alpha1 * kappa1 * exp(eta1)
    h2_const = alpha2 * kappa2 * exp(eta2)
    h3_const = alpha3 * kappa3 * exp(eta3)
    H1_const = kappa1 * exp(as.vector(eta1))
    H2_const = kappa2 * exp(as.vector(eta2))
    H3_const = kappa3 * exp(as.vector(eta3))
    alpha1_m1 = alpha1 - 1
    alpha2_m1 = alpha2 - 1
    alpha3_m1 = alpha3 - 1

  ##*****************************************************************##
  ## Calculating posterior predictive density for WB baseline hazard ##
  ##*****************************************************************##

  #First, write functions that compute the integrand, which we feed into integration function
  #the univariate function when T1=infinity
  f_t2 <- function(t2, index){
    h2_const[index] * (t2)^alpha2_m1 * exp(-H1_const[index] * (t2)^alpha1 - H2_const[index] * (t2)^alpha2)
  }

  #next, the different regions of the joint density on the upper triangle
  if(tolower(model) == "semi-markov"){
    #this is the old code for the original 2-dim integral on the upper triangle. We avoid that now.
    # f.joint <- function(time.pt.vec){
    #   ##
    #   time.pt1 = time.pt.vec[1]
    #   time.pt2 = time.pt.vec[2]
    #   if(time.pt1 == 0 | time.pt2 == 0 | (time.pt1 >= time.pt2)){
    #     return(0)
    #   }
    #   return((h1.const * time.pt1^alpha1.m1) * (h3.const * (time.pt2-time.pt1)^alpha3.m1) *
    #            exp(-H1.const * time.pt1^alpha1 - H2.const * time.pt1^alpha2 - H3.const * (time.pt2 - time.pt1)^alpha3))
    # }

    #if we pre-integrate t2 from t_cutoff to infinity
  # \int \lambda_1(t_1)\exp(-\Lambda_1(t_1)-\Lambda_2(t_1)) \left[\exp(-\Lambda_3(u-t_1)) - \exp(-\Lambda_3(v-t_1))\right] dt_1

    f_joint_t1_nonTerm <- function(time_pt1,index){
      return( ( h1_const[index] * time_pt1^alpha1_m1 *
                exp( -H1_const[index] * time_pt1^alpha1 - H2_const[index] * time_pt1^alpha2 ) ) *
                exp( -H3_const[index] * (t_cutoff - time_pt1)^alpha3))
    }
    #if we pre-integrate t2 from t1 to t_cutoff
    f_joint_t1_both <- function(time_pt1,index){
      return( ( h1_const[index] * time_pt1^alpha1_m1 *
                exp( -H1_const[index] * time_pt1^alpha1 - H2_const[index] * time_pt1^alpha2 ) ) *
                ( 1 - exp( -H3_const[index] * (t_cutoff - time_pt1)^alpha3)) )
    }
    #if we pre-integrate t2 from t1 to infinity
    f_joint_t1_neither <- function(time_pt1,index){
      return( ( h1_const[index] * time_pt1^alpha1_m1 *
                exp( -H1_const[index] * time_pt1^alpha1 - H2_const[index] * time_pt1^alpha2 ) ) )
    }

  } else{
    stop("model must be 'semi-markov'")
  }

  p_ntonly <- sapply(1:n,function(x){tryCatch(integrate(f_joint_t1_nonTerm, lower = 0, upper = t_cutoff, index=x)$value,
                                error = function(cnd){
                                  # message(cnd)
                                  # cat("\n")
                                  return(NA)})
    })

  p_tonly <- sapply(1:n,function(x){tryCatch(integrate(f_t2, lower = 0, upper = t_cutoff, index=x)$value,
                                error = function(cnd){
                                  # message(cnd)
                                  # cat("\n")
                                  return(NA)})
    })

  p_neither_t2 <- sapply(1:n,function(x){tryCatch(integrate(f_t2, lower = t_cutoff, upper = Inf,index=x)$value,
                                 error = function(cnd){
                                   # message(cnd)
                                   # cat("\n")
                                   return(NA)})
    })

  p_neither_joint <-  sapply(1:n,function(x){tryCatch(integrate(f_joint_t1_neither, lower = t_cutoff, upper = Inf,index=x)$value,
                               error = function(cnd){
                                 # message(cnd)
                                 # cat("\n")
                                 return(NA)})
    })

  p_both <- sapply(1:n,function(x){tryCatch(integrate(f_joint_t1_both, lower = 0, upper = t_cutoff, index=x)$value,
                     error = function(cnd){
                       # message(cnd)
                       # cat("\n")
                       return(NA)})
    })
  # end <- Sys.time()
  # print(end-begin)
  # print("did both part of joint")
  # return(p_ntonly + p_both)
  return(cbind(p_ntonly = p_ntonly,
              p_both = p_both,
              p_tonly = p_tonly,
              p_neither = p_neither_t2 + p_neither_joint))
}



calc_risk_PW <- function(para, X1, X2, X3, knots_list,
                         t_cutoff, tol=1e-3, frailty=TRUE,
                         model = "semi-markov"){
  #notice reduced default tolerance
  # browser()
  ##TO START, EXTRACT POINT ESTIMATES OF ALL PARAMETERS FROM MODEL OBJECT##
  ##*********************************************************************##

  #left pad with a zero if it is not already present
  if(knots_list[[1]][1] != 0){knots_list[[1]] <- c(0,knots_list[[1]])}
  if(knots_list[[2]][1] != 0){knots_list[[2]] <- c(0,knots_list[[2]])}
  if(knots_list[[3]][1] != 0){knots_list[[3]] <- c(0,knots_list[[3]])}

  nP01 <- length(knots_list[[1]])
  nP02 <- length(knots_list[[2]])
  nP03 <- length(knots_list[[3]])
  nP0 <- if(frailty) nP01 + nP02 + nP03 + 1 else nP01 + nP02 +nP03


  if(!is.null(X1) && !(ncol(X1)==0)){
    nP1 <- ncol(X1)
    beta1 <- para[(1+nP0):(nP0+nP1)]
    eta1 <- as.vector(X1 %*% beta1)
  } else{
    nP1 <- 0
    eta1 <- 0
  }
  if(!is.null(X2) && !(ncol(X2)==0)){
    nP2 <- ncol(X2)
    beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    eta2 <- as.vector(X2 %*% beta2)
  } else{
    nP2 <- 0
    eta2 <- 0
  }
  if(!is.null(X3) && !(ncol(X3)==0)){
    nP3 <- ncol(X3)
    beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    eta3 <- as.vector(X3 %*% beta3)
  } else{
    nP3 <- 0
    eta3 <- 0
  }


  n <- max(nrow(X1),nrow(X2),nrow(X3))

  stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss

  phi1 <- as.numeric(para[(1):(nP01)])
  phi2 <- as.numeric(para[(1+nP01):(nP01+nP02)])
  phi3 <- as.numeric(para[(1+nP01+nP02):(nP01+nP02+nP03)])

  ##NOW, USE THE ESTIMATES TO COMPUTE THE PREDICTIONS##
  ##*************************************************##

  haz <- function(t,phi,knots){
    exp(phi)[findInterval(x = t, vec = knots, left.open = TRUE)]
  }

  Haz <- function(t,phi,knots){
    rowSums(sweep(x = pw_cum_mat(t,knots),MARGIN = 2,STATS=exp(phi),FUN ="*"))
  }

  ##*****************************************************************##
  ## Calculating posterior predictive density for WB baseline hazard ##
  ##*****************************************************************##


  #First, write functions that compute the integrand, which we feed into integration function

  #first, the univariate function when T1=infinity
  f_t2 <- function(t2,index){
    haz(t=t2,phi=phi2,knots = knots_list[[2]]) * exp(eta2[index]) *
      exp(-Haz(t=t2,phi=phi1,knots=knots_list[[1]])* exp(eta1[index]) -
            Haz(t=t2,phi=phi2,knots=knots_list[[2]])* exp(eta2[index]))
  }


  #next, the different regions of the joint density on the upper triangle
  if(tolower(model) == "semi-markov"){
    #this is the old code for the original 2-dim integral on the upper triangle. We avoid that now.
    # f.joint <- function(time.pt.vec){
    #   ##
    #   time.pt1 = time.pt.vec[1]
    #   time.pt2 = time.pt.vec[2]
    #   if(time.pt1 == 0 | time.pt2 == 0 | (time.pt1 >= time.pt2)){
    #     return(0)
    #   }
    #   return((h1.const * time.pt1^alpha1.m1) * (h3.const * (time.pt2-time.pt1)^alpha3.m1) *
    #            exp(-H1.const * time.pt1^alpha1 - H2.const * time.pt1^alpha2 - H3.const * (time.pt2 - time.pt1)^alpha3))
    # }

    #if we pre-integrate t2 from t_cutoff to infinity
    f_joint_t1_nonTerm <- function(time_pt1,index){
      return( ( haz(t=time_pt1,phi=phi1,knots = knots_list[[1]]) * exp(eta1[index]) *
                  exp(-Haz(t=time_pt1,phi=phi1,knots = knots_list[[1]]) * exp(eta1[index]) -
                        Haz(t=time_pt1,phi=phi2,knots = knots_list[[2]])* exp(eta2[index]) ) ) *
                exp( -Haz(t=t_cutoff-time_pt1,phi=phi3,knots = knots_list[[3]]) * exp(eta3[index])))
    }
    #if we pre-integrate t2 from t1 to t_cutoff
    f_joint_t1_both <- function(time_pt1,index){
      return( ( haz(t=time_pt1,phi=phi1,knots = knots_list[[1]]) * exp(eta1[index]) *
                  exp(-Haz(t=time_pt1,phi=phi1,knots = knots_list[[1]]) * exp(eta1[index]) -
                        Haz(t=time_pt1,phi=phi2,knots = knots_list[[2]])* exp(eta2[index]) ) ) *
                ( 1 - exp( -Haz(t=t_cutoff-time_pt1,phi=phi3,knots = knots_list[[3]]) * exp(eta3[index]))) )
    }
    #if we pre-integrate t2 from t1 to infinity
    f_joint_t1_neither <- function(time_pt1,index){
      return( ( haz(t=time_pt1,phi=phi1,knots = knots_list[[1]]) * exp(eta1[index]) *
                  exp(-Haz(t=time_pt1,phi=phi1,knots = knots_list[[1]]) * exp(eta1[index]) -
                        Haz(t=time_pt1,phi=phi2,knots = knots_list[[2]])* exp(eta2[index]) ) ) )
    }

  } else{
    stop("model must be 'semi-markov'")
  }

  p_ntonly <- sapply(1:n,function(x){tryCatch(integrate(f_joint_t1_nonTerm, lower = 0, upper = t_cutoff, index=x)$value,
                                              error = function(cnd){
                                                # message(cnd)
                                                # cat("\n")
                                                return(NA)})
  })

  p_tonly <- sapply(1:n,function(x){tryCatch(integrate(f_t2, lower = 0, upper = t_cutoff, index=x)$value,
                                             error = function(cnd){
                                               # message(cnd)
                                               # cat("\n")
                                               return(NA)})
  })

  p_neither_t2 <- sapply(1:n,function(x){tryCatch(integrate(f_t2, lower = t_cutoff, upper = Inf,index=x)$value,
                                                  error = function(cnd){
                                                    # message(cnd)
                                                    # cat("\n")
                                                    return(NA)})
  })

  p_neither_joint <-  sapply(1:n,function(x){tryCatch(integrate(f_joint_t1_neither, lower = t_cutoff, upper = Inf,index=x)$value,
                                                      error = function(cnd){
                                                        # message(cnd)
                                                        # cat("\n")
                                                        return(NA)})
  })

  p_both <- sapply(1:n,function(x){tryCatch(integrate(f_joint_t1_both, lower = 0, upper = t_cutoff, index=x)$value,
                                            error = function(cnd){
                                              # message(cnd)
                                              # cat("\n")
                                              return(NA)})
  })
  # end <- Sys.time()
  # print(end-begin)
  # print("did both part of joint")
  # return(p_ntonly + p_both)
  return(cbind(p_ntonly = p_ntonly,
               p_both = p_both,
               p_tonly = p_tonly,
               p_neither = p_neither_t2 + p_neither_joint))
}
