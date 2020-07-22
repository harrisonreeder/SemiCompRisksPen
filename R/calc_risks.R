#' Calculate absolute risk profiles
#'
#' This function calculates absolute risk profiles under weibull baseline hazard specifications.
#'
#' @inheritParams proximal_gradient_descent
#' @param t_cutoff Numeric vector indicating the time(s) to compute the risk profile.
#' @param tol Numeric value for the tolerance of the numerical integration procedure.
#' @param type String either indicating 'marginal' for population-averaged probabilities,
#'   or 'conditional' for probabilities computed at the specified gamma
#' @param gamma Numeric value indicating the fixed level of the frailty assumed for predicted probabilities,
#'   if 'type' is set to 'conditional'
#'
#' @return if Xmat has only one row, and t_cutoff is a scalar, then returns a 4 element row matrix
#'   of probabilities. If Xmat has \code{n} rows, then returns an \code{n} by 4 matrix of probabilities.
#'   If Xmat has \code{n} rows and t_cutoff is a vector of length \code{s}, then returns an \code{s} by 4 by \code{n} array.
#' @export
calc_risk_WB <- function(para, Xmat1, Xmat2, Xmat3,
                         t_cutoff, tol=1e-3, frailty=TRUE,
                         type="marginal", gamma=1,
                         model="semi-markov"){
  #notice reduced default tolerance
  # browser()
  ##TO START, EXTRACT POINT ESTIMATES OF ALL PARAMETERS FROM MODEL OBJECT##
  ##*********************************************************************##

  n <- max(1,nrow(Xmat1),nrow(Xmat2),nrow(Xmat3))
  t_length <- length(t_cutoff)

  if(frailty){
    nP0 <- 7
    theta <- exp(para[7])
    if(type=="conditional" & length(gamma)==1){
      gamma <- rep(gamma,n)
    }
  } else{
    nP0 <- 6
    type <- "conditional"
    gamma <- rep(1,n)
  }

  if(!is.null(Xmat1) && !(ncol(Xmat1)==0)){
    nP1 <- ncol(Xmat1)
    beta1 <- para[(1+nP0):(nP0+nP1)]
    eta1 <- as.vector(Xmat1 %*% beta1)
  } else{
    nP1 <- 0
    eta1 <- 0
  }
  if(!is.null(Xmat2) && !(ncol(Xmat2)==0)){
    nP2 <- ncol(Xmat2)
    beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    eta2 <- as.vector(Xmat2 %*% beta2)
  } else{
    nP2 <- 0
    eta2 <- 0
  }
  if(!is.null(Xmat3) && !(ncol(Xmat3)==0)){
    nP3 <- ncol(Xmat3)
    beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    eta3 <- as.vector(Xmat3 %*% beta3)
  } else{
    nP3 <- 0
    eta3 <- 0
  }


  stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
  alpha1=exp(para[2])
  alpha2=exp(para[4])
  alpha3=exp(para[6])
  kappa1=exp(para[1])
  kappa2=exp(para[3])
  kappa3=exp(para[5])

  ##NOW, USE THE ESTIMATES TO COMPUTE THE PREDICTIONS##
  ##*************************************************##

  #first, compute some helper quantities
  #note this might crap out if there are no covariates in an arm?
  h1_const=alpha1 * kappa1 * exp(eta1)
  h2_const=alpha2 * kappa2 * exp(eta2)
  h3_const=alpha3 * kappa3 * exp(eta3)
  H1_const=kappa1 * exp(as.vector(eta1))
  H2_const=kappa2 * exp(as.vector(eta2))
  H3_const=kappa3 * exp(as.vector(eta3))
  alpha1_m1=alpha1 - 1
  alpha2_m1=alpha2 - 1
  alpha3_m1=alpha3 - 1

  ##*****************************************************************##
  ## Calculating posterior predictive density for WB baseline hazard ##
  ##*****************************************************************##

  #First, write functions that compute the integrand, which we feed into integration function
  #the univariate function when T1=infinity
  if(tolower(type) %in% c("marginal","m")){
    f_t2 <- function(t2, index){
      h2_const[index] * (t2)^alpha2_m1 *
        (1 + theta*(H1_const[index] * (t2)^alpha1 +
                      H2_const[index] * (t2)^alpha2) )^(-theta^(-1) - 1)
    }
  } else{
    f_t2 <- function(t2, index){
      gamma[index] * h2_const[index] * (t2)^alpha2_m1 *
        exp(-gamma[index]*(H1_const[index] * (t2)^alpha1 +
                             H2_const[index] * (t2)^alpha2))
    }
  }

  #next, the different regions of the joint density on the upper triangle
  # here is the generic formula if we pre-integrate t2 from u to v
  # \int \lambda_1(t_1)\exp(-\Lambda_1(t_1)-\Lambda_2(t_1)) \left[\exp(-\Lambda_3(u-t_1)) - \exp(-\Lambda_3(v-t_1))\right] dt_1
  if(tolower(model) == "semi-markov"){

    #if we pre-integrate t2 from t1 to t_cutoff
    if(tolower(type) %in% c("marginal","m")){
      f_joint_t1_both <- function(time_pt1,t_cutoff,index){
        h1_const[index] * time_pt1^alpha1_m1 * (
          (1 + theta*(H1_const[index] * time_pt1^alpha1 +
                        H2_const[index] * time_pt1^alpha2))^(-theta^(-1)-1) -
            (1 + theta*(H1_const[index] * time_pt1^alpha1 +
                          H2_const[index] * time_pt1^alpha2 +
                          H3_const[index] * (t_cutoff - time_pt1)^alpha3))^(-theta^(-1)-1)
        )}
    } else{
      f_joint_t1_both <- function(time_pt1,t_cutoff,index){
        gamma[index] * h1_const[index] * time_pt1^alpha1_m1 *
                    exp( -gamma[index]*(H1_const[index] * time_pt1^alpha1 + H2_const[index] * time_pt1^alpha2 ) ) *
                  ( 1 - exp( -gamma[index] * H3_const[index] * (t_cutoff - time_pt1)^alpha3))
      }
    }

    #if we pre-integrate t2 from t_cutoff to infinity
    if(tolower(type) %in% c("marginal","m")){
      f_joint_t1_nonTerm <- function(time_pt1,t_cutoff,index){
        h1_const[index] * time_pt1^alpha1_m1 *
                    (1 + theta * (H1_const[index] * time_pt1^alpha1 +
                                           H2_const[index] * time_pt1^alpha2 +
                                           H3_const[index] * (t_cutoff - time_pt1)^alpha3))^(-theta^(-1)-1)
      }
    } else{
      f_joint_t1_nonTerm <- function(time_pt1,t_cutoff,index){
        gamma[index] * h1_const[index] * time_pt1^alpha1_m1 *
                    exp(-gamma[index] * (H1_const[index] * time_pt1^alpha1 +
                                           H2_const[index] * time_pt1^alpha2 +
                                           H3_const[index] * (t_cutoff - time_pt1)^alpha3))
      }
    }

    #if we pre-integrate t2 from t1 to infinity
    if(tolower(type) %in% c("marginal","m")){
      f_joint_t1_neither <- function(time_pt1,index){
        h1_const[index] * time_pt1^alpha1_m1 *
                  (1 + theta*(H1_const[index] * time_pt1^alpha1 +
                                        H2_const[index] * time_pt1^alpha2) )^(-theta^(-1) - 1)
      }
    } else{
      f_joint_t1_neither <- function(time_pt1,index){
        gamma[index]*h1_const[index] * time_pt1^alpha1_m1 *
                    exp( -gamma[index]*(H1_const[index] * time_pt1^alpha1 +
                                          H2_const[index] * time_pt1^alpha2) )
      }
    }

  } else{
    stop("model must be 'semi-markov'")
  }

  if(n > 1){
    if(t_length > 1){
      out_mat <- array(dim=c(t_length,4,n),dimnames = list(paste0("t",t_cutoff),c("p_ntonly","p_both","p_tonly","p_neither"),paste0("i",1:n)))
    } else{
      out_mat <- matrix(nrow=n,ncol=4,dimnames = list(paste0("i",1:n),c("p_ntonly","p_both","p_tonly","p_neither")))
    }
  } else{
      out_mat <- matrix(nrow=t_length,ncol=4,dimnames = list(paste0("t",t_cutoff),c("p_ntonly","p_both","p_tonly","p_neither")))
  }

  for(t_ind in 1:t_length){
    p_ntonly <- sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_nonTerm, lower=0, upper=t_cutoff[t_ind],
                                                          t_cutoff=t_cutoff[t_ind], index=x)$value,
                                                error=function(cnd){return(NA)}) })

    p_tonly <- sapply(1:n,function(x){tryCatch(stats::integrate(f_t2, lower=0, upper=t_cutoff[t_ind], index=x)$value,
                                               error=function(cnd){return(NA)}) })

    p_neither_t2 <- sapply(1:n,function(x){tryCatch(stats::integrate(f_t2, lower=t_cutoff[t_ind], upper=Inf, index=x)$value,
                                                    error=function(cnd){return(NA)}) })

    p_neither_joint <-  sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_neither, lower=t_cutoff[t_ind], upper=Inf, index=x)$value,
                                                        error=function(cnd){return(NA)}) })

    p_both <- sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_both, lower=0, upper=t_cutoff[t_ind],
                                                        t_cutoff=t_cutoff[t_ind], index=x)$value,
                                              error=function(cnd){return(NA)}) })

    out_temp <- cbind(p_ntonly=p_ntonly,
                      p_both=p_both,
                      p_tonly=p_tonly,
                      p_neither=p_neither_t2 + p_neither_joint)

    #I noticed that sometimes, if exactly one category has an NA, then we could back out the value
    #from the other categories. However, we don't want any to be negative.
    #I will also supply a warning if any of the rows are way off from 1.
    out_temp <- t(apply(out_temp,1,
                        function(x){
                          if(sum(is.na(x))==1){
                            x[is.na(x)]<- max(1-sum(x,na.rm=TRUE),0)
                          }
                          return(x)}))
    if(any(is.na(out_temp)) | any(abs(1-rowSums(out_temp))>tol)){
      warning(paste0("some predicted probabilities do not sum to within",tol,"of 1."))
    }



    if(n > 1){
      if(t_length > 1){
        out_mat[t_ind,,] <- t(out_temp)
      } else{
        out_mat <- out_temp
      }
    } else{
      out_mat[t_ind,] <- out_temp
    }
  }

  return(out_mat)
}



#' Calculate absolute risk profiles
#'
#' This function calculates absolute risk profiles under weibull baseline hazard specifications.
#'
#' @inheritParams calc_risk_WB
#' @inheritParams simID_PW
#'
#' @return if Xmat has only one row, and t_cutoff is a scalar, then returns a 4 element row matrix
#'   of probabilities. If Xmat has \code{n} rows, then returns an \code{n} by 4 matrix of probabilities.
#'   If Xmat has \code{n} rows and t_cutoff is a vector of length \code{s}, then returns an \code{s} by 4 by \code{n} array.
#' @export
calc_risk_PW <- function(para, Xmat1, Xmat2, Xmat3, knots_list,
                         t_cutoff, tol=1e-3, frailty=TRUE,
                         type="marginal", gamma=1, model="semi-markov"){
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

  n <- max(1,nrow(Xmat1),nrow(Xmat2),nrow(Xmat3))
  t_length <- length(t_cutoff)

  if(frailty){
    nP0 <- nP01 + nP02 + nP03 + 1
    theta <- exp(para[nP0])
    if(type=="conditional" & length(gamma)==1){
      gamma <- rep(gamma,n)
    }
  } else{
    nP0 <- nP01 + nP02 +nP03
    type <- "conditional"
    gamma <- rep(1,n)
  }

  if(!is.null(Xmat1) && !(ncol(Xmat1)==0)){
    nP1 <- ncol(Xmat1)
    beta1 <- para[(1+nP0):(nP0+nP1)]
    eta1 <- as.vector(Xmat1 %*% beta1)
  } else{
    nP1 <- 0
    eta1 <- 0
  }
  if(!is.null(Xmat2) && !(ncol(Xmat2)==0)){
    nP2 <- ncol(Xmat2)
    beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    eta2 <- as.vector(Xmat2 %*% beta2)
  } else{
    nP2 <- 0
    eta2 <- 0
  }
  if(!is.null(Xmat3) && !(ncol(Xmat3)==0)){
    nP3 <- ncol(Xmat3)
    beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    eta3 <- as.vector(Xmat3 %*% beta3)
  } else{
    nP3 <- 0
    eta3 <- 0
  }

  stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss

  phi1 <- as.numeric(para[(1):(nP01)])
  phi2 <- as.numeric(para[(1+nP01):(nP01+nP02)])
  phi3 <- as.numeric(para[(1+nP01+nP02):(nP01+nP02+nP03)])

  ##NOW, USE THE ESTIMATES TO COMPUTE THE PREDICTIONS##
  ##*************************************************##

  haz <- function(t,phi,knots){
    exp(phi)[findInterval(x=t, vec=knots, left.open=TRUE)]
  }

  Haz <- function(t,phi,knots){
    rowSums(sweep(x=pw_cum_mat(t,knots),MARGIN=2,STATS=exp(phi),FUN ="*"))
  }

  ##*****************************************************************##
  ## Calculating posterior predictive density for WB baseline hazard ##
  ##*****************************************************************##


  #First, write functions that compute the integrand, which we feed into integration function

  #first, the univariate function when T1=infinity
  if(tolower(type) %in% c("marginal","m")){
    f_t2 <- function(t2,index){
      haz(t=t2,phi=phi2,knots=knots_list[[2]]) * exp(eta2[index]) *
        (1 + theta * (Haz(t=t2,phi=phi1,knots=knots_list[[1]])* exp(eta1[index]) +
                        Haz(t=t2,phi=phi2,knots=knots_list[[2]])* exp(eta2[index])))^(-theta^(-1)-1)
    }
  } else{
    f_t2 <- function(t2,index){
      gamma[index] * haz(t=t2,phi=phi2,knots=knots_list[[2]]) * exp(eta2[index]) *
        exp(-gamma[index]*(Haz(t=t2,phi=phi1,knots=knots_list[[1]])* exp(eta1[index]) +
                             Haz(t=t2,phi=phi2,knots=knots_list[[2]])* exp(eta2[index])))
    }
  }


  #next, the different regions of the joint density on the upper triangle
  if(tolower(model) == "semi-markov"){

    #if we pre-integrate t2 from t1 to t_cutoff
    if(tolower(type) %in% c("marginal","m")){
      f_joint_t1_both <- function(time_pt1,t_cutoff,index){
        haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) * (
          (1 + theta*(Haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) +
                        Haz(t=time_pt1,phi=phi2,knots=knots_list[[2]])* exp(eta2[index])))^(-theta^(-1)-1) -
            (1 + theta*(Haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) +
                          Haz(t=time_pt1,phi=phi2,knots=knots_list[[2]])* exp(eta2[index]) +
                          Haz(t=t_cutoff-time_pt1,phi=phi3,knots=knots_list[[3]]) * exp(eta3[index])))^(-theta^(-1)-1))
      }
    } else{
      f_joint_t1_both <- function(time_pt1,t_cutoff,index){
        gamma[index] * haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) *
                  exp(-gamma[index] * (Haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) +
                                         Haz(t=time_pt1,phi=phi2,knots=knots_list[[2]])* exp(eta2[index]) )) *
                  ( 1 - exp(-gamma[index] * Haz(t=t_cutoff-time_pt1,phi=phi3,knots=knots_list[[3]]) * exp(eta3[index])))
      }
    }


    #if we pre-integrate t2 from t_cutoff to infinity
    if(tolower(type) %in% c("marginal","m")){
      f_joint_t1_nonTerm <- function(time_pt1,t_cutoff,index){
        haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) *
                  (1 + theta * (Haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) +
                                  Haz(t=time_pt1,phi=phi2,knots=knots_list[[2]])* exp(eta2[index]) +
                                  Haz(t=t_cutoff-time_pt1,phi=phi3,knots=knots_list[[3]]) * exp(eta3[index])))^(-theta^(-1)-1)
      }
    } else{
      f_joint_t1_nonTerm <- function(time_pt1,t_cutoff,index){
        gamma[index] * haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) *
                  exp(-gamma[index] * (Haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) +
                                         Haz(t=time_pt1,phi=phi2,knots=knots_list[[2]])* exp(eta2[index]) +
                                         Haz(t=t_cutoff-time_pt1,phi=phi3,knots=knots_list[[3]]) * exp(eta3[index])))
      }
    }

    #if we pre-integrate t2 from t1 to infinity
    if(tolower(type) %in% c("marginal","m")){
      f_joint_t1_neither <- function(time_pt1,index){
        haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) *
                    (1 + theta*(Haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) +
                                  Haz(t=time_pt1,phi=phi2,knots=knots_list[[2]])* exp(eta2[index]) ) )^(-theta^(-1)-1)
      }
    } else{
      f_joint_t1_neither <- function(time_pt1,index){
        gamma[index] * haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) *
                  exp(-gamma[index] * (Haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) +
                                         Haz(t=time_pt1,phi=phi2,knots=knots_list[[2]])* exp(eta2[index]) ) )
      }
    }

  } else{
    stop("model must be 'semi-markov'")
  }



  if(n > 1){
    if(t_length > 1){
      out_mat <- array(dim=c(t_length,4,n),dimnames = list(paste0("t",t_cutoff),c("p_ntonly","p_both","p_tonly","p_neither"),paste0("i",1:n)))
    } else{
      out_mat <- matrix(nrow=n,ncol=4,dimnames = list(paste0("i",1:n),c("p_ntonly","p_both","p_tonly","p_neither")))
    }
  } else{
    out_mat <- matrix(nrow=t_length,ncol=4,dimnames = list(paste0("t",t_cutoff),c("p_ntonly","p_both","p_tonly","p_neither")))
  }

  for(t_ind in 1:t_length){

    p_ntonly <- sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_nonTerm, lower=0, upper=t_cutoff[t_ind],
                                                          t_cutoff=t_cutoff[t_ind], index=x)$value,
                                                error=function(cnd){return(NA)}) })

    p_tonly <- sapply(1:n,function(x){tryCatch(stats::integrate(f_t2, lower=0, upper=t_cutoff[t_ind], index=x)$value,
                                               error=function(cnd){return(NA)}) })

    p_neither_t2 <- sapply(1:n,function(x){tryCatch(stats::integrate(f_t2, lower=t_cutoff[t_ind], upper=Inf,index=x)$value,
                                                    error=function(cnd){return(NA)}) })

    p_neither_joint <- sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_neither, lower=t_cutoff[t_ind], upper=Inf, index=x)$value,
                                                       error=function(cnd){return(NA)}) })

    p_both <- sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_both, lower=0, upper=t_cutoff[t_ind],
                                                        t_cutoff=t_cutoff[t_ind], index=x)$value,
                                              error=function(cnd){return(NA)}) })

    out_temp <- cbind(p_ntonly=p_ntonly,
                      p_both=p_both,
                      p_tonly=p_tonly,
                      p_neither=p_neither_t2 + p_neither_joint)


    #I noticed that sometimes, if exactly one category has an NA, then we could back out the value
    #from the other categories. However, we don't want any to be negative.
    #I will also supply a warning if any of the rows are way off from 1.
    out_temp <- t(apply(out_temp,1,
                        function(x){
                          if(sum(is.na(x))==1){
                            x[is.na(x)]<- max(1-sum(x,na.rm=TRUE),0)
                          }
                          return(x)}))
    if(any(is.na(out_temp)) | any(abs(1-rowSums(out_temp))>tol)){
      warning(paste0("some predicted probabilities do not sum to within",tol,"of 1."))
    }


    if(n > 1){
      if(t_length > 1){
        out_mat[t_ind,,] <- t(out_temp)
      } else{
        out_mat <- out_temp
      }
    } else{
      out_mat[t_ind,] <- out_temp
    }
  }

  return(out_mat)
}



#' Get matrix of observed outcome categories
#'
#' This function returns a matrix giving the observed outcome categories of each observation at various
#'   time cutoffs.
#'
#' @inheritParams proximal_gradient_descent
#' @inheritParams calc_risk_WB
#'
#' @return a matrix or array.
#' @export
get_outcome_mat <- function(y1, y2, delta1, delta2, t_cutoff){

  n <- length(y1)
  t_length <- length(t_cutoff)

  if(n > 1){
    if(t_length > 1){
      out_mat <- array(dim=c(t_length,4,n),dimnames = list(paste0("t",t_cutoff),c("ntonly","both","tonly","neither"),paste0("i",1:n)))
    } else{
      out_mat <- matrix(nrow=n,ncol=4,dimnames = list(paste0("i",1:n),c("ntonly","both","tonly","neither")))
    }
  } else{
    out_mat <- matrix(nrow=t_length,ncol=4,dimnames = list(paste0("t",t_cutoff),c("ntonly","both","tonly","neither")))
  }

  for(t_ind in 1:t_length){

    #For cases where y=t_cutoff, I consider events that happened exactly at t_cutoff in categorization.
    neither <- t_cutoff[t_ind] < y1 | #neither
                y2 <= t_cutoff[t_ind] & delta1 == 0 & delta2 == 0 #neither
    ntonly <- y1 <= t_cutoff[t_ind] & t_cutoff[t_ind] < y2 | #ntonly
                      y2 <= t_cutoff[t_ind] & delta1 == 1 & delta2 == 0 #ntonly
    tonly <- y2 <= t_cutoff[t_ind] & delta1 == 0 & delta2 == 1 #tonly
    both <- y2 <= t_cutoff[t_ind] & delta1 == 1 & delta2 == 1 #both



    out_temp <- cbind(ntonly=as.numeric(ntonly),
                      both=as.numeric(both),
                      tonly=as.numeric(tonly),
                      neither=as.numeric(neither))

    if(n > 1){
      if(t_length > 1){
        out_mat[t_ind,,] <- t(out_temp)
      } else{
        out_mat <- out_temp
      }
    } else{
      out_mat[t_ind,] <- out_temp
    }
  }

  return(out_mat)
}


#' Get inverse probability of censoring weights
#'
#' This function returns a vector of inverse probability of censoring weights from an unadjusted Cox model
#'   for censoring times.
#'
#' @inheritParams proximal_gradient_descent
#' @inheritParams calc_risk_WB
#'
#' @return a vector.
#' @export
get_ipcw_mat <- function(y2,delta2,t_cutoff){

  # browser()
  n <- length(y2)
  t_length <- length(t_cutoff)

  #this is Ghat, a non-parametric model of the 'survival distribution' of censoring var C.
  sfcens <- survival::survfit(survival::Surv(y2, delta2==0) ~ 1)

  #* want to get Ghat(z) where
  #* z=y2- if y2<=s and delta2=1,
  #* z=Inf if y2<s and delta2=0, (aka, 1/Ghat(z)=0, I know they're subtly different)
  #* z=s if s=y2 and delta2=0,
  #* z=s if s<y2
  #*
  #* so, we define
  #* z = min(y2,s)
  #* then change z=y2- if y2<=s and delta2=1
  #* then change z=Inf if y2< s and delta2=0 (again, technically change 1/Ghat(z)=0 for these obs but still)
  #* and that should do it! (I hope)

  ipcw_mat <- matrix(nrow=n,ncol=t_length,dimnames = list(paste0("i",1:n),paste0("t",t_cutoff)))
  for(t_ind in 1:t_length){
    #vector of min(ttilde,t)
    y_last <- pmin(y2,t_cutoff[t_ind])
    if(sum(y2 <= t_cutoff[t_ind] & delta2==1)>0){
      y_last[y2 <= t_cutoff[t_ind] & delta2==1] <- y_last[y2 <= t_cutoff[t_ind] & delta2==1] - 1e-8
    }
    y_last_cens <- rep(NA,n)
    y_last_cens[order(y_last)] <- summary(sfcens, times = y_last, extend=TRUE)$surv
    if(sum(y_last < t_cutoff[t_ind] & delta2==0) > 0){
      y_last_cens[y2 < t_cutoff[t_ind] & delta2==0] <- Inf
    }
    ipcw_mat[,t_ind] <- 1/y_last_cens
  }

  ipcw_mat

}


#' Compute prediction performance score
#'
#' This function takes in all of the ingredients needed for prediction validation,
#'   and returns the corresponding scores.
#'
#' @param outcome_mat Output from get_outcome_mat function
#' @param pred_mat Output from calc_risks function
#' @param ipcw_mat Output from get_ipcw_mat function
#' @param score String indicating whether 'brier' score, or 'entropy' should be computed.
#'
#' @return a vector.
#' @export
compute_score <- function(outcome_mat, pred_mat, ipcw_mat, score="brier"){
  #this function is for brier and kl scores, while the below function is for the AUC
  # browser()
  if(length(dim(outcome_mat))==3){
    if(tolower(score) %in% c("brier")){
      out <- apply( t(apply((outcome_mat - pred_mat)^2,
                            MARGIN = c(1,3),
                            FUN = sum)) *
                      ipcw_mat, MARGIN = 2, FUN = mean)
    } else{
      out <- apply( t(apply(-outcome_mat*log(pred_mat),
                            MARGIN = c(1,3),
                            FUN = sum)) *
                      ipcw_mat, MARGIN = 2, FUN = mean)
    }
  } else{ #this must mean that there is only a single time cutoff, so the input mats are n by 4 and weights is just a vector
    if(tolower(score) %in% c("brier")){
      out <- colMeans( apply((outcome_mat - pred_mat)^2,
                            MARGIN = c(1),
                            FUN = sum) * ipcw_mat)
    } else{
      out <- colMeans( apply(-outcome_mat*log(pred_mat),
                         MARGIN = c(1),
                         FUN = sum) * ipcw_mat)
    }
  }
  return(out)

}


compute_auc <- function(simData,t_cutoff, pred_mat){
  #For now, this one is only implemented for when these are just matrices, aka for a single choice of t.
  #need pred_mat to be from same value as t_cutoff, however!

  # browser()

  outcomes <- simData %>% select(y1,delta1,y2,delta2) %>%
    mutate(nonterm_comp_risk_time = ifelse(y1 < y2, y1, y2),
           comp_risk_event = ifelse( (y1 == y2) & delta2==1,2,ifelse(y1 == y2 & delta2==0,0,1))
    )

  #treats terminal outcome as a competing risk
  ROC_nonterm <- timeROC::timeROC(T=outcomes$nonterm_comp_risk_time,
                           delta=outcomes$comp_risk_event,
                           marker=pred_mat[,"p_ntonly"] + pred_mat[,"p_both"],
                           cause=1,weighting="marginal",
                           times=t_cutoff,
                           iid=TRUE)

  ROC_term <- timeROC::timeROC(T=outcomes$y2,
                        delta=outcomes$delta2,
                        marker=pred_mat[,"p_tonly"] + pred_mat[,"p_both"],
                        cause=1,weighting="marginal",
                        times=t_cutoff,
                        iid=TRUE)

  return(c(AUC_nonterm = ROC_nonterm$AUC_1[2],
           AUC_term = ROC_term$AUC[2]))
}
