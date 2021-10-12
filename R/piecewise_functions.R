#' Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under piecewise constant baseline hazard, gamma subject-specific frailty,
#'   and semi-Markov transition assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#'
#' @inheritParams proximal_gradient_descent
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nlogLikPW_ID_frail_SM <- function(para, y1, y2, delta1, delta2,
                                  Xmat1, Xmat2, Xmat3,
                                  basis1, basis2, basis3){
  # browser()

  #count of the number of intervals in each list
  num_int <- c(ncol(basis1),ncol(basis2),ncol(basis3))
  n <- length(y1)

  phi1 <- para[1:num_int[1]]
  phi2 <- para[(1+num_int[1]):(num_int[1]+num_int[2])]
  phi3 <- para[(1+num_int[1]+num_int[2]):(num_int[1]+num_int[2]+num_int[3])]
  theta <- exp(para[(1+num_int[1]+num_int[2]+num_int[3]):(num_int[1]+num_int[2]+num_int[3]+1)])

  if(!is.null(Xmat1) && !(ncol(Xmat1)==0)){
    nP1 <- ncol(Xmat1)
    beta1 <- para[(num_int[1]+num_int[2]+num_int[3]+1+1):(num_int[1]+num_int[2]+num_int[3]+1+nP1)]
    eta1 <- as.vector(Xmat1 %*% beta1)
  } else{
    nP1 <- 0
    eta1 <- 0
  }
  if(!is.null(Xmat2) && !(ncol(Xmat2)==0)){
    nP2 <- ncol(Xmat2)
    beta2 <- para[(1+num_int[1]+num_int[2]+num_int[3]+nP1+1):(num_int[1]+num_int[2]+num_int[3]+1+nP1+nP2)]
    eta2 <- as.vector(Xmat2 %*% beta2)
  } else{
    nP2 <- 0
    eta2 <- 0
  }
  if(!is.null(Xmat3) && !(ncol(Xmat3)==0)){
    nP3 <- ncol(Xmat3)
    beta3 <- para[(1+num_int[1]+num_int[2]+num_int[3]+nP1+nP2+1):(num_int[1]+num_int[2]+num_int[3]+1+nP1+nP2+nP3)]
    eta3 <- as.vector(Xmat3 %*% beta3)
  } else{
    nP3 <- 0
    eta3 <- 0
  }

  #matrix with each row being an observation, and each column being an interval, capturing how much time is accrued in each interval
  # knots1_mat <- matrix(c(knots1_diff,0),nrow=n,ncol=num_int[1],byrow = TRUE)
  # knots2_mat <- matrix(c(knots2_diff,0),nrow=n,ncol=num_int[2],byrow = TRUE)

  #an n length vector of integers indicating which interval the observation is in
  cut_cats1 <- rowSums(basis1!=0)
  cut_cats2 <- rowSums(basis2!=0)
  cut_cats3 <- rowSums(basis3!=0) #ASSUMES SEMIMARKOV FOR THE MOMENT
  cut_cats3[cut_cats3==0] <- 1 #for bookkeeping, put people who never accrue time into interval 1

  #computes n length vector of the hazard at the failure time of each observation
  haz01 <- exp(phi1)[cut_cats1]
  haz02 <- exp(phi2)[cut_cats2]
  haz03 <- exp(phi3)[cut_cats3]

  #fill in the spots that accrue no time in the third hazard with NAs
  # haz03 <- rep(1,n) #placeholder, never gets used
  # haz03[delta1==1] <- exp(phi3)[cut_cats3]

  #computes n length vector of the cumulative hazard at the failure time of each observation
  Lambda01 <- basis1 %*% exp(phi1)
  Lambda02 <- basis2 %*% exp(phi2)
  Lambda03 <- basis3 %*% exp(phi3)

  # browser()

  ll <- sum(
    delta1 * (log(haz01) + eta1) +
      (1-delta1)*delta2*(log(haz02) + eta2) +
      delta1*delta2*(log(haz03) + eta3 + log1p(theta)) -
      (theta^(-1) + delta1 + delta2)* as.vector(log1p(theta*(Lambda01*exp(eta1) +
                                                                 Lambda02*exp(eta2) +
                                                                 Lambda03*exp(eta3)
      )
      )
      )
  )
  return(-ll)
}


#' Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under piecewise constant baseline hazard, gamma subject-specific frailty,
#'   and Markov transition assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#'
#' @inheritParams proximal_gradient_descent
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nlogLikPW_ID_frail_M <- function(para, y1, y2, delta1, delta2,
                                 Xmat1 = NULL, Xmat2 = NULL, Xmat3 = NULL,
                                 basis1, basis2, basis3, basis3_y1){
  # browser()

  #count of the number of intervals in each list
  num_int <- c(ncol(basis1),ncol(basis2),ncol(basis3))
  n <- length(y1)

  phi1 <- para[1:num_int[1]]
  phi2 <- para[(1+num_int[1]):(num_int[1]+num_int[2])]
  phi3 <- para[(1+num_int[1]+num_int[2]):(num_int[1]+num_int[2]+num_int[3])]
  theta <- exp(para[(1+num_int[1]+num_int[2]+num_int[3]):(num_int[1]+num_int[2]+num_int[3]+1)])

  if(!is.null(Xmat1) && !(ncol(Xmat1)==0)){
    nP1 <- ncol(Xmat1)
    beta1 <- para[(num_int[1]+num_int[2]+num_int[3]+1+1):(num_int[1]+num_int[2]+num_int[3]+1+nP1)]
    eta1 <- as.vector(Xmat1 %*% beta1)
  } else{
    nP1 <- 0
    eta1 <- 0
  }
  if(!is.null(Xmat2) && !(ncol(Xmat2)==0)){
    nP2 <- ncol(Xmat2)
    beta2 <- para[(1+num_int[1]+num_int[2]+num_int[3]+nP1+1):(num_int[1]+num_int[2]+num_int[3]+1+nP1+nP2)]
    eta2 <- as.vector(Xmat2 %*% beta2)
  } else{
    nP2 <- 0
    eta2 <- 0
  }
  if(!is.null(Xmat3) && !(ncol(Xmat3)==0)){
    nP3 <- ncol(Xmat3)
    beta3 <- para[(1+num_int[1]+num_int[2]+num_int[3]+nP1+nP2+1):(num_int[1]+num_int[2]+num_int[3]+1+nP1+nP2+nP3)]
    eta3 <- as.vector(Xmat3 %*% beta3)
  } else{
    nP3 <- 0
    eta3 <- 0
  }

  #matrix with each row being an observation, and each column being an interval, capturing how much time is accrued in each interval
  # knots1_mat <- matrix(c(knots1_diff,0),nrow=n,ncol=num_int[1],byrow = TRUE)
  # knots2_mat <- matrix(c(knots2_diff,0),nrow=n,ncol=num_int[2],byrow = TRUE)

  #an n length vector of integers indicating which interval the observation is in
  cut_cats1 <- rowSums(basis1!=0)
  cut_cats2 <- rowSums(basis2!=0)
  cut_cats3 <- rowSums(basis3!=0) #ASSUMES MARKOV FOR THE MOMENT

  #computes n length vector of the hazard at the failure time of each observation
  haz01 <- exp(phi1)[cut_cats1]
  haz02 <- exp(phi2)[cut_cats2]
  haz03 <- exp(phi3)[cut_cats3]

  #computes n length vector of the cumulative hazard at the failure time of each observation
  Lambda01 <- basis1 %*% exp(phi1)
  Lambda02 <- basis2 %*% exp(phi2)
  Lambda03 <- (basis3 - basis3_y1) %*% exp(phi3)

  #browser()
  ll <- sum(
    delta1 * (log(haz01) + eta1) +
      (1-delta1)*delta2*(log(haz02) + eta2) +
      delta1*delta2*(log(haz03) + eta3 + log1p(theta)) -
      (theta^(-1) + delta1 + delta2)* as.vector(log1p(theta*(Lambda01*exp(eta1) +
                                                                 Lambda02*exp(eta2) +
                                                                 Lambda03*exp(eta3)
      )
      )
      )
  )
  return(-ll)
}


#' Gradient of Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the gradient of the negative log-likelihood for the illness-death model,
#'   under piecewise constant baseline hazard, gamma subject-specific frailty,
#'   and semi-Markov transition assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#'
#' @inheritParams proximal_gradient_descent
#'
#' @return Returns numeric vector of same length as \code{para} with sum of gradient contributions
#'   for the negative log likelihood.
#' @export
ngradPW_ID_frail_SM <- function(para, y1, y2, delta1, delta2,
                                Xmat1 = NULL, Xmat2 = NULL, Xmat3 = NULL,
                                basis1, basis2, basis3){
  # browser()

  #count of the number of intervals in each list
  num_int <- c(ncol(basis1),ncol(basis2),ncol(basis3))
  n <- length(y1)

  phi1 <- para[1:num_int[1]]
  phi2 <- para[(1+num_int[1]):(num_int[1]+num_int[2])]
  phi3 <- para[(1+num_int[1]+num_int[2]):(num_int[1]+num_int[2]+num_int[3])]
  theta <- exp(para[(1+num_int[1]+num_int[2]+num_int[3]):(num_int[1]+num_int[2]+num_int[3]+1)])
  h <- log(theta)

  if(!is.null(Xmat1) && !(ncol(Xmat1)==0)){
    nP1 <- ncol(Xmat1)
    beta1 <- para[(num_int[1]+num_int[2]+num_int[3]+1+1):(num_int[1]+num_int[2]+num_int[3]+1+nP1)]
    eta1 <- as.vector(Xmat1 %*% beta1)
  } else{
    nP1 <- 0
    eta1 <- 0
  }
  if(!is.null(Xmat2) && !(ncol(Xmat2)==0)){
    nP2 <- ncol(Xmat2)
    beta2 <- para[(1+num_int[1]+num_int[2]+num_int[3]+nP1+1):(num_int[1]+num_int[2]+num_int[3]+1+nP1+nP2)]
    eta2 <- as.vector(Xmat2 %*% beta2)
  } else{
    nP2 <- 0
    eta2 <- 0
  }
  if(!is.null(Xmat3) && !(ncol(Xmat3)==0)){
    nP3 <- ncol(Xmat3)
    beta3 <- para[(1+num_int[1]+num_int[2]+num_int[3]+nP1+nP2+1):(num_int[1]+num_int[2]+num_int[3]+1+nP1+nP2+nP3)]
    eta3 <- as.vector(Xmat3 %*% beta3)
  } else{
    nP3 <- 0
    eta3 <- 0
  }

  #matrix with each row being an observation, and each column being an interval, capturing how much time is accrued in each interval
  # knots1_mat <- matrix(c(knots1_diff,0),nrow=n,ncol=num_int[1],byrow = TRUE)
  # knots2_mat <- matrix(c(knots2_diff,0),nrow=n,ncol=num_int[2],byrow = TRUE)

  #an n length vector of integers indicating which interval the observation is in
  cut_cats1 <- rowSums(basis1!=0)
  cut_cats2 <- rowSums(basis2!=0)
  cut_cats3 <- rowSums(basis3!=0) #ASSUMES SEMIMARKOV FOR THE MOMENT
  cut_cats3[cut_cats3==0] <- 1

  #computes n length vector of the cumulative hazard at the failure time of each observation
  Lambda01 <- basis1 %*% exp(phi1)
  Lambda02 <- basis2 %*% exp(phi2)
  Lambda03 <- basis3 %*% exp(phi3)

  #code to multiply each row of basis matrix elementwise by the vector exp(phig)
  #this gives us exactly what we want for Lambda0gPhi (faster version of 'sweep')
  #https://stackoverflow.com/questions/49462591/r-multiply-every-row-of-df-or-matrix-with-a-vector/49462658

  Lambda01Phi <- basis1*t(exp(phi1))[rep_len(1,nrow(Lambda01)),]
  Lambda02Phi <- basis2*t(exp(phi2))[rep_len(1,nrow(Lambda02)),]
  Lambda03Phi <- basis3*t(exp(phi3))[rep_len(1,nrow(Lambda03)),]

  #browser()

  AVec <- Lambda01 * exp(eta1) + Lambda02 * exp(eta2) + Lambda03 * exp(eta3)
  commonVec <- (exp(-h) + delta1 + delta2) / (1+exp(h)*AVec)

  #first term is just counting how many people experience the nonterminal event in each interval
  temp_val <- tapply(delta1,cut_cats1,sum)
  if(length(temp_val) != num_int[1]){stop("It looks like you do not observe events in every interval of the first transition")}
  score_phi1 <- temp_val - as.vector( t(Lambda01Phi) %*% (commonVec * exp(h + eta1)) )

  #first term is just counting how many people experience the terminal event without the nonterminal in each interval
  temp_val <- tapply((1-delta1)*delta2,cut_cats2,sum)
  if(length(temp_val) != num_int[2]){stop("It looks like you do not observe events in every interval of the second transition")}
  score_phi2 <- temp_val - as.vector( t(Lambda02Phi) %*% (commonVec * exp(h + eta2)) )

  #first term is just counting how many people experience the terminal event after the nonterminal in each interval
  temp_val <- tapply((delta1)*delta2,cut_cats3,sum)
  if(length(temp_val) != num_int[3]){stop("It looks like you do not observe events in every interval of the third transition")}
  score_phi3 <- temp_val - as.vector(  t(Lambda03Phi) %*% (commonVec * exp(h + eta3)) )

  #h (what ina calls u1)
  score_h <- sum(exp(h)*(delta1 * delta2/(1+exp(h)) + log1p(exp(h) * AVec)/exp(2*h) - commonVec * AVec))

  #beta1 (what ina calls u2)
  if(nP1 == 0){ score_beta1 <- NULL } else{
    score_beta1 <-  t(Xmat1) %*% (delta1 - commonVec * Lambda01 * exp(h + eta1))
  }

  #beta2 (what ina calls u3)
  if(nP2 == 0){ score_beta2 <- NULL } else{
    score_beta2 <- t(Xmat2) %*% ( (1-delta1) * delta2 - commonVec * Lambda02 * exp(h + eta2))
  }

  #beta3 (what ina calls u4)
  if(nP3 == 0){ score_beta3 <- NULL } else{
    score_beta3 <-  t(Xmat3) %*% (delta1 * delta2 - commonVec * Lambda03 * exp(h + eta3))
  }

  return(-c(score_phi1,score_phi2,score_phi3,score_h,score_beta1,score_beta2,score_beta3))

}



#' Gradient of Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the gradient of the negative log-likelihood for the illness-death model,
#'   under piecewise constant baseline hazard, gamma subject-specific frailty,
#'   and Markov transition assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#'
#' @inheritParams proximal_gradient_descent
#'
#' @return Returns numeric vector of same length as \code{para} with sum of gradient contributions
#'   for the negative log likelihood.
#' @export
ngradPW_ID_frail_M <- function(para, y1, y2, delta1, delta2,
                               Xmat1 = NULL, Xmat2 = NULL, Xmat3 = NULL,
                               basis1, basis2, basis3, basis3_y1){
  # browser()

  #count of the number of intervals in each list
  num_int <- c(ncol(basis1),ncol(basis2),ncol(basis3))
  n <- length(y1)

  phi1 <- para[1:num_int[1]]
  phi2 <- para[(1+num_int[1]):(num_int[1]+num_int[2])]
  phi3 <- para[(1+num_int[1]+num_int[2]):(num_int[1]+num_int[2]+num_int[3])]
  theta <- exp(para[(1+num_int[1]+num_int[2]+num_int[3]):(num_int[1]+num_int[2]+num_int[3]+1)])
  h <- log(theta)

  if(!is.null(Xmat1) && !(ncol(Xmat1)==0)){
    nP1 <- ncol(Xmat1)
    beta1 <- para[(num_int[1]+num_int[2]+num_int[3]+1+1):(num_int[1]+num_int[2]+num_int[3]+1+nP1)]
    eta1 <- as.vector(Xmat1 %*% beta1)
  } else{
    nP1 <- 0
    eta1 <- 0
  }
  if(!is.null(Xmat2) && !(ncol(Xmat2)==0)){
    nP2 <- ncol(Xmat2)
    beta2 <- para[(1+num_int[1]+num_int[2]+num_int[3]+nP1+1):(num_int[1]+num_int[2]+num_int[3]+1+nP1+nP2)]
    eta2 <- as.vector(Xmat2 %*% beta2)
  } else{
    nP2 <- 0
    eta2 <- 0
  }
  if(!is.null(Xmat3) && !(ncol(Xmat3)==0)){
    nP3 <- ncol(Xmat3)
    beta3 <- para[(1+num_int[1]+num_int[2]+num_int[3]+nP1+nP2+1):(num_int[1]+num_int[2]+num_int[3]+1+nP1+nP2+nP3)]
    eta3 <- as.vector(Xmat3 %*% beta3)
  } else{
    nP3 <- 0
    eta3 <- 0
  }

  #matrix with each row being an observation, and each column being an interval, capturing how much time is accrued in each interval
  # knots1_mat <- matrix(c(knots1_diff,0),nrow=n,ncol=num_int[1],byrow = TRUE)
  # knots2_mat <- matrix(c(knots2_diff,0),nrow=n,ncol=num_int[2],byrow = TRUE)

  #an n length vector of integers indicating which interval the observation is in
  cut_cats1 <- rowSums(basis1!=0)
  cut_cats2 <- rowSums(basis2!=0)
  cut_cats3 <- rowSums(basis3!=0) #ASSUMES MARKOV FOR THE MOMENT

  #computes n length vector of the cumulative hazard at the failure time of each observation
  Lambda01 <- basis1 %*% exp(phi1)
  Lambda02 <- basis2 %*% exp(phi2)
  Lambda03 <- (basis3 - basis3_y1) %*% exp(phi3)

  #code to multiply each row of Lambda0g elementwise by the vector exp(phig)
  #this gives us exactly what we want for Lambda0gPhi (faster version of 'sweep')
  #https://stackoverflow.com/questions/49462591/r-multiply-every-row-of-df-or-matrix-with-a-vector/49462658

  Lambda01Phi <- basis1*t(exp(phi1))[rep_len(1,nrow(Lambda01)),]
  Lambda02Phi <- basis2*t(exp(phi2))[rep_len(1,nrow(Lambda02)),]
  Lambda03Phi <- (basis3 - basis3_y1)*t(exp(phi3))[rep_len(1,nrow(Lambda03)),]

  #browser()

  AVec <- Lambda01 * exp(eta1) + Lambda02 * exp(eta2) + Lambda03 * exp(eta3)
  commonVec <- (exp(-h) + delta1 + delta2) / (1+exp(h)*AVec)

  #first term is just counting how many people experience the nonterminal event in each interval
  temp_val <- tapply(delta1,cut_cats1,sum)
  if(length(temp_val) != num_int[1]){stop("It looks like you do not observe events in every interval of the first transition")}
  score_phi1 <- temp_val - as.vector( t(Lambda01Phi) %*% (commonVec * exp(h + eta1)) )

  #first term is just counting how many people experience the terminal event without the nonterminal in each interval
  temp_val <- tapply((1-delta1)*delta2,cut_cats2,sum)
  if(length(temp_val) != num_int[2]){stop("It looks like you do not observe events in every interval of the second transition")}
  score_phi2 <- temp_val - as.vector( t(Lambda02Phi) %*% (commonVec * exp(h + eta2)) )

  #first term is just counting how many people experience the terminal event after the nonterminal in each interval
  temp_val <- tapply((delta1)*delta2,cut_cats3,sum)
  if(length(temp_val) != num_int[3]){stop("It looks like you do not observe events in every interval of the third transition")}
  score_phi3 <- temp_val - as.vector(  t(Lambda03Phi) %*% (commonVec * exp(h + eta3)) )

  #h (what ina calls u1)
  score_h <- sum(exp(h)*(delta1 * delta2/(1+exp(h)) + log1p(exp(h) * AVec)/exp(2*h) - commonVec * AVec))

  #beta1 (what ina calls u2)
  if(nP1 == 0){ score_beta1 <- NULL } else{
    score_beta1 <-  t(Xmat1) %*% (delta1 - commonVec * Lambda01 * exp(h + eta1))
  }

  #beta2 (what ina calls u3)
  if(nP2 == 0){ score_beta2 <- NULL } else{
    score_beta2 <- t(Xmat2) %*% ( (1-delta1) * delta2 - commonVec * Lambda02 * exp(h + eta2))
  }

  #beta3 (what ina calls u4)
  if(nP3 == 0){ score_beta3 <- NULL } else{
    score_beta3 <-  t(Xmat3) %*% (delta1 * delta2 - commonVec * Lambda03 * exp(h + eta3))
  }

  return(-c(score_phi1,score_phi2,score_phi3,score_h,score_beta1,score_beta2,score_beta3))

}

# FUNCTION TO GET 'BASIS' MATRIX OF TIME ACCRUED IN EACH ARM ------------

#' Generate Piecewise Constant Cumulative Hazard Matrix
#'
#' This helper function takes in a vector of event times, and
#'   generates a matrix dividing each time into a vector giving the time accrued within each
#'   interval between consecutive elements of a vector of knots.
#'
#' @param y vector of event times.
#' @param knots increasing vector of cutpoints
#'
#' @return a numeric matrix, with rows corresponding to elements of y.
#' @export
pw_cum_mat <- function(y, knots){
  # browser()

  if(knots[1] != 0){knots <- c(0,knots)}
  knots_diff <- diff(knots)

  #count of the number of intervals in each list
  num_int <- c(length(knots))
  n <- length(y)

  #matrix with each row being an observation, and each column being an interval, capturing how much time is accrued in each interval
  knots_mat <- matrix(c(knots_diff,0),nrow=n,ncol=num_int,byrow = TRUE)

  #an n length vector of integers indicating which interval the observation is in
  cut_cats <- findInterval(x = y, vec = knots)

  #an n length vector capturing the residual time accrued in the final interval of each observation
  y_last <- y-knots[cut_cats]

  #a loop through each observation to finalize the matrix of time intervals
  for(i in 1:n){
    knots_mat[i,cut_cats[i]] <- y_last[i]
    knots_mat[i,-c(1:cut_cats[i])] <- 0
  }
  return(knots_mat)
}

