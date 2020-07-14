#' The function that simulates independent/cluster-correlated semi-competing
#'   risks data under semi-Markov PEM/PEM-MVN models.
#'
#' @param id A vector of cluster information for \code{n} subjects.
#'   The cluster membership must be set to consecutive positive integers, \eqn{1:J}.
#'   Required only when generating clustered data.
#' @param x1,x2,x3 Covariate matrices with \code{n} rows.
#' @param beta1.true,beta2.true,beta3.true Vectors of true regression parameter values.
#'   The length of each vector should equal the number of columns in the corresponding covariate matrix.
#' @param phi1.true,phi2.true,phi3.true Vectors of true baseline parameter values.
#' @param theta.true True value for \eqn{\theta}.
#' @param SigmaV.true True value for covariance matrix of MVN cluster-level random effects.
#'   Required only when generating clustered data. Should be a numeric \eqn{J\times J} matrix.
#' @param cens A numeric vector of two elements. The right censoring times are generated from Uniform(\eqn{cens[1]}, \eqn{cens[2]}).
#' @param knots_list A list containing three numeric vectors, representing the breakpoints
#'   of the piecewise specification for each transition hazard (excluding 0).
#'
#' @return returns a data.frame containing semi-competing risks outcomes from \code{n} subjects.
#'   It is of dimension \eqn{n\times 4}: the columns correspond to \eqn{y_1}, \eqn{\delta_1}, \eqn{y_2}, \eqn{\delta_2}. \cr
#'   \itemize{
#'   \item{y1}{a vector of \code{n} times to the non-terminal event}
#'   \item{y2}{a vector of \code{n} times to the terminal event}
#'   \item{delta1}{a vector of \code{n} censoring indicators for the non-terminal event time (1=event occurred, 0=censored)}
#'   \item{delta2}{a vector of \code{n} censoring indicators for the terminal event time (1=event occurred, 0=censored)}
#'   }
#'
#' @export
simID_PW <- function (id = NULL, x1, x2, x3, beta1.true, beta2.true, beta3.true,
                      phi1.true, phi2.true, phi3.true, theta.true, SigmaV.true = NULL, cens, knots_list)
{
  if (!is.null(id) & is.null(SigmaV.true)) {
    stop("SigmaV.true must be given to simulate correlated data")
  }
  else {
    n <- dim(x1)[1]
    p1 <- dim(x1)[2]
    p2 <- dim(x2)[2]
    p3 <- dim(x3)[2]
    if (theta.true > 0) {
      gamma.true <- stats::rgamma(n, 1/theta.true, 1/theta.true)
    }
    if (theta.true == 0) {
      gamma.true <- rep(1, n)
    }
    if (is.null(id)) {
      #      LP1 <- LP2 <- LP3 <- rep(0,n)
      LP1 <- as.vector(beta1.true %*% t(x1))
      LP2 <- as.vector(beta2.true %*% t(x2))
      LP3 <- as.vector(beta3.true %*% t(x3))
    }
    if (!is.null(id)) {
      J <- length(unique(id))
      nj <- as.vector(table(id))
      Vmat <- MASS::mvrnorm(J, rep(0, 3), SigmaV.true)
      LP1 <- as.vector(beta1.true %*% t(x1) + rep(Vmat[,
                                                       1], nj))
      LP2 <- as.vector(beta2.true %*% t(x2) + rep(Vmat[,
                                                       2], nj))
      LP3 <- as.vector(beta3.true %*% t(x3) + rep(Vmat[,
                                                       3], nj))
    }

    knots_list <- as.list(knots_list) #make it so that index can be made rowwise whether there is 1 col or 3
    if(length(knots_list) == 1){

      #consider testing and adding this change eventuallyyyy
      # if(knots_list[[1]][1] != 0){knots_list[[1]] <- c(0,knots_list[[1]])}
      # knots03 <- knots02 <- knots01 <- knots_list[[1]]

      knots03 <- knots02 <- knots01 <- c(0,knots_list[[1]])
      knots3_diff <- knots2_diff <- knots1_diff <- diff(knots01)
    } else if(length(knots_list) == 3){
      knots01 <- c(0,knots_list[[1]])
      knots1_diff <- diff(knots01)
      knots02 <- c(0,knots_list[[2]])
      knots2_diff <- diff(knots02)
      knots03 <- c(0,knots_list[[3]])
      knots3_diff <- diff(knots03)
    } else{
      stop("knots must be either vector of knots lambdas,
           or matrix with three columns corresponding to hazard-specific knots")
    }


    Rind <- NULL
    R <- sapply(X = 1:n, FUN = function(x){rpwexp(n=1,
                                                  rate = exp(phi1.true + LP1[x] + log(gamma.true[x])),
                                                  intervals = knots1_diff,cumulative = FALSE)})
    D <- sapply(X = 1:n, FUN = function(x){rpwexp(n=1,
                                                  rate = exp(phi2.true + LP2[x] + log(gamma.true[x])),
                                                  intervals = knots2_diff,cumulative = FALSE)})
    yesR <- R < D
    D[yesR] <- R[yesR] + sapply(X = which(yesR), FUN = function(x){rpwexp(n=1,
                                                                          rate = exp(phi3.true + LP3[x] + log(gamma.true[x])),
                                                                          intervals = knots3_diff,cumulative = FALSE)})
    delta1 <- rep(NA, n)
    delta2 <- rep(NA, n)
    y1 <- R
    y2 <- D
    Cen <- stats::runif(n, cens[1], cens[2])
    ind01 <- which(D < R & D < Cen)
    y1[ind01] <- D[ind01]
    delta1[ind01] <- 0
    delta2[ind01] <- 1
    ind10 <- which(R < D & R < Cen & D >= Cen)
    y2[ind10] <- Cen[ind10]
    delta1[ind10] <- 1
    delta2[ind10] <- 0
    ind00 <- which(R >= Cen & D >= Cen)
    y1[ind00] <- Cen[ind00]
    y2[ind00] <- Cen[ind00]
    delta1[ind00] <- 0
    delta2[ind00] <- 0
    ind11 <- which(R < Cen & D < Cen & R < D)
    delta1[ind11] <- 1
    delta2[ind11] <- 1
    ret <- data.frame(cbind(y1, delta1, y2, delta2))
    return(ret)
  }
}





#' Simulate univariate time-to-event data under piecewise exponential hazard
#'
#' \code{rpwexp} generates \code{n} independent, identically distributed event times under
#'   a piecewise exponential model.
#'
#' @param n Number of observations.
#' If length(n) > 1, the length is taken to be the number required.
#' @param rate Vector containing exponential failure rates in intervals described by
#' \code{intervals}
#' @param intervals Vector containing positive values indicating interval lengths where
#' the exponential rates provided in \code{rate} apply. Note that the length of
#' \code{intervals} is 1 less than that of \code{rate} and that the final value rate
#' in \code{rate} applies after time `sum(intervals)`.
#' @param cumulative \code{FALSE} (the default) generates \code{n} independent,
#' identically distributed piecewise exponential failure rates according to the distribution
#' specified by \code{intervals} and \code{rate}. \code{TRUE} generates independent
#' inter-arrival times with the rates of arrival in each interval specified by
#' \code{intervals} determined by \code{rate}.
#'
#' @export
rpwexp <- function(n, rate=1, intervals=NULL, cumulative=FALSE){
  if(is.null(intervals)){
    if (cumulative){return(cumsum(stats::rexp(n,rate[1])))}else
      return(stats::rexp(n,rate[1]))}
  k <- length(rate)
  if (k==1){
    if(cumulative){return(cumsum(stats::rexp(n,rate)))}else
      return(stats::rexp(n,rate))
  }
  if (length(intervals) < k-1) stop("length(intervals) must be at least length(rate) - 1")
  tx <- 0
  j <- 1
  times <- array(0,n)
  timex <- cumsum(intervals)
  indx <- array(TRUE,n)
  for(i in 1:k){
    nindx <- sum(indx)
    if (nindx==0) break
    increment <- stats::rexp(nindx,rate[i])
    if (cumulative) times[indx] <- tx + cumsum(increment)
    else times[indx] <- tx + increment
    if (i<k){
      tx <- timex[i]
      indx <- (times > timex[i])
    }
  }
  return(times)
}


