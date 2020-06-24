#' Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#'
#' @param para A numeric vector of parameters, arranged as follows:
#'   the first \eqn{k_1+k_2+k_3} elements correspond to the baseline hazard parameters,
#'   then the \eqn{k_1+k_2+k_3+1} element corresponds to the gamma frailty log-variance parameter,
#'   then the last\eqn{q_1+q_2+q_3} elements correspond with the regression parameters.
#' @param y1,y2 Numeric vectors of length \eqn{n} with (possibly censored) non-terminal and terminal event times
#' @param delta1,delta2 Numeric vectors of length \eqn{n}  with indicators of 1 if the event was observed and 0 otherwise
#' @param Xmat1,Xmat2,Xmat3 Numeric matrices with \eqn{n} rows and \eqn{q_1,q_2,q_3} columns containing covariates.
#' @param hazard String specifying the form of the baseline hazard.
#' @param frailty Boolean indicating whether a gamma distributed subject-specific frailty should
#'   be included. Currently this must be set to TRUE.
#' @param model String specifying the transition assumption
#' @param basis1,basis2,basis3,basis3_y1 Numeric matrices with \eqn{n} rows and \eqn{k_1,k_2,k_3} columns
#'   with piecewise/spline basis function values at the corresponding \code{y1} and \code{y2} values.
#'   Under semi-Markov model, basis3 represents basis derived from \eqn{y_2-y_1} and \code{basis3_y1} is unused,
#'   while under Markov model, basis3 represents basis derived from \eqn{y_2} and \code{basis3_y1} is from \eqn{y_1}
#'   Not used under Weibull model.
#' @param dbasis1,dbasis2,dbasis3 Numeric matrices with \eqn{n} rows and \eqn{k_1,k_2,k_3} columns
#'   with piecewise/spline basis function derivative values at the corresponding \code{y1} and \code{y2} values.
#'   Used only under Royston-Parmar model.
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nll_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                     hazard, frailty, model,
                     basis1, basis2, basis3, basis3_y1,
                     dbasis1, dbasis2, dbasis3){

  # browser()

  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  nP01 <- if(!is.null(basis1)) ncol(basis1) else 0
  nP02 <- if(!is.null(basis2)) ncol(basis2) else 0
  nP03 <- if(!is.null(basis3)) ncol(basis3) else 0
  n <- length(y1)

  if(tolower(hazard) == "weibull"){
    if(frailty){
      nP0 <- 7
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
      nll <- nlogLikWB_ID_frail(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                                X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                                X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0),
                                model = tolower(model))
    } else{stop("non-frailty weibull model not yet implemented")}
  }
  # else if(tolower(hazard) ==  "bspline"){
  #   if(frailty){
  #     nP0 <- nP01 + nP02 + nP03 + 1
  #     stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
  #     if(tolower(model)=="semi-markov"){
  #       nll <- nlogLikBS_ID_frail_SM(para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
  #                                    G_X1= if(nP1>0) as.matrix(Xmat1) else matrix(0,nrow = n, ncol=1),
  #                                    G_X2= if(nP2>0) as.matrix(Xmat2) else matrix(0,nrow = n, ncol=1),
  #                                    G_X3= if(nP3>0) as.matrix(Xmat3) else matrix(0,nrow = n, ncol=1),
  #                                    nP1 = nP1, nP2 = nP2, nP3 = nP3,
  #                                    G_basis1=basis1, G_basis2=basis2, G_basis3=basis3, wts=rep(1,n))
  #     } else{
  #       stop("markov bspline not yet implemented")
  #     }
  #   } else{
  #     stop("non-frailty bspline model not yet implemented")
  #   }
  # }
  else if(tolower(hazard)=="piecewise"){
    if(frailty){
      nP0 <- nP01 + nP02 + nP03 + 1
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
      if(tolower(model)=="semi-markov"){
        nll <- nlogLikPW_ID_frail_SM(para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                     Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                     basis1=basis1, basis2=basis2, basis3=basis3)
      } else{
        nll <- nlogLikPW_ID_frail_M(para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1)
      }
    } else{
      stop("non-frailty piecewise constant model not yet implemented")
    }
  } else if(tolower(hazard)=="royston-parmar"){
    if(frailty){
      nP0 <- nP01 + nP02 + nP03 + 1
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
      if(tolower(model)=="semi-markov"){
        nll <- nlogLikRP_ID_frail_SM(para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                     X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                                     X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                                     X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0),
                                     basis1=basis1, basis2=basis2, basis3=basis3,
                                     dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
      } else{
        nll <- nlogLikRP_ID_frail_M(para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                                    X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                                    X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0),
                                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1 = basis3_y1,
                                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
      }
    } else {"non-frailty royston-parmar model not yet implemented"}
  } else{ stop("please choose hazard of 'weibull', 'bspline', or 'piecewise'")}

  return(nll)
}



#' Gradient of Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the gradient of the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#'
#' @inheritParams nll_func
#'
#' @return Returns numeric vector of same length as \code{para} with sum of gradient contributions
#'   for the negative log likelihood.
#' @export
ngrad_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                       hazard, frailty, model,
                       basis1, basis2, basis3, basis3_y1,
                       dbasis1, dbasis2, dbasis3){

  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  nP01 <- if(!is.null(basis1)) ncol(basis1) else 0
  nP02 <- if(!is.null(basis2)) ncol(basis2) else 0
  nP03 <- if(!is.null(basis3)) ncol(basis3) else 0
  n <- length(y1)

  if(tolower(hazard) == "weibull"){
    if(frailty){
      nP0 <- 7
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
      if(tolower(model)=="semi-markov"){
        ngrad <- ngradWB_ID_frail_SM(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                     X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                                     X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                                     X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0))
      } else{ #markov model
        ngrad <- ngradWB_ID_frail_M(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                                    X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                                    X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0))
      }
    } else{stop("non-frailty not yet implemented")}
  }
  # else if(tolower(hazard) == "bspline"){
  #   if(frailty){
  #     nP0 <- nP01 + nP02 + nP03 + 1
  #     stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
  #     if(tolower(model)=="semi-markov"){
  #       ngrad <- ngradBS_ID_frail_SM(para=para,y1=y1, y2=y2, delta1=delta1, delta2=delta2, wts = rep(1,n),
  #                                    nP1 = nP1, nP2 = nP2, nP3 = nP3,
  #                                    G_X1=if(nP1>0) as.matrix(Xmat1) else matrix(0,nrow = n, ncol=1),
  #                                    G_X2=if(nP2>0) as.matrix(Xmat2) else matrix(0,nrow = n, ncol=1),
  #                                    G_X3=if(nP3>0) as.matrix(Xmat3) else matrix(0,nrow = n, ncol=1),
  #                                    G_basis1 = basis1, G_basis2 = basis2, G_basis3 = basis3)
  #     } else{stop("markov bspline model not yet implemented")}
  #   } else{
  #     stop("non-frailty bspline model not yet implemented.")
  #   }
  # }
  else if(tolower(hazard) == "piecewise"){
    if(frailty){
      nP0 <- nP01 + nP02 + nP03 + 1
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
      if(tolower(model)=="semi-markov"){
        ngrad <- ngradPW_ID_frail_SM(para=para,y1=y1, y2=y2, delta1=delta1, delta2=delta2,  #knots list is list with 1 or 3 elements, containing vector of internal knot locations
                                     Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                     basis1=basis1, basis2=basis2, basis3=basis3)
      } else{
        ngrad <- ngradPW_ID_frail_M(para=para,y1=y1, y2=y2, delta1=delta1, delta2=delta2,  #knots list is list with 1 or 3 elements, containing vector of internal knot locations
                                    Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1)
      }
    } else{
      stop("non-frailty piecewise model not yet implemented.")
    }
  } else if(tolower(hazard)=="royston-parmar"){
    if(frailty){
      nP0 <- nP01 + nP02 + nP03 + 1
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
      if(tolower(model)=="semi-markov"){
        ngrad <- ngradRP_ID_frail_SM(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                     X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                                     X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                                     X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0),
                                     basis1=basis1, basis2=basis2, basis3=basis3,
                                     dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
      } else{
        ngrad <- ngradRP_ID_frail_M(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                                    X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                                    X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0),
                                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
      }
    } else {"non-frailty royston-parmar model not yet implemented"}
  } else {stop("please choose hazard of 'weibull', 'bspline', or 'piecewise'")}
  return( ngrad )
}




#' Hessian of Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the Hessian of the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#'
#' @inheritParams nll_func
#'
#' @return Returns numeric square matrix with dimensions the same length as \code{para}
#'   with sum of gradient contributions for the negative log likelihood.
#' @export
nhess_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                       frailty, hazard, model){

  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  if(tolower(hazard) == "weibull"){
    if(frailty){
      nP0 <- 7
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
      if(tolower(model)=="semi-markov"){
        nhess <- nhessWB_ID_frail_SM(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                     X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                                     X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                                     X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0))
      } else{ #markov model
        nhess <- nhessWB_ID_frail_M(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                                    X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                                    X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0))
      }
    } else{stop("non-frailty not yet implemented")}
  } else{stop("non-Weibull not yet implemented")}

  return(nhess)

}
