####HELPER KNOT AND BASIS FUNCTIONS####

#' Generate List of Knot Location Vectors from Event Times
#'
#' This function creates a list containing three numeric vectors. Each numeric vector
#'   is a sequence of increasing integers giving the location of knots used for
#'   spline and piecewise baseline hazard specifications. This function
#'   generates this list according to the conventions and recommended locations
#'   of these knots, which depends on the choice of hazard specification, number
#'   of knots requested, and distribution of observed event times.
#'
#' @inheritParams proximal_gradient_descent
#' @param p01,p02,p03 Integers indicating how many parameters the model for each
#'   transition baseline hazard should specify.
#'
#' @return A list of three increasing sequences of integers, each corresponding to
#'   the knots for the flexible model on the corresponding transition baseline hazard.
#' @export
get_default_knots_list <- function(y1,y2,delta1,delta2,p01,p02,p03,hazard,model){

  if(tolower(hazard) %in% c("bspline","bs")){

    quantile_seq1 <- seq(from = 0,to = 1, length.out = p01-2)[-c(1,p01-2)]
    quantile_seq2 <- seq(from = 0,to = 1, length.out = p02-2)[-c(1,p02-2)]
    quantile_seq3 <- seq(from = 0,to = 1, length.out = p03-2)[-c(1,p03-2)]

    knots1 <- c(0,stats::quantile(y1[delta1==1],quantile_seq1),max(y1))
    knots2 <- c(0,stats::quantile(y1[(1-delta1)*delta2==1],quantile_seq2),max(y1))
    if(tolower(model)=="semi-markov"){
      knots3 <- c(0,stats::quantile((y2-y1)[delta1*delta2==1],quantile_seq3),max(y2-y1))
    } else {
      knots3 <- c(0,stats::quantile(y2[delta1*delta2==1],quantile_seq3),max(y2))
    }
  } else if(tolower(hazard) %in% c("piecewise","pw")){

    quantile_seq1 <- seq(from = 0,to = 1, length.out = p01+1)[-c(1,p01+1)]
    quantile_seq2 <- seq(from = 0,to = 1, length.out = p02+1)[-c(1,p02+1)]
    quantile_seq3 <- seq(from = 0,to = 1, length.out = p03+1)[-c(1,p03+1)]

    knots1 <- c(0,stats::quantile(y1[delta1==1],quantile_seq1))
    knots2 <- c(0,stats::quantile(y1[(1-delta1)*delta2==1],quantile_seq2))
    if(tolower(model)=="semi-markov"){
      knots3 <- c(0,stats::quantile((y2-y1)[delta1*delta2==1],quantile_seq3))
    } else {
      knots3 <- c(0,stats::quantile(y2[delta1*delta2==1],quantile_seq3))
    }

  } else if(tolower(hazard) %in% c("royston-parmar","rp")){

    quantile_seq1 <- seq(from = 0.05,to = 0.95, length.out = p01)[-c(1,p01)]
    quantile_seq2 <- seq(from = 0.05,to = 0.95, length.out = p02)[-c(1,p02)]
    quantile_seq3 <- seq(from = 0.05,to = 0.95, length.out = p03)[-c(1,p03)]

    knots1 <- c(0,stats::quantile(log(y1)[delta1==1],quantile_seq1),max(log(y1)))
    knots2 <- c(0,stats::quantile(log(y1)[(1-delta1)*delta2==1],quantile_seq2),max(log(y1)))
    if(tolower(model)=="semi-markov"){
      knots3 <- c(0,stats::quantile(log(y2-y1)[delta1*delta2==1],quantile_seq3),max(log(y2-y1)))
    } else {
      knots3 <- c(0,stats::quantile(log(y2)[delta1*delta2==1],quantile_seq3),max(y2))
    }
  } else {return(NULL)}

  return(list(knots1,knots2,knots3))

}

#' Get Basis Function Values for Flexible Hazard Specifications
#'
#' @param x Numeric vector of event times (e.g., \code{y1} or \code{y2}) at which
#'   to generate basis function values.
#' @param knots Increasing vector of integers corresponding to the knots used
#'   in the desired spline or piecewise specification. Often an element of
#'   list generated from \code{\link{get_default_knots_list}}.
#' @inheritParams proximal_gradient_descent
#' @param deriv Boolean for whether returned values should be from derivatives of
#'   basis functions if \code{TRUE}, or original basis functions if \code{FALSE}.
#'
#' @return A matrix with each row corresponding to an element of the input, and
#'   each column giving the corresponding basis function value.
#' @export
get_basis <- function(x,knots,hazard,deriv=FALSE){

  #the exact form of the knots passed into this function come from the above function


  if(tolower(hazard) %in% c("bspline","bs")){
    if(deriv){return(NULL)}
    basis_out <- splines::bs(x = x,knots = knots[-c(1,length(knots))],Boundary.knots = knots[c(1,length(knots))],intercept = TRUE)
  } else if(tolower(hazard) %in% c("piecewise","pw")){
    if(deriv){return(NULL)}
    basis_out <- pw_cum_mat(y = x,knots = knots)
  } else if(tolower(hazard) %in% c("royston-parmar","rp")){
    temp_log <- log(x)
    temp_log[is.infinite(temp_log)] <- NA
    if(deriv){
#      basis_out <- dbasis(x=temp_log,knots=knots)
      basis_out <- ns_d(x = temp_log,knots = knots[-c(1,length(knots))],Boundary.knots = knots[c(1,length(knots))],intercept = TRUE)
      basis_out[is.na(basis_out)] <- 1 #can't set this to 0, because it is then logged and that causes a mess even when it multiplies with delta1delta2 and would otherwise be 0
    } else{
#      basis_out <- basis(x = temp_log,knots = knots)
      basis_out <- splines::ns(x = temp_log,knots = knots[-c(1,length(knots))],Boundary.knots = knots[c(1,length(knots))],intercept = TRUE)
      basis_out[is.na(basis_out)] <- 0
    }
  } else {return(NULL)}

  return(basis_out)
}


#' Get Starting Parameter Values
#'
#' A function to generate principled starting values for optimization, based on
#'   model specifications.
#'
#' @inheritParams proximal_gradient_descent
#' @param sparse_start Boolean indicating whether to set all regression parameters
#'   to 0 if \code{TRUE}, or to pre-estimate them using univariate models if \code{FALSE}.
#'
#' @return A vector of starting parameter values.
#' @export
get_start <- function(y1,y2,delta1,delta2,
                      Xmat1,Xmat2,Xmat3,
                      basis1,basis2,basis3,basis3_y1,
                      dbasis1,dbasis2,dbasis3,
                      hazard,frailty,model,sparse_start=FALSE){
  #generate starting values based on the chosen form of the baseline hazard.

  # browser()
  #number of parameters in each arm dictated by number of covariate columns in each matrix
  p1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  p2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  p3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  #For all methods, start by fitting three univariate models. If this were too big, consider a ridge approach
  if(p1 == 0 | sparse_start){
    fit_survreg_1 <- survival::survreg(survival::Surv(y1,delta1) ~ 1,
                                       dist="weibull")
  } else{
    fit_survreg_1 <- survival::survreg(survival::Surv(y1,delta1) ~ Xmat1,
                                       dist="weibull")
  }
  if(p2 == 0 | sparse_start){
    fit_survreg_2 <- survival::survreg(survival::Surv(y2,delta2) ~ 1,
                                       dist="weibull")
  } else{
    fit_survreg_2 <- survival::survreg(survival::Surv(y2,delta2) ~ Xmat2,
                                       dist="weibull")
  }
  if(p3 == 0 | sparse_start){
    fit_survreg_3 <- survival::survreg(survival::Surv(y2,delta2) ~ 1,
                                       dist="weibull")
  } else{
    fit_survreg_3 <- survival::survreg(survival::Surv(y2,delta2) ~ Xmat3,
                                       dist="weibull")
  }

  alpha1 <- 1/fit_survreg_1$scale
  alpha2 <- 1/fit_survreg_2$scale
  alpha3 <- 1/fit_survreg_3$scale
  kappa1 <- exp(-alpha1 * stats::coef(fit_survreg_1)[1])
  kappa2 <- exp(-alpha2 * stats::coef(fit_survreg_2)[1])
  kappa3 <- exp(-alpha3 * stats::coef(fit_survreg_3)[1])
  if(sparse_start){
    beta1 <- numeric(p1)
    beta2 <- numeric(p2)
    beta3 <- numeric(p3)
  } else{
    beta1 <- -stats::coef(fit_survreg_1)[-1] * alpha1
    beta2 <- -stats::coef(fit_survreg_2)[-1] * alpha2
    beta3 <- -stats::coef(fit_survreg_3)[-1] * alpha3
  }

  #help create a useful naming convention (varname_1)
  if(p1 > 0){
    names(beta1) <- if(!is.null(colnames(Xmat1))) paste0(colnames(Xmat1),"_1") else NULL
  }
  if(p2 > 0){
    names(beta2) <- if(!is.null(colnames(Xmat2))) paste0(colnames(Xmat2),"_2") else NULL
  }
  if(p3 > 0){
    names(beta3) <- if(!is.null(colnames(Xmat3))) paste0(colnames(Xmat3),"_3") else NULL
  }

  #for weibull, basically carry over the same stuff from originally.
  if(tolower(hazard) %in% c("weibull","wb")){
    #assign starting values
    startVals <- c(log(kappa1), #h1 intercept (k1)
                   log(alpha1), #a1
                   log(kappa2), #h2 intercept (k2)
                   log(alpha2), #a2
                   log(kappa3), #h3 intercept (k3)
                   log(alpha3)) #a3
    names(startVals) <- c("lkappa1","lalpha1","lkappa2","lalpha2","lkappa3","lalpha3")
    if (frailty == TRUE) {
      startVals <- c(startVals, "ltheta"=log(0.5)) #ltheta starting value of log(0.5) just because
    }
    startVals <- c(startVals,
                   beta1, #h1 remaining covariates
                   beta2, #h2 remaining covariates
                   beta3) #h3 remaining covariates
    return(startVals)
  }

  p01<- ncol(basis1)
  p02<- ncol(basis2)
  p03<- ncol(basis3)

  if(tolower(hazard) %in% c("bspline", "bs", "piecewise", "pw")){
    #run the model with no covariates to get possible start values
    startVals <- stats::optim(par = rep(0,p01+p02+p03+1),fn = nll_func, gr = ngrad_func,
                       y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                       Xmat1=NULL, Xmat2=NULL, Xmat3=NULL,
                       hazard=hazard, frailty=frailty, model=model,
                       basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                       dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                       method="L-BFGS-B",control = list(maxit=300))$par
    names(startVals) <- c(paste0("phi",1:p01,"_1"),paste0("phi",1:p02,"_2"),paste0("phi",1:p03,"_3"),"ltheta")
    startVals <- c(startVals,beta1,beta2,beta3)
    return(startVals)
  }

  if(tolower(hazard) %in% c("royston-parmar","rp")){
    #following discussion in Royston-Parmar (2002), use weibull to fit log(H0) and then regress spline basis
    #generates starting estimates corresponding to weibull baseline hazard
    log_Haz01 <- log(kappa1) + alpha1 * log(y1[delta1==1])
    phi1 <- stats::lm.fit(x = basis1[delta1==1,], y = log_Haz01 )$coef

    log_Haz02 <- log(kappa2) + alpha2 * log(y1[(1-delta1)*delta2==1])
    phi2 <- stats::lm.fit(x = basis2[(1-delta1)*delta2==1,], y = log_Haz02)$coef

    if(tolower(model)=="semi-markov"){
      log_Haz03 <- log(kappa3) + alpha3 * log((y2-y1)[delta1*delta2==1])
      phi3 <- stats::lm.fit(x = basis3[delta1*delta2==1,], y = log_Haz03)$coef
    } else {
      log_Haz03 <- log(kappa3) + alpha3 * log(y2[delta1*delta2==1])
      phi3 <- stats::lm.fit(x = basis3[delta1*delta2==1,], y = log_Haz03)$coef
    }

    startVals <- c(phi1,phi2,phi3,log(0.5))
    names(startVals) <- c(paste0("phi",1:p01,"_1"),paste0("phi",1:p02,"_2"),paste0("phi",1:p03,"_3"),"ltheta")
    startVals <- c(startVals,beta1,beta2,beta3)

    #do a final check that the resulting values are feasible in the royston parmar model
    getout <- FALSE
    for(i in c(seq(0,2,0.5),seq(-0.5,-2,-0.5))){
      for(j in c(seq(0,2,0.5),seq(-0.5,-2,-0.5))){
        # for(k in c(seq(0,2,0.5),seq(-0.5,-2,-0.5))){
        #we will just consider toggling the first two parameter values
        startVals[c(1,1+p01,1+p01+p02)] <- startVals[c(1,1+p01,1+p01+p02)] + i
        startVals[c(2,2+p01,2+p01+p02)] <- startVals[c(2,2+p01,2+p01+p02)] + j
        ll_temp <- nll_func(para=startVals,y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                            Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                            hazard=hazard, frailty=frailty, model=model,
                            basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                            dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
        if(!is.nan(ll_temp) ){
          print('viable royston parmar starting value found!')
          getout <- TRUE
        }
        # }
        if(getout){break}
      }
      if(getout){break}
    }

    if(!getout){stop("starting royston-parmar value is not finite. Input valid starting values.")}
  }

  return(startVals)
}



####CHECKING FUNCTIONS####


#function to check that penalty parameters meet basic sanity checks
check_pen_params <- function(penalty,penalty_fusedcoef,penalty_fusedbaseline,
                             lambda,lambda_fusedcoef,lambda_fusedbaseline,
                             a){
  if(!(penalty %in% c("scad","mcp","adalasso","lasso"))){
    stop("penalty must be either 'scad', 'mcp', 'adalasso', 'lasso'")
  }
  if(!(penalty_fusedcoef %in% c("none","fusedlasso","adafusedlasso"))){
    stop("penalty_fusedcoef must be 'none', 'fusedlasso' or 'adafusedlasso'")
  }
  if(!(penalty_fusedbaseline %in% c("none","fusedlasso","adafusedlasso"))){
    stop("penalty_fusedbaseline must be 'none', 'fusedlasso' or 'adafusedlasso'")
  }

  if (a <= 0)
    stop("a must be greater than 0")
  if (a <= 1 & penalty == "mcp")
    stop("a must be greater than 1 for the MC penalty")
  if (a <= 2 & penalty == "scad")
    stop("a must be greater than 2 for the SCAD penalty")

  if (any(lambda < 0)){
    stop("lambda must be non-negative")
  }
  if (any(lambda_fusedcoef < 0)){
    stop("lambda must be non-negative")
  }
  if (any(lambda_fusedbaseline < 0)){
    stop("lambda must be non-negative")
  }

  invisible(penalty) #it's best practice to invisibly return first argument, according to hadley wickham
}
