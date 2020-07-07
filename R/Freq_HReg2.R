#' Fit Parametric Frailty Illness-Death Model for Semi-Competing Risks Data
#'
#' @inheritParams FreqID_HReg_R
#' @param hessian Boolean indicating whether the hessian (aka, the inverse covariance matrix)
#'   should be computed and returned.
#' @param control a list of control attributes passed directly into the \code{optim} function.
#'
#' @return A list.
#' @export
FreqID_HReg2 <- function(Formula, data, na.action="na.fail", subset=NULL,
                        hazard=c("weibull"),frailty=TRUE, model, knots_list = NULL,
                        startVals=NULL, hessian=TRUE, control=NULL){

  ####WRAPPER FUNCTION FOR REGULARIZED ESTIMATION####
  ##***********************************************##

  #wrapper function for fitting a model for a single choice of lambdas, and taking the best of a number of different methods


  ##INITIALIZE OPTIONS##
  ##******************##
  #checks on the penalties and on 'a' happen in the underlying functions to minimize duplication

  na.action <- match.arg(na.action) #technically na.fail I think is the only one currently implemented
  if (na.action != "na.fail" & na.action != "na.omit") {
    stop("na.action should be either na.fail or na.omit")
  }

  con=list(maxit=500)
  nmsC <- names(con)
  namc <- names(control)
  con[namc] <- control


  ##DATA PREPROCESSING##
  ##******************##
  ##MAKE THIS MORE EFFICIENT BY GOING DIRECTLY INTO NUMERICS

  #rearrange input Formula object (which stores the different pieces of the input formula)
  #this seems to ensure it's in the correct form, doesn't seem to make a diff
  form2 <- Formula::as.Formula(paste(Formula[2], Formula[1], Formula[3],
                                     sep=""))
  #arrange data into correct form, extracting subset of variables
  #and observations needed in model
  data <- stats::model.frame(form2, data=data, na.action=na.action,
                             subset=subset)
  #create matrices storing two outcomes, and then component vectors
  time1 <- Formula::model.part(form2, data=data, lhs=1)
  time2 <- Formula::model.part(form2, data=data, lhs=2)
  # Y <- cbind(time1[1], time1[2], time2[1], time2[2])
  y1 <- time1[[1]]
  delta1 <- time1[[2]]
  y2 <- time2[[1]]
  delta2 <- time2[[2]]
  #Create covariate matrices for each of three transition hazards
  Xmat1 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=1),
                                        data=data))
  Xmat2 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=2),
                                        data=data))
  Xmat3 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=3),
                                        data=data))

  ##PREPARE KNOTS AND BASIS FUNCTIONS FOR FLEXIBLE MODELS##
  ##*****************************************************##

  if(tolower(hazard) %in% c("bspline","royston-parmar","piecewise","pw","rp","bs")){
    if(is.null(knots_list)){
      p01 <- p02 <- p03 <- if(tolower(hazard) %in% c("bspline","bs")) 5 else 3
      knots_list <- get_default_knots_list(y1,y2,delta1,delta2,p01,p02,p03,hazard,model)
    }
    basis1 <- get_basis(x = y1,knots = knots_list[[1]],hazard = hazard)
    basis2 <- get_basis(x = y1,knots = knots_list[[2]],hazard = hazard)
    dbasis1 <- get_basis(x = y1,knots = knots_list[[1]],hazard = hazard,deriv = TRUE)
    dbasis2 <- get_basis(x = y1,knots = knots_list[[2]],hazard = hazard,deriv = TRUE)
    if(tolower(model)=="semi-markov"){
      basis3 <- get_basis(x = y2-y1,knots = knots_list[[3]],hazard = hazard)
      basis3_y1 <- NULL
      dbasis3 <- get_basis(x = y2-y1,knots = knots_list[[3]],hazard = hazard,deriv = TRUE)
    } else{
      basis3 <- get_basis(x = y2,knots = knots_list[[3]],hazard = hazard)
      basis3_y1 <- get_basis(x = y1,knots = knots_list[[3]],hazard = hazard)
      dbasis3 <- get_basis(x = y2,knots = knots_list[[3]],hazard = hazard,deriv = TRUE)
    }
    p01 <- ncol(basis1)
    p02 <- ncol(basis2)
    p03 <- ncol(basis3)
  } else{
    basis1 <- basis2 <- basis3 <- basis3_y1 <- dbasis1 <- dbasis2 <- dbasis3 <- NULL
    p01 <- p02 <- p03 <- 2 #must be weibull
  }


  if(is.null(startVals)){
    startVals <- get_start(y1=y1,y2=y2,delta1=delta1,delta2=delta2,Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,
                           basis1=basis1,basis2=basis2,basis3=basis3,basis3_y1=basis3_y1,
                           dbasis1=dbasis1,dbasis2=dbasis2,dbasis3=dbasis3,
                           hazard=hazard,frailty=frailty,model=model)
  }

  ##GET MLE##
  ##*******##

  fit0 <- stats::optim(par = startVals, fn = nll_func, gr = ngrad_func,
                              y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                              Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                              hazard=hazard,frailty=frailty,model=model,
                              basis1 = basis1, basis2 = basis2, basis3 = basis3, basis3_y1 = basis3_y1,
                              dbasis1 = dbasis1, dbasis2 = dbasis2, dbasis3 = basis3,
                              control=con, hessian=hessian,
                              method = if(tolower(hazard) %in% c("royston-parmar","rp")) "BFGS" else "L-BFGS")$par

  if (fit0$code == 0 | fit0$code == 1) {
    myLabels <- names(startVals)
    # myLabels <- c("log(kappa1)", "log(alpha1)", "log(kappa2)",
    #               "log(alpha2)", "log(kappa3)", "log(alpha3)")
    # if (frailty == TRUE)
    #   myLabels <- c(myLabels, "log(theta)")
    # myLabels <- c(myLabels, colnames(Xmat1), colnames(Xmat2),
    #               colnames(Xmat3))
    nP <- c(ncol(Xmat1), ncol(Xmat2), ncol(Xmat3))
    value <- list(estimate = fit0$par,
                  Finv = if(hessian) solve(fit0$hessian) else NA,
                  logLike = -fit0$value,
                  knots_list=knots_list,
                  myLabels = myLabels,
                  frailty = frailty,
                  nP = nP)#, Xmat = list(Xmat1, Xmat2, Xmat3))
    value$class <- c("Freq_HReg",
                     "ID",
                     "Ind",
                     switch(tolower(hazard),
                            weibull="WB",
                            wb="WB",
                            bspline="BS",
                            bs="BS",
                            "royston-parmar"="RP",
                            rp="RP",
                            piecewise="PW",
                            pw="PW"),
                     switch(tolower(model),
                            "semi-markov"="semi-Markov",
                            "markov"="Markov"))
    class(value) <- "Freq_HReg"
    return(value)
  }

}
