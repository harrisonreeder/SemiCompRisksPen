#' Fit Parametric Frailty Illness-Death Model for Semi-Competing Risks Data
#'
#' @inheritParams FreqID_HReg_R
#' @param p0_vec vector of length three of integers indicating how many baseline hazard parameters
#'   should be specified for each of the three transition hazards. This input is only relevant when
#'   hazard is something other than "weibull" and is superceded by knots_list.
#' @param hessian Boolean indicating whether the hessian (aka, the inverse of the covariance matrix)
#'   should be computed and returned.
#' @param control a list of control attributes passed directly into the \code{optim} function.
#' @param optim_method a string naming which \code{optim} method should be used.
#'
#' @return \code{FreqID_HReg2} returns an object of class \code{Freq_HReg}.
#' @export
FreqID_HReg2 <- function(Formula, data, na.action="na.fail", subset=NULL,
                        hazard=c("weibull"),frailty=TRUE, model, knots_list = NULL,
                        p0_vec = rep(4,3),
                        startVals=NULL, hessian=TRUE, control=NULL,
                        optim_method = if(tolower(hazard) %in% c("royston-parmar","rp")) "BFGS" else "L-BFGS-B"){

  #browser()

  ##Check that chosen hazard is among available options
  if(!(tolower(hazard) %in% c("weibull","royston-parmar","piecewise","wb","pw","rp"))){
    stop("valid choices of hazard are 'weibull', 'royston-parmar', or 'piecewise'")
  }

  ##INITIALIZE OPTIONS##
  ##******************##

  na.action <- match.arg(na.action) #technically na.fail I think is the only one currently implemented
  if (na.action != "na.fail" & na.action != "na.omit") {
    stop("na.action should be either na.fail or na.omit")
  }

  con=list(maxit=1000)
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
      if(length(p0_vec) != 3){stop("If hazard not equal to 'weibull' and knots_list set to NULL, then p0_vec must be vector of 3 integers.")}
      p01 <- p0_vec[1]
      p02 <- p0_vec[2]
      p03 <- p0_vec[3]
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
  fit0 <- tryCatch(stats::optim(par = startVals, fn = nll_func, gr = ngrad_func,
                              y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                              Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                              hazard=hazard,frailty=frailty,model=model,
                              basis1 = basis1, basis2 = basis2, basis3 = basis3, basis3_y1 = basis3_y1,
                              dbasis1 = dbasis1, dbasis2 = dbasis2, dbasis3 = dbasis3,
                              control=con, hessian=hessian,
                              method = optim_method),
                   error=function(cnd){
                     message(cnd)
                     cat("\n")
                     return(list(fail=TRUE,hazard=hazard,frailty=frailty,model=model,
                                 startVals=startVals,knots_list=knots_list,
                                 basis1=basis1,basis2=basis2,basis3=basis3,
                                 dbasis1=dbasis1,dbasis2=dbasis2,dbasis3=dbasis3,
                                 control=control))
                     })

  if(!is.null(fit0$fail)){
    return(fit0)
    }

  if (fit0$convergence == 0 | fit0$convergence == 1) {
    # myLabels <- names(startVals)
    # myLabels <- gsub(pattern = "_1|_2|_3",replacement = "",x = myLabels)
    if(tolower(hazard) %in% c("weibull","wb")){
      myLabels <- c("log(kappa1)", "log(alpha1)", "log(kappa2)",
                    "log(alpha2)", "log(kappa3)", "log(alpha3)")
    } else{
      myLabels <- c(paste0("phi1",1:p01),paste0("phi2",1:p02),paste0("phi3",1:p03))
    }
    if (frailty == TRUE){
      myLabels <- c(myLabels, "log(theta)")
    }
    myLabels <- c(myLabels, colnames(Xmat1), colnames(Xmat2),
                  colnames(Xmat3))
    nP <- c(ncol(Xmat1), ncol(Xmat2), ncol(Xmat3))
    nP0 <- c(p01,p02,p03)
    value <- list(estimate = fit0$par,
                  Finv = if(hessian) solve(fit0$hessian) else NA,
                  logLike = -fit0$value,
                  knots_list=knots_list,
                  myLabels = myLabels,
                  frailty = frailty,
                  optim_details = list(counts=fit0$counts,convergence=fit0$convergence,message=fit0$message,startVals=startVals),
                  nP = nP, nP0=nP0)#, Xmat = list(Xmat1, Xmat2, Xmat3))
    value$class <- c("Freq_HReg2",
                     "ID",
                     "Ind",
                     switch(tolower(hazard),
                            weibull="Weibull",
                            wb="Weibull",
                            bspline="B-Spline",
                            bs="B-Spline",
                            "royston-parmar"="Royston-Parmar",
                            rp="Royston-Parmar",
                            piecewise="Piecewise Constant",
                            pw="Piecewise Constant"),
                     switch(tolower(model),
                            "semi-markov"="semi-Markov",
                            "markov"="Markov"))
    class(value) <- "Freq_HReg2"
    return(value)
  }
}





print.Freq_HReg2 <- function (x, digits = 3, alpha = 0.05, ...)
{
  conf.level = alpha
  obj <- x
  logEst <- obj$estimate
  logSE <- sqrt(diag(obj$Finv))
  value <- cbind(logEst, logSE, logEst - abs(stats::qnorm(conf.level/2,
                                                   0, 1)) * logSE, logEst + abs(stats::qnorm(conf.level/2, 0, 1)) *
                   logSE)
  dimnames(value) <- list(obj$myLabels, c("Estimate", "SE",
                                          "LL", "UL"))
  if (obj$class[2] == "Surv") {
    cat("\nAnalysis of independent univariate time-to-event data \n")
    if (length(obj$myLabels) >= 3) {
      cat("Confidence level: ", conf.level, "\n", sep = "")
      cat("\nRegression coefficients:\n")
      cat(round(value[-c(1:2), ], digits = digits))
    }
  }
  if (obj$class[2] == "ID") {
    cat("\nAnalysis of independent semi-competing risks data \n")
    cat(obj$class[4], "baseline hazard specification\n")
    cat(obj$class[5], "specification for h3\n")
    value_theta <- matrix(exp(value[7, ]), ncol = 4)
    dimnames(value_theta) <- list("", c("Estimate", "SE",
                                        "LL", "UL"))
    value_theta[1, 2] <- value[7, 2] * exp(value[7, 1])
    cat("Confidence level: ", conf.level, "\n", sep = "")
    cat("\nVariance of frailties, theta:\n")
    if (obj$frailty == TRUE)
      print(round(value_theta, digits = digits))
    if (obj$frailty == FALSE)
      cat("NA")
    if (sum(obj$nP) != 0) {
      cat("\nRegression coefficients:\n")
      if (obj$frailty == TRUE)
        print(round(value[-c(1:(sum(obj$nP0)+1)), ], digits = digits))
      if (obj$frailty == FALSE)
        print(round(value[-c(1:(sum(obj$nP0))), ], digits = digits))
    }
    cat("\nNote: Covariates are arranged in order of transition number, 1->3.\n")
  }
  invisible()
}


summary.Freq_HReg2 <- function (object, digits = 3, alpha = 0.05, ...)
{
  # browser()
  conf.level = alpha
  obj <- object
  logEst <- obj$estimate
  logSE <- sqrt(diag(obj$Finv))
  results <- cbind(logEst, logEst - abs(stats::qnorm(conf.level/2,
                                              0, 1)) * logSE, logEst + abs(stats::qnorm(conf.level/2, 0, 1)) *
                     logSE)
  if (obj$class[2] == "Surv") { #this is unchanged
    if (length(obj$myLabels) >= 3) {
      output.coef <- results[-c(1:2), ]
      dimnames(output.coef) <- list(unique(obj$myLabels[-c(1:2)]),
                                    c("beta", "LL", "UL"))
    }
    else {
      output.coef <- NULL
    }
    output.h0 <- results[c(1:2), ]
    dimnames(output.h0) <- list(c("Weibull: log-kappa", "Weibull: log-alpha"),
                                c("h-PM", "LL", "UL"))
    value <- list(coef = output.coef, h0 = output.h0, code = obj$code,
                  logLike = obj$logLike, nP = nrow(results), class = obj$class)
  }
  if (obj$class[2] == "ID") {
    nP.0 <- ifelse(obj$frailty, sum(obj$nP0)+1, sum(obj$nP0))
    nP.1 <- obj$nP[1]
    nP.2 <- obj$nP[2]
    nP.3 <- obj$nP[3]
    beta.names <- unique(obj$myLabels[-c(1:nP.0)])
    nP <- length(beta.names)
    output <- matrix(NA, nrow = nP, ncol = 9)
    dimnames(output) <- list(beta.names, c("beta1", "LL",
                                           "UL", "beta2", "LL", "UL", "beta3", "LL", "UL"))
    for (i in 1:nP) {
      if (nP.1 != 0) {
        for (j in 1:nP.1) if (obj$myLabels[nP.0 + j] ==
                              beta.names[i])
          output[i, 1:3] <- results[nP.0 + j, ]
      }
      if (nP.2 != 0) {
        for (j in 1:nP.2) if (obj$myLabels[nP.0 + nP.1 +
                                           j] == beta.names[i])
          output[i, 4:6] <- results[nP.0 + nP.1 + j,
          ]
      }
      if (nP.3 != 0) {
        for (j in 1:nP.3) if (obj$myLabels[nP.0 + nP.1 +
                                           nP.2 + j] == beta.names[i])
          output[i, 7:9] <- results[nP.0 + nP.1 + nP.2 +
                                      j, ]
      }
    }
    output.coef <- output
    output <- matrix(NA, nrow = 1, ncol = 3)
    dimnames(output) <- list(c("theta"), c("Estimate", "LL",
                                           "UL"))
    if (obj$frailty == TRUE)
      output[1, ] <- exp(results[nP.0, ])
    if (obj$frailty == FALSE)
      output[1, ] <- rep(NA, 3)
    output.theta <- output

    knots_mat <- NULL
    if(obj$class[4]=="Weibull"){
      output <- matrix(NA, nrow = 2, ncol = 9)
      dimnames(output) <- list(c("Weibull: log-kappa", "Weibull: log-alpha"),
                               c("h1-PM", "LL", "UL", "h2-PM", "LL", "UL", "h3-PM",
                                 "LL", "UL"))
      output[1, 1:3] <- results[1, ]
      output[1, 4:6] <- results[3, ]
      output[1, 7:9] <- results[5, ]
      output[2, 1:3] <- results[2, ]
      output[2, 4:6] <- results[4, ]
      output[2, 7:9] <- results[6, ]
    } else{
      p01 <- obj$nP0[1]; p02 <- obj$nP0[2]; p03 <- obj$nP0[3]
      p0max <- max(obj$nP0)

      #generate "wide" matrix of baseline parameters by padding with 0s so all are same height
      output <- cbind(
        rbind(results[1:p01,],matrix(data=0,ncol=3,nrow=(p0max-p01))),
        rbind(results[(1+p01):(p01+p02),],matrix(data=0,ncol=3,nrow=(p0max-p02))),
        rbind(results[(1+p01):(p01+p02),],matrix(data=0,ncol=3,nrow=(p0max-p03)))
      )
      dimnames(output) <- list(paste0(obj$class[4],": phi",1:p0max),
                               c("h1-PM", "LL", "UL", "h2-PM", "LL", "UL", "h3-PM",
                                 "LL", "UL"))

      #lastly, make a matrix with the knot locations, padded with NAs
      knots_mat <- sapply(obj$knots_list,FUN = function(x){c(x,rep(NA,p0max-length(x)))})
      dimnames(knots_mat) <- list(paste0("knot",1:p0max),
                               c("h1","h2","h3"))
    }

    output.h0 <- output
    value <- list(coef = output.coef, theta = output.theta,
                  h0 = output.h0, code = obj$code, logLike = obj$logLike,
                  nP = nrow(results), class = obj$class, conf.level = conf.level,
                  knots_mat=knots_mat) #important to include the knots in the output
  }
  class(value) <- "summ.Freq_HReg2"
  return(value)
}



print.summ.Freq_HReg2 <- function (x, digits = 3, ...)
{
  obj <- x
  if (obj$class[2] == "Surv") {
    cat("\nAnalysis of independent univariate time-to-event data \n")
  }
  if (obj$class[2] == "ID") {
    cat("\nAnalysis of independent semi-competing risks data \n")
    cat(obj$class[4], "baseline hazard specification\n")
    cat(obj$class[5], "specification for h3\n")
  }
  cat("Confidence level: ", x$conf.level, "\n", sep = "")
  if (!is.null(obj$coef)) {
    cat("\nHazard ratios:\n")
    print(round(exp(obj$coef), digits = digits))
  }
  if (obj$class[2] == "ID") {
    cat("\nVariance of frailties:\n")
    print(round(obj$theta, digits = digits))
  }
  cat("\nBaseline hazard function components:\n")
  print(round(obj$h0, digits = digits))
  if (obj$class[4] != "Weibull") {
    cat("\nKnots:\n")
    print(round(obj$knots_mat, digits = digits))
  }
  invisible()
}



# These are prediction functions Kyu Ha has written in his current package, that I'll try and update
# For the new models.
#
#
# predict.Freq_HReg2 <- function (object, xnew = NULL, x1new = NULL, x2new = NULL, x3new = NULL,
#           tseq = c(0, 5, 10), alpha = 0.05, ...)
# {
#   conf.level = alpha
#   obj <- object
#   T2seq <- tseq
#   yLim <- NULL
#   if (obj$class[2] == "Surv") {
#     if (!is.null(x1new)) {
#       stop("'x1new' is for semi-competing risks models so it must be specified as NULL for univariate models")
#     }
#     if (!is.null(x2new)) {
#       stop("'x2new' is for semi-competing risks models so it must be specified as NULL for univariate models")
#     }
#     if (!is.null(x3new)) {
#       stop("'x3new' is for semi-competing risks models so it must be specified as NULL for univariate models")
#     }
#     T2 <- seq(from = min(T2seq), to = max(T2seq), length = 100)
#     kappa <- exp(obj$estimate[1])
#     alpha <- exp(obj$estimate[2])
#     log_kappa <- obj$estimate[1]
#     log_alpha <- obj$estimate[2]
#     S0 <- exp(-(kappa * (T2)^alpha))
#     if (!is.null(xnew)) {
#       J <- cbind(1, exp(log_alpha) * log(T2), matrix(xnew,
#                                                      nrow = 100, ncol = length(xnew), byrow = T))
#       Var.loglogS0 <- J %*% obj$Finv %*% t(J)
#     }
#     else {
#       J <- cbind(1, exp(log_alpha) * log(T2))
#       Var.loglogS0 <- J %*% obj$Finv[1:2, 1:2] %*% t(J)
#     }
#     se.loglogS0 <- sqrt(diag(Var.loglogS0))
#     se.loglogS0[is.na(se.loglogS0)] <- 0
#     LL <- S0^exp(-qnorm(conf.level/2) * se.loglogS0)
#     UL <- S0^exp(qnorm(conf.level/2) * se.loglogS0)
#     BS_tbl <- data.frame(time = T2, S = S0, LL = LL, UL = UL)
#     h0 <- alpha * kappa * (T2)^(alpha - 1)
#     if (!is.null(xnew)) {
#       J <- cbind(h0, h0 * (1 + alpha * log(T2)), h0 * matrix(xnew,
#                                                              nrow = 100, ncol = length(xnew), byrow = T))
#       Var.h0 <- J %*% obj$Finv %*% t(J)
#     }
#     else {
#       J <- cbind(h0, h0 * (1 + alpha * log(T2)))
#       Var.h0 <- J %*% obj$Finv[1:2, 1:2] %*% t(J)
#     }
#     se.h0 <- sqrt(diag(Var.h0))
#     se.h0[is.nan(se.h0)] <- 0
#     LLh0 <- h0 - qnorm(conf.level/2) * se.h0
#     ULh0 <- h0 + qnorm(conf.level/2) * se.h0
#     LLh0[LLh0 < 0] <- 0
#     T2h <- T2
#     if (T2[1] == 0) {
#       T2h <- T2h[-1]
#       h0 <- h0[-1]
#       LLh0 <- LLh0[-1]
#       ULh0 <- ULh0[-1]
#     }
#     BH_tbl <- data.frame(time = T2h, h = h0, LL = LLh0, UL = ULh0)
#     value <- list(h = BH_tbl, S = BS_tbl)
#   }
#   if (obj$class[2] == "ID") {
#     if (!is.null(xnew)) {
#       stop("'xnew' is for univariate models so it must be specified as NULL for semi-competing risks models")
#     }
#     nP = obj$nP
#     T2 <- seq(from = min(T2seq), to = max(T2seq), length = 100)
#     kappa <- exp(obj$estimate[1])
#     alpha <- exp(obj$estimate[2])
#     log_alpha <- obj$estimate[2]
#     S0.1 <- exp(-kappa * (T2)^alpha)
#     if (!is.null(x1new) & nP[1] > 0) {
#       J <- cbind(1, exp(log_alpha) * log(T2), matrix(x1new,
#                                                      nrow = 100, ncol = length(x1new), byrow = T))
#       if (obj$frailty == TRUE) {
#         Var.loglogS0 <- J %*% obj$Finv[c(1:2, 8:(8 +
#                                                    nP[1] - 1)), c(1:2, 8:(8 + nP[1] - 1))] %*%
#           t(J)
#       }
#       else if (obj$frailty == FALSE) {
#         Var.loglogS0 <- J %*% obj$Finv[c(1:2, 7:(7 +
#                                                    nP[1] - 1)), c(1:2, 7:(7 + nP[1] - 1))] %*%
#           t(J)
#       }
#     }
#     else if (is.null(x1new) | nP[1] == 0) {
#       J <- cbind(1, exp(log_alpha) * log(T2))
#       Var.loglogS0 <- J %*% obj$Finv[1:2, 1:2] %*% t(J)
#     }
#     se.loglogS0 <- sqrt(diag(Var.loglogS0))
#     LL.1 <- S0.1^exp(-qnorm(conf.level/2) * se.loglogS0)
#     UL.1 <- S0.1^exp(qnorm(conf.level/2) * se.loglogS0)
#     h0.1 <- alpha * kappa * (T2)^(alpha - 1)
#     if (!is.null(x1new) & nP[1] > 0) {
#       J <- cbind(h0.1, h0.1 * (1 + alpha * log(T2)), h0.1 *
#                    matrix(x1new, nrow = 100, ncol = length(x1new),
#                           byrow = T))
#       if (obj$frailty == TRUE) {
#         Var.h0.1 <- J %*% obj$Finv[c(1:2, 8:(8 + nP[1] -
#                                                1)), c(1:2, 8:(8 + nP[1] - 1))] %*% t(J)
#       }
#       else if (obj$frailty == FALSE) {
#         Var.h0.1 <- J %*% obj$Finv[c(1:2, 7:(7 + nP[1] -
#                                                1)), c(1:2, 7:(7 + nP[1] - 1))] %*% t(J)
#       }
#     }
#     else if (is.null(x1new) | nP[1] == 0) {
#       J <- cbind(h0.1, h0.1 * (1 + alpha * log(T2)))
#       Var.h0.1 <- J %*% obj$Finv[1:2, 1:2] %*% t(J)
#     }
#     se.h0.1 <- sqrt(diag(Var.h0.1))
#     se.h0.1[is.nan(se.h0.1)] <- 0
#     LLh0.1 <- h0.1 - qnorm(conf.level/2) * se.h0.1
#     ULh0.1 <- h0.1 + qnorm(conf.level/2) * se.h0.1
#     LLh0.1[LLh0.1 < 0] <- 0
#     kappa <- exp(obj$estimate[3])
#     alpha <- exp(obj$estimate[4])
#     log_alpha <- obj$estimate[4]
#     S0.2 <- exp(-kappa * (T2)^alpha)
#     if (!is.null(x2new) & nP[2] > 0) {
#       J <- cbind(1, exp(log_alpha) * log(T2), matrix(x2new,
#                                                      nrow = 100, ncol = length(x2new), byrow = T))
#       if (obj$frailty == TRUE) {
#         Var.loglogS0 <- J %*% obj$Finv[c(3:4, (8 + nP[1]):(8 +
#                                                              nP[1] + nP[2] - 1)), c(3:4, (8 + nP[1]):(8 +
#                                                                                                         nP[1] + nP[2] - 1))] %*% t(J)
#       }
#       else if (obj$frailty == FALSE) {
#         Var.loglogS0 <- J %*% obj$Finv[c(3:4, (7 + nP[1]):(7 +
#                                                              nP[1] + nP[2] - 1)), c(3:4, (7 + nP[1]):(7 +
#                                                                                                         nP[1] + nP[2] - 1))] %*% t(J)
#       }
#     }
#     else if (is.null(x2new) | nP[2] == 0) {
#       J <- cbind(1, exp(log_alpha) * log(T2))
#       Var.loglogS0 <- J %*% obj$Finv[3:4, 3:4] %*% t(J)
#     }
#     se.loglogS0 <- sqrt(diag(Var.loglogS0))
#     LL.2 <- S0.2^exp(-qnorm(conf.level/2) * se.loglogS0)
#     UL.2 <- S0.2^exp(qnorm(conf.level/2) * se.loglogS0)
#     h0.2 <- alpha * kappa * (T2)^(alpha - 1)
#     if (!is.null(x2new) & nP[2] > 0) {
#       J <- cbind(h0.2, h0.2 * (1 + alpha * log(T2)), h0.2 *
#                    matrix(x2new, nrow = 100, ncol = length(x2new),
#                           byrow = T))
#       if (obj$frailty == TRUE) {
#         Var.h0.2 <- J %*% obj$Finv[c(3:4, (8 + nP[1]):(8 +
#                                                          nP[1] + nP[2] - 1)), c(3:4, (8 + nP[1]):(8 +
#                                                                                                     nP[1] + nP[2] - 1))] %*% t(J)
#       }
#       else if (obj$frailty == FALSE) {
#         Var.h0.2 <- J %*% obj$Finv[c(3:4, (7 + nP[1]):(7 +
#                                                          nP[1] + nP[2] - 1)), c(3:4, (7 + nP[1]):(7 +
#                                                                                                     nP[1] + nP[2] - 1))] %*% t(J)
#       }
#     }
#     else if (is.null(x2new) | nP[2] == 0) {
#       J <- cbind(h0.2, h0.2 * (1 + alpha * log(T2)))
#       Var.h0.2 <- J %*% obj$Finv[3:4, 3:4] %*% t(J)
#     }
#     se.h0.2 <- sqrt(diag(Var.h0.2))
#     se.h0.2[is.nan(se.h0.2)] <- 0
#     LLh0.2 <- h0.2 - qnorm(conf.level/2) * se.h0.2
#     ULh0.2 <- h0.2 + qnorm(conf.level/2) * se.h0.2
#     LLh0.2[LLh0.2 < 0] <- 0
#     kappa <- exp(obj$estimate[5])
#     alpha <- exp(obj$estimate[6])
#     log_alpha <- obj$estimate[6]
#     S0.3 <- exp(-kappa * (T2)^alpha)
#     if (!is.null(x3new) & nP[3] > 0) {
#       J <- cbind(1, exp(log_alpha) * log(T2), matrix(x3new,
#                                                      nrow = 100, ncol = length(x3new), byrow = T))
#       if (obj$frailty == TRUE) {
#         Var.loglogS0 <- J %*% obj$Finv[c(5:6, (8 + nP[1] +
#                                                  nP[2]):(8 + nP[1] + nP[2] + +nP[3] - 1)), c(5:6,
#                                                                                              (8 + nP[1] + nP[3]):(8 + nP[1] + nP[2] + nP[3] -
#                                                                                                                     1))] %*% t(J)
#       }
#       else if (obj$frailty == FALSE) {
#         Var.loglogS0 <- J %*% obj$Finv[c(5:6, (7 + nP[1] +
#                                                  nP[2]):(7 + nP[1] + nP[2] + +nP[3] - 1)), c(5:6,
#                                                                                              (7 + nP[1] + nP[3]):(7 + nP[1] + nP[2] + nP[3] -
#                                                                                                                     1))] %*% t(J)
#       }
#     }
#     else if (is.null(x3new) | nP[3] == 0) {
#       J <- cbind(1, exp(log_alpha) * log(T2))
#       Var.loglogS0 <- J %*% obj$Finv[5:6, 5:6] %*% t(J)
#     }
#     se.loglogS0 <- sqrt(diag(Var.loglogS0))
#     LL.3 <- S0.3^exp(-qnorm(conf.level/2) * se.loglogS0)
#     UL.3 <- S0.3^exp(qnorm(conf.level/2) * se.loglogS0)
#     h0.3 <- alpha * kappa * (T2)^(alpha - 1)
#     if (!is.null(x3new) & nP[3] > 0) {
#       J <- cbind(h0.3, h0.3 * (1 + alpha * log(T2)), h0.3 *
#                    matrix(x3new, nrow = 100, ncol = length(x3new),
#                           byrow = T))
#       if (obj$frailty == TRUE) {
#         Var.h0.3 <- J %*% obj$Finv[c(5:6, (8 + nP[1] +
#                                              nP[2]):(8 + nP[1] + nP[2] + +nP[3] - 1)), c(5:6,
#                                                                                          (8 + nP[1] + nP[3]):(8 + nP[1] + nP[2] + nP[3] -
#                                                                                                                 1))] %*% t(J)
#       }
#       else if (obj$frailty == FALSE) {
#         Var.h0.3 <- J %*% obj$Finv[c(5:6, (7 + nP[1] +
#                                              nP[2]):(7 + nP[1] + nP[2] + +nP[3] - 1)), c(5:6,
#                                                                                          (7 + nP[1] + nP[3]):(7 + nP[1] + nP[2] + nP[3] -
#                                                                                                                 1))] %*% t(J)
#       }
#     }
#     else if (is.null(x3new) | nP[3] == 0) {
#       J <- cbind(h0.3, h0.3 * (1 + alpha * log(T2)))
#       Var.h0.3 <- J %*% obj$Finv[5:6, 5:6] %*% t(J)
#     }
#     se.h0.3 <- sqrt(diag(Var.h0.3))
#     se.h0.3[is.nan(se.h0.3)] <- 0
#     LLh0.3 <- h0.3 - qnorm(conf.level/2) * se.h0.3
#     ULh0.3 <- h0.3 + qnorm(conf.level/2) * se.h0.3
#     LLh0.3[LLh0.3 < 0] <- 0
#     T2h <- T2
#     if (T2[1] == 0) {
#       T2h <- T2h[-1]
#       h0.1 <- h0.1[-1]
#       LLh0.1 <- LLh0.1[-1]
#       ULh0.1 <- ULh0.1[-1]
#       h0.2 <- h0.2[-1]
#       LLh0.2 <- LLh0.2[-1]
#       ULh0.2 <- ULh0.2[-1]
#       h0.3 <- h0.3[-1]
#       LLh0.3 <- LLh0.3[-1]
#       ULh0.3 <- ULh0.3[-1]
#     }
#     BH1_tbl <- data.frame(time = T2h, h.1 = h0.1, LL.1 = LLh0.1,
#                           UL.1 = ULh0.1)
#     BH2_tbl <- data.frame(time = T2h, h.2 = h0.2, LL.2 = LLh0.2,
#                           UL.2 = ULh0.2)
#     BH3_tbl <- data.frame(time = T2h, h.3 = h0.3, LL.3 = LLh0.3,
#                           UL.3 = ULh0.3)
#     BS1_tbl <- data.frame(time = T2, S.1 = S0.1, LL.1 = LL.1,
#                           UL.1 = UL.1)
#     BS2_tbl <- data.frame(time = T2, S.2 = S0.2, LL.2 = LL.2,
#                           UL.2 = UL.2)
#     BS3_tbl <- data.frame(time = T2, S.3 = S0.3, LL.3 = LL.3,
#                           UL.3 = UL.3)
#     value <- list(h.1 = BH1_tbl, h.2 = BH2_tbl, h.3 = BH3_tbl,
#                   S.1 = BS1_tbl, S.2 = BS2_tbl, S.3 = BS3_tbl)
#   }
#   value$xnew <- xnew
#   value$x1new <- x1new
#   value$x2new <- x2new
#   value$x3new <- x3new
#   value$tseq <- tseq
#   value$setup$model <- obj$setup$model
#   value$class <- obj$class
#   class(value) <- "pred.Freq_HReg2"
#   return(value)
# }
#
#
#
#
#
# plot.pred.Freq_HReg <- function (x, plot.est = "Haz", xlab = NULL, ylab = NULL, ...)
# {
#   obj <- x
#   T2seq <- x$tseq
#   yLim <- NULL
#   if (obj$class[2] == "Surv") {
#     if (is.null(yLim)) {
#       if (plot.est == "Surv") {
#         yLim <- seq(from = 0, to = 1, by = 0.2)
#       }
#       if (plot.est == "Haz") {
#         grid <- (max(obj$h$UL) - min(obj$h$LL))/5
#         yLim <- seq(from = min(obj$h$LL), to = max(obj$h$UL),
#                     by = grid)
#       }
#     }
#     if (is.null(ylab)) {
#       if (plot.est == "Surv") {
#         ylab <- "Survival"
#       }
#       if (plot.est == "Haz") {
#         ylab <- "Hazard"
#       }
#     }
#     if (is.null(xlab))
#       xlab <- "Time"
#     if (plot.est == "Surv") {
#       plot(range(T2seq), range(yLim), xlab = xlab, ylab = ylab,
#            type = "n", main = expression(paste("Estimated ",
#                                                S(t), "")), axes = FALSE)
#       axis(1, at = T2seq)
#       axis(2, at = yLim)
#       lines(obj$S$time, obj$S$S, col = "red", lwd = 3)
#       lines(obj$S$time, obj$S$LL, col = "red", lwd = 3,
#             lty = 3)
#       lines(obj$S$time, obj$S$UL, col = "red", lwd = 3,
#             lty = 3)
#     }
#     if (plot.est == "Haz") {
#       plot(range(T2seq), range(yLim), xlab = xlab, ylab = ylab,
#            type = "n", main = expression(paste("Estimated ",
#                                                h(t), "")), axes = FALSE)
#       axis(1, at = T2seq)
#       axis(2, at = round(yLim, 4))
#       lines(obj$h$time, obj$h$h, col = "red", lwd = 3)
#       lines(obj$h$time, obj$h$LL, col = "red", lwd = 3,
#             lty = 3)
#       lines(obj$h$time, obj$h$UL, col = "red", lwd = 3,
#             lty = 3)
#     }
#   }
#   if (obj$class[2] == "ID") {
#     if (is.null(ylab)) {
#       if (plot.est == "Surv") {
#         ylab <- "Survival"
#       }
#       if (plot.est == "Haz") {
#         ylab <- "Hazard"
#       }
#     }
#     if (is.null(xlab)) {
#       xlab <- c("Time", "Time", "Time")
#       if (obj$class[5] == "semi-Markov") {
#         xlab[3] <- "Time since non-terminal event"
#       }
#     }
#     if (is.null(yLim)) {
#       if (plot.est == "Surv") {
#         yLim <- seq(from = 0, to = 1, by = 0.2)
#       }
#       if (plot.est == "Haz") {
#         ygrid <- (max(x$h.1$UL.1, x$h.2$UL.2, x$h.3$UL.3) -
#                     0)/5
#         yLim <- seq(from = 0, to = max(x$h.1$UL.1, x$h.2$UL.2,
#                                        x$h.3$UL.3), by = ygrid)
#       }
#     }
#     if (plot.est == "Surv") {
#       par(mfrow = c(1, 3))
#       plot(range(T2seq), range(yLim), xlab = xlab[1], ylab = ylab,
#            type = "n", main = expression(paste("Estimated ",
#                                                S[1](t), "")), axes = FALSE)
#       axis(1, at = T2seq)
#       axis(2, at = yLim)
#       lines(obj$S.1$time, obj$S.1$S.1, col = "blue", lwd = 3)
#       lines(obj$S.1$time, obj$S.1$LL.1, col = "blue", lwd = 3,
#             lty = 3)
#       lines(obj$S.1$time, obj$S.1$UL.1, col = "blue", lwd = 3,
#             lty = 3)
#       plot(range(T2seq), range(yLim), xlab = xlab[2], ylab = ylab,
#            type = "n", main = expression(paste("Estimated ",
#                                                S[2](t), "")), axes = FALSE)
#       axis(1, at = T2seq)
#       axis(2, at = yLim)
#       lines(obj$S.2$time, obj$S.2$S.2, col = "red", lwd = 3)
#       lines(obj$S.2$time, obj$S.2$LL.2, col = "red", lwd = 3,
#             lty = 3)
#       lines(obj$S.2$time, obj$S.2$UL.2, col = "red", lwd = 3,
#             lty = 3)
#       plot(range(T2seq), range(yLim), xlab = xlab[3], ylab = ylab,
#            type = "n", main = expression(paste("Estimated ",
#                                                S[3](t), "")), axes = FALSE)
#       axis(1, at = T2seq)
#       axis(2, at = yLim)
#       lines(obj$S.3$time, obj$S.3$S.3, col = "red", lwd = 3)
#       lines(obj$S.3$time, obj$S.3$LL.3, col = "red", lwd = 3,
#             lty = 3)
#       lines(obj$S.3$time, obj$S.3$UL.3, col = "red", lwd = 3,
#             lty = 3)
#     }
#     if (plot.est == "Haz") {
#       par(mfrow = c(1, 3))
#       plot(range(T2seq), range(yLim), xlab = xlab[1], ylab = ylab,
#            type = "n", main = expression(paste("Estimated ",
#                                                h[1](t), "")), axes = FALSE)
#       axis(1, at = T2seq)
#       axis(2, at = round(yLim, 4))
#       lines(obj$h.1$time, obj$h.1$h.1, col = "blue", lwd = 3)
#       lines(obj$h.1$time, obj$h.1$LL.1, col = "blue", lwd = 3,
#             lty = 3)
#       lines(obj$h.1$time, obj$h.1$UL.1, col = "blue", lwd = 3,
#             lty = 3)
#       plot(range(T2seq), range(yLim), xlab = xlab[2], ylab = ylab,
#            type = "n", main = expression(paste("Estimated ",
#                                                h[2](t), "")), axes = FALSE)
#       axis(1, at = T2seq)
#       axis(2, at = round(yLim, 4))
#       lines(obj$h.2$time, obj$h.2$h.2, col = "red", lwd = 3)
#       lines(obj$h.2$time, obj$h.2$LL.2, col = "red", lwd = 3,
#             lty = 3)
#       lines(obj$h.2$time, obj$h.2$UL.2, col = "red", lwd = 3,
#             lty = 3)
#       plot(range(T2seq), range(yLim), xlab = xlab[3], ylab = ylab,
#            type = "n", main = expression(paste("Estimated ",
#                                                h[3](t), "")), axes = FALSE)
#       axis(1, at = T2seq)
#       axis(2, at = round(yLim, 4))
#       lines(obj$h.3$time, obj$h.3$h.3, col = "red", lwd = 3)
#       lines(obj$h.3$time, obj$h.3$LL.3, col = "red", lwd = 3,
#             lty = 3)
#       lines(obj$h.3$time, obj$h.3$UL.3, col = "red", lwd = 3,
#             lty = 3)
#     }
#   }
#   invisible()
# }
