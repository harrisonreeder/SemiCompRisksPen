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
                     return(list(fail=TRUE,
                                 formula = form2,hazard=hazard,frailty=frailty,model=model,
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
                  Finv = if(hessian) MASS::ginv(fit0$hessian) else NA,
                  logLike = -fit0$value,
                  knots_list=knots_list,
                  myLabels = myLabels,
                  frailty = frailty,
                  formula = form2,
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
  logSE <- if(all(is.na(obj$Finv))) NA else sqrt(diag(obj$Finv))
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
  logSE <- if(all(is.na(obj$Finv))) NA else sqrt(diag(obj$Finv))
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
    } else{ #this covers piecewise as the main other case
      p01 <- obj$nP0[1]; p02 <- obj$nP0[2]; p03 <- obj$nP0[3]
      p0max <- max(obj$nP0)

      #generate "wide" matrix of baseline parameters by padding with 0s so all are same height
      output <- cbind(
        rbind(results[1:p01,],matrix(data=0,ncol=3,nrow=(p0max-p01))),
        rbind(results[(1+p01):(p01+p02),],matrix(data=0,ncol=3,nrow=(p0max-p02))),
        rbind(results[(1+p01+p02):(p01+p02+p03),],matrix(data=0,ncol=3,nrow=(p0max-p03)))
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
predict.Freq_HReg2 <- function (object, xnew = NULL, x1new = NULL, x2new = NULL, x3new = NULL,
          tseq = c(0, 5, 10), alpha = 0.05, ...)
{
  # browser()
  conf.level = alpha
  obj <- object
  T2seq <- tseq
  yLim <- NULL
  if (obj$class[2] == "ID") {
    if (!is.null(xnew)) {
      stop("'xnew' is for univariate models so it must be specified as NULL for semi-competing risks models")
    }
    nP = obj$nP
    nP0 = obj$nP0
    nP0_tot <- if (obj$frailty == TRUE) sum(nP0) + 1 else sum(nP0)
    #WHY IS THIS LIKE THIS, AND NOT TSEQ ITSELF?
    # T2 <- seq(from = min(T2seq), to = max(T2seq), length = 100)
    T2 <- tseq

    eta1 <- if (!is.null(x1new) & nP[1] > 0) {
      x1new %*% obj$estimate[(1 + nP0_tot):(nP0_tot + nP[1])]
    } else 0
    eta2 <- if (!is.null(x2new) & nP[2] > 0) {
      x2new %*% obj$estimate[(1 + nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])]
    } else 0
    eta3 <- if (!is.null(x3new) & nP[3] > 0) {
      x3new %*% obj$estimate[(1 + nP0_tot + nP[1] + nP[2]):(nP0_tot + nP[1] + nP[2] + nP[3])]
    } else 0

    #FIRST TRANSITION
    if(obj$class[4] == "Weibull"){
      kappa <- exp(obj$estimate[1])
      alpha <- exp(obj$estimate[2])
      log_alpha <- obj$estimate[2]
      S.1 <- exp(-kappa * (T2)^alpha) * exp(eta1)
      #Use delta method to compute baseline survival confidence intervals
      #but, it would only be problematic at the exact knots I think?
      if(all(is.na(obj$Finv))){
        Var.loglogS.1 <- NA
      } else if (!is.null(x1new) & nP[1] > 0) {
        #matrix with as many columns as parameters in S1, and as many rows as times in T2
        #under weibull, log(-log(S(t))) = log_alpha + log_kappa + xtbeta + alpha*log(t)
        #so, each row of J is the gradient at a particular T2 wrt logalpha, logkappa, beta1, ..., betap
        J <- cbind(1, exp(log_alpha) * log(T2),
                   matrix(x1new, nrow = length(T2), ncol = length(x1new), byrow = T))
        #vector with as many elements as times in T2
        Var.loglogS.1 <- J %*% obj$Finv[c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1])),
                                      c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1]))] %*% t(J)
      } else if (is.null(x1new) | nP[1] == 0) {
        J <- cbind(1, exp(log_alpha) * log(T2))
        Var.loglogS.1 <- J %*% obj$Finv[1:nP0[1], 1:nP0[1]] %*% t(J)
      }

      h.1 <- alpha * kappa * (T2)^(alpha - 1) * exp(eta1)
      #hazard is exp(logalpha) * exp(logkappa) * t^(exp(logalpha)-1) * exp(xtbeta)
      #so, each row of J is gradient wrt logkappa, logalpha, beta1, ..., betap
      #WHY ISNT EVERYTHING MULTIPLIED BY exp(xtbeta) ? seems like it should be. I'm gonna add it
      if(all(is.na(obj$Finv))){
        Var.h.1 <- NA
      } else if (!is.null(x1new) & nP[1] > 0) {
        J <- cbind(h.1, h.1 * (1 + alpha * log(T2)), h.1 *
                     matrix(x1new, nrow = length(T2), ncol = length(x1new), byrow = T))
        Var.h.1 <- J %*% obj$Finv[c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1])),
                                  c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1]))] %*% t(J)
      }
      else if (is.null(x1new) | nP[1] == 0) {
        J <- cbind(h.1, h.1 * (1 + alpha * log(T2)))
        Var.h.1 <- J %*% obj$Finv[1:nP0[1], 1:nP0[1]] %*% t(J)
      }
    } else{ #for now the only alternative we entertain is piecewise constant
      stopifnot(obj$class[4] == "Piecewise Constant")
      basis1 <- get_basis(x = T2,knots = obj$knots_list[[1]],hazard = "piecewise")
      phi1 <- obj$estimate[1:nP0[1]]
      Lambda01 <- as.vector(basis1 %*% exp(phi1))
      S.1 <- exp(-Lambda01*exp(eta1))
      #Use delta method to compute baseline survival confidence intervals
      if(all(is.na(obj$Finv))){
        Var.loglogS.1 <- NA
      } else if (!is.null(x1new) & nP[1] > 0) {
        #matrix with as many columns as parameters in S1, and as many rows as times in T2
        #under piecewise, log(-log(S(t))) = log( basis1 %*% exp(phi1)) + log(xtbeta)
        #so, each row of J is the gradient at a particular T2 wrt phi1, ..., phiK, beta1, ..., betap
        J <- cbind(
          #first, this multiplies every row of basis1 by exp(phi1),
          #then, this divides every column of the resulting matrix by the vector Lambda01
          t(t(basis1) * exp(phi1)) / Lambda01,
                    matrix(x1new, nrow = length(T2), ncol = length(x1new), byrow = T))
        #vector with as many elements as times in T2
        Var.loglogS.1 <- J %*% obj$Finv[c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1])),
                                        c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1]))] %*%
          t(J)
      } else if (is.null(x1new) | nP[1] == 0) {
        J <- t(t(basis1) * exp(phi1)) / Lambda01
        Var.loglogS.1 <- J %*% obj$Finv[1:nP0[1], 1:nP0[1]] %*% t(J)
      }

      #Use delta method compute baseline hazard confidence intervals
      #vector saying which interval each time falls into
      cut_cats1 <- rowSums(basis1!=0)
      if(T2[1]==0){cut_cats1[1] <- 1}
      stopifnot(length(cut_cats1)==length(T2))
      h.1 <- exp(phi1)[cut_cats1] * exp(eta1)
      #build a matrix with 0's everywhere, except on ith row, set column that T2i falls in to 1
      temp_mat <- matrix(data=0,nrow=length(T2),ncol=nP0[1])
      temp_mat[cbind(1:length(T2),cut_cats1)] <- 1
      temp_mat <- t(t(temp_mat) * exp(phi1))

      if(all(is.na(obj$Finv))){
        Var.h.1 <- NA
      } else if (!is.null(x1new) & nP[1] > 0) {
        J <- cbind(temp_mat * exp(eta1), h.1 *
                     matrix(x1new, nrow = length(T2), ncol = length(x1new), byrow = T))
        Var.h.1 <- J %*% obj$Finv[c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1])),
                                  c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1]))] %*% t(J)
      }
      else if (is.null(x1new) | nP[1] == 0) {
        J <- temp_mat * exp(eta1)
        Var.h.1 <- J %*% obj$Finv[1:nP0[1], 1:nP0[1]] %*% t(J)
      }
    }
    if(all(is.na(Var.loglogS.1))){
      se.loglogS.1 <- NA
    } else{
      se.loglogS.1 <- sqrt(diag(Var.loglogS.1))
      se.loglogS.1[is.nan(se.loglogS.1)] <- 0
    }
    LL.1 <- S.1^exp(-qnorm(conf.level/2) * se.loglogS.1)
    UL.1 <- S.1^exp(qnorm(conf.level/2) * se.loglogS.1)
    if(all(is.na(Var.h.1))){
      se.h.1 <- NA
    } else{
      se.h.1 <- sqrt(diag(Var.h.1))
      se.h.1[is.nan(se.h.1)] <- 0
    }
    LLh.1 <- h.1 + qnorm(conf.level/2) * se.h.1 #sign reversed because 0.025 quantile is negative
    ULh.1 <- h.1 - qnorm(conf.level/2) * se.h.1
    LLh.1[LLh.1 < 0] <- 0

    #SECOND TRANSITION
    if(obj$class[4] == "Weibull"){
      kappa <- exp(obj$estimate[3])
      alpha <- exp(obj$estimate[4])
      log_alpha <- obj$estimate[4]
      S.2 <- exp(-kappa * (T2)^alpha) * exp(eta2)
      #Use delta method to compute baseline survival confidence intervals
      if(all(is.na(obj$Finv))){
        Var.loglogS.2 <- NA
      } else if (!is.null(x2new) & nP[2] > 0) {
        #matrix with as many columns as parameters in S2, and as many rows as times in T2
        #under weibull, log(-log(S(t))) = log_alpha + log_kappa + xtbeta + alpha*log(t)
        #so, each row of J is the gradient at a particular T2 wrt logalpha, logkappa, beta1, ..., betap
        J <- cbind(1, exp(log_alpha) * log(T2),
                   matrix(x2new, nrow = length(T2), ncol = length(x2new), byrow = T))
        #vector with as many elements as times in T2
        Var.loglogS.2 <- J %*% obj$Finv[c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])),
                                        c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2]))] %*% t(J)
      } else if (is.null(x2new) | nP[2] == 0) {
        J <- cbind(1, exp(log_alpha) * log(T2))
        Var.loglogS.2 <- J %*% obj$Finv[(1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0[1]):(nP0[1]+nP0[2])] %*% t(J)
      }

      h.2 <- alpha * kappa * (T2)^(alpha - 1) * exp(eta2)
      #hazard is exp(logalpha) * exp(logkappa) * t^(exp(logalpha)-1) * exp(xtbeta)
      #so, each row of J is gradient wrt logkappa, logalpha, beta1, ..., betap
      #WHY ISNT EVERYTHING MULTIPLIED BY exp(xtbeta) ? seems like it should be. I'm gonna add it
      if(all(is.na(obj$Finv))){
        Var.h.2 <- NA
      } else if (!is.null(x2new) & nP[2] > 0) {
        J <- cbind(h.2, h.2 * (1 + alpha * log(T2)), h.2 *
                     matrix(x2new, nrow = length(T2), ncol = length(x2new), byrow = T))
        Var.h.2 <- J %*% obj$Finv[c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])),
                                  c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2]))] %*% t(J)
      }
      else if (is.null(x2new) | nP[2] == 0) {
        J <- cbind(h.2, h.2 * (1 + alpha * log(T2)))
        Var.h.2 <- J %*% obj$Finv[(1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0[1]):(nP0[1]+nP0[2])] %*% t(J)
      }
    } else{ #for now the only alternative we entertain is piecewise constant
      stopifnot(obj$class[4] == "Piecewise Constant")
      basis2 <- get_basis(x = T2,knots = obj$knots_list[[2]],hazard = "piecewise")
      phi2 <- obj$estimate[(1+nP0[1]):(nP0[1]+nP0[2])]
      Lambda02 <- as.vector(basis2 %*% exp(phi2))
      S.2 <- exp(-Lambda02*exp(eta2))
      #Use delta method to compute baseline survival confidence intervals
      if(all(is.na(obj$Finv))){
        Var.loglogS.2 <- NA
      } else if (!is.null(x2new) & nP[2] > 0) {
        #matrix with as many columns as parameters in S1, and as many rows as times in T2
        #under piecewise, log(-log(S(t))) = log( basis1 %*% exp(phi1)) + log(xtbeta)
        #so, each row of J is the gradient at a particular T2 wrt phi1, ..., phiK, beta1, ..., betap
        J <- cbind(
          #first, this multiplies every row of basis1 by exp(phi1),
          #then, this divides every column of the resulting matrix by the vector Lambda01
          t(t(basis2) * exp(phi2)) / Lambda02,
          matrix(x2new, nrow = length(T2), ncol = length(x2new), byrow = T))
        #vector with as many elements as times in T2
        Var.loglogS.2 <- J %*% obj$Finv[c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])),
                                        c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2]))] %*%
          t(J)
      } else if (is.null(x2new) | nP[2] == 0) {
        J <- t(t(basis2) * exp(phi2)) / Lambda02
        Var.loglogS.2 <- J %*% obj$Finv[(1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0[1]):(nP0[1]+nP0[2])] %*% t(J)
      }

      #Use delta method compute baseline hazard confidence intervals
      #vector saying which interval each time falls into
      cut_cats2 <- rowSums(basis2!=0)
      if(T2[1]==0){cut_cats2[1] <- 1}
      stopifnot(length(cut_cats2)==length(T2))
      h.2 <- exp(phi2)[cut_cats2] * exp(eta2)
      #build a matrix with 0's everywhere, except on ith row, set column that T2i falls in to 1
      temp_mat <- matrix(data=0,nrow=length(T2),ncol=nP0[2])
      temp_mat[cbind(1:length(T2),cut_cats2)] <- 1
      temp_mat <- t(t(temp_mat) * exp(phi2))

      if(all(is.na(obj$Finv))){
        Var.h.2 <- NA
      } else if (!is.null(x2new) & nP[2] > 0) {
        J <- cbind(temp_mat * exp(eta2), h.2 *
                     matrix(x2new, nrow = length(T2), ncol = length(x2new), byrow = T))
        Var.h.2 <- J %*% obj$Finv[c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])),
                                  c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2]))] %*% t(J)
      }
      else if (is.null(x2new) | nP[2] == 0) {
        J <- temp_mat * exp(eta2)
        Var.h.2 <- J %*% obj$Finv[(1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0[1]):(nP0[1]+nP0[2])] %*% t(J)
      }
    }
    if(all(is.na(Var.loglogS.2))){
      se.loglogS.2 <- NA
    } else{
      se.loglogS.2 <- sqrt(diag(Var.loglogS.2))
      se.loglogS.2[is.nan(se.loglogS.2)] <- 0
    }
    LL.2 <- S.2^exp(-qnorm(conf.level/2) * se.loglogS.2)
    UL.2 <- S.2^exp(qnorm(conf.level/2) * se.loglogS.2)
    if(all(is.na(Var.h.2))){
      se.h.2 <- NA
    } else{
      se.h.2 <- sqrt(diag(Var.h.2))
      se.h.2[is.nan(se.h.2)] <- 0
    }
    LLh.2 <- h.2 + qnorm(conf.level/2) * se.h.2 #sign reversed because 0.025 quantile is negative
    ULh.2 <- h.2 - qnorm(conf.level/2) * se.h.2
    LLh.2[LLh.2 < 0] <- 0

    #THIRD TRANSITION
    if(obj$class[4] == "Weibull"){
      kappa <- exp(obj$estimate[5])
      alpha <- exp(obj$estimate[6])
      log_alpha <- obj$estimate[6]
      S.3 <- exp(-kappa * (T2)^alpha) * exp(eta3)
      #Use delta method to compute baseline survival confidence intervals
      if(all(is.na(obj$Finv))){
        Var.loglogS.3 <- NA
      } else if(!is.null(x3new) & nP[3] > 0) {
        #matrix with as many columns as parameters in S3, and as many rows as times in T2
        #under weibull, log(-log(S(t))) = log_alpha + log_kappa + xtbeta + alpha*log(t)
        #so, each row of J is the gradient at a particular T2 wrt logalpha, logkappa, beta1, ..., betap
        J <- cbind(1, exp(log_alpha) * log(T2),
                   matrix(x3new, nrow = length(T2), ncol = length(x3new), byrow = T))
        #vector with as many elements as times in T2
        Var.loglogS.3 <- J %*% obj$Finv[c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3])),
                                        c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3]))] %*% t(J)
      } else if (is.null(x3new) | nP[2] == 0) {
        J <- cbind(1, exp(log_alpha) * log(T2))
        Var.loglogS.3 <- J %*% obj$Finv[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])] %*% t(J)
      }

      h.3 <- alpha * kappa * (T2)^(alpha - 1) * exp(eta3)
      #hazard is exp(logalpha) * exp(logkappa) * t^(exp(logalpha)-1) * exp(xtbeta)
      #so, each row of J is gradient wrt logkappa, logalpha, beta1, ..., betap
      #kyu has code does not multiply the variances by exp(xtbeta), which it seems like it should be. I'm gonna add it
      if(all(is.na(obj$Finv))){
        Var.h.3 <- NA
      } else if (!is.null(x3new) & nP[3] > 0) {
        J <- cbind(h.3, h.3 * (1 + alpha * log(T2)), h.3 *
                     matrix(x3new, nrow = length(T2), ncol = length(x3new), byrow = T))
        Var.h.3 <- J %*% obj$Finv[c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3])),
                                  c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3]))] %*% t(J)
      }
      else if (is.null(x3new) | nP[2] == 0) {
        J <- cbind(h.3, h.3 * (1 + alpha * log(T2)))
        Var.h.3 <- J %*% obj$Finv[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])] %*% t(J)
      }
    } else{ #for now the only alternative we entertain is piecewise constant
      stopifnot(obj$class[4] == "Piecewise Constant")
      basis3 <- get_basis(x = T2,knots = obj$knots_list[[3]],hazard = "piecewise")
      phi3 <- obj$estimate[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])]
      Lambda03 <- as.vector(basis3 %*% exp(phi3))
      S.3 <- exp(-Lambda03*exp(eta3))
      #Use delta method to compute baseline survival confidence intervals
      if(all(is.na(obj$Finv))){
        Var.loglogS.3 <- NA
      } else if (!is.null(x3new) & nP[3] > 0) {
        #matrix with as many columns as parameters in S1, and as many rows as times in T2
        #under piecewise, log(-log(S(t))) = log( basis1 %*% exp(phi1)) + log(xtbeta)
        #so, each row of J is the gradient at a particular T2 wrt phi1, ..., phiK, beta1, ..., betap
        J <- cbind(
          #first, this multiplies every row of basis1 by exp(phi1),
          #then, this divides every column of the resulting matrix by the vector Lambda01
          t(t(basis3) * exp(phi3)) / Lambda03,
          matrix(x3new, nrow = length(T2), ncol = length(x3new), byrow = T))
        #vector with as many elements as times in T2
        Var.loglogS.3 <- J %*% obj$Finv[c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3])),
                                        c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3]))] %*%
          t(J)
      } else if (is.null(x3new) | nP[3] == 0) {
        J <- t(t(basis3) * exp(phi3)) / Lambda03
        Var.loglogS.3 <- J %*% obj$Finv[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])] %*% t(J)
      }

      #Use delta method compute baseline hazard confidence intervals
      #vector saying which interval each time falls into
      cut_cats3 <- rowSums(basis3!=0)
      if(T2[1]==0){cut_cats3[1] <- 1}
      stopifnot(length(cut_cats3)==length(T2))
      h.3 <- exp(phi3)[cut_cats3] * exp(eta3)
      #build a matrix with 0's everywhere, except on ith row, set column that T2i falls in to 1
      temp_mat <- matrix(data=0,nrow=length(T2),ncol=nP0[3])
      temp_mat[cbind(1:length(T2),cut_cats3)] <- 1
      temp_mat <- t(t(temp_mat) * exp(phi3))

      if(all(is.na(obj$Finv))){
        Var.h.3 <- NA
      } else if (!is.null(x3new) & nP[3] > 0) {
        J <- cbind(temp_mat * exp(eta3), h.3 *
                     matrix(x3new, nrow = length(T2), ncol = length(x3new), byrow = T))
        Var.h.3 <- J %*% obj$Finv[c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3])),
                                  c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3]))] %*% t(J)
      }
      else if (is.null(x3new) | nP[3] == 0) {
        J <- temp_mat * exp(eta3)
        Var.h.3 <- J %*% obj$Finv[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])] %*% t(J)
      }
    }
    if(all(is.na(Var.loglogS.3))){
      se.loglogS.3 <- NA
    } else{
      se.loglogS.3 <- sqrt(diag(Var.loglogS.3))
      se.loglogS.3[is.nan(se.loglogS.3)] <- 0
    }
    LL.3 <- S.3^exp(-qnorm(conf.level/2) * se.loglogS.3)
    UL.3 <- S.3^exp(qnorm(conf.level/2) * se.loglogS.3)
    if(all(is.na(Var.h.3))){
      se.h.3 <- NA
    } else{
      se.h.3 <- sqrt(diag(Var.h.3))
      se.h.3[is.nan(se.h.3)] <- 0
    }
    LLh.3 <- h.3 + qnorm(conf.level/2) * se.h.3 #sign reversed because 0.025 quantile is negative
    ULh.3 <- h.3 - qnorm(conf.level/2) * se.h.3
    LLh.3[LLh.3 < 0] <- 0

    T2h <- T2
    if (T2[1] == 0) {
      T2h <- T2h[-1]
      h.1 <- h.1[-1]
      LLh.1 <- LLh.1[-1]
      ULh.1 <- ULh.1[-1]
      h.2 <- h.2[-1]
      LLh.2 <- LLh.2[-1]
      ULh.2 <- ULh.2[-1]
      h.3 <- h.3[-1]
      LLh.3 <- LLh.3[-1]
      ULh.3 <- ULh.3[-1]
    }
    BH1_tbl <- data.frame(time = T2h, h.1 = h.1, LL.1 = LLh.1,
                          UL.1 = ULh.1)
    BH2_tbl <- data.frame(time = T2h, h.2 = h.2, LL.2 = LLh.2,
                          UL.2 = ULh.2)
    BH3_tbl <- data.frame(time = T2h, h.3 = h.3, LL.3 = LLh.3,
                          UL.3 = ULh.3)
    BS1_tbl <- data.frame(time = T2, S.1 = S.1, LL.1 = LL.1,
                          UL.1 = UL.1)
    BS2_tbl <- data.frame(time = T2, S.2 = S.2, LL.2 = LL.2,
                          UL.2 = UL.2)
    BS3_tbl <- data.frame(time = T2, S.3 = S.3, LL.3 = LL.3,
                          UL.3 = UL.3)
    value <- list(h.1 = BH1_tbl, h.2 = BH2_tbl, h.3 = BH3_tbl,
                  S.1 = BS1_tbl, S.2 = BS2_tbl, S.3 = BS3_tbl)
  }
  value$xnew <- xnew
  value$x1new <- x1new
  value$x2new <- x2new
  value$x3new <- x3new
  value$tseq <- tseq
  value$setup$model <- obj$setup$model
  value$class <- obj$class
  class(value) <- "pred.Freq_HReg2"
  return(value)
}

plot.pred.Freq_HReg2 <- function (x, plot.est = "Haz", xlab = NULL, ylab = NULL, ...)
{
  obj <- x
  T2seq <- x$tseq
  yLim <- NULL
  if (obj$class[2] == "ID") {
    if (is.null(ylab)) {
      if (plot.est == "Surv") {
        ylab <- "Survival"
      }
      if (plot.est == "Haz") {
        ylab <- "Hazard"
      }
    }
    if (is.null(xlab)) {
      xlab <- c("Time", "Time", "Time")
      if (obj$class[5] == "semi-Markov") {
        xlab[3] <- "Time since non-terminal event"
      }
    }
    if (is.null(yLim)) {
      if (plot.est == "Surv") {
        yLim <- seq(from = 0, to = 1, by = 0.2)
      }
      if (plot.est == "Haz") {
        ygrid <- (max(x$h.1$h.1,x$h.2$h.2,x$h.3$h.3,
                      x$h.1$UL.1, x$h.2$UL.2, x$h.3$UL.3,na.rm = TRUE) -
                    0)/5
        yLim <- seq(from = 0, to = max(x$h.1$h.1,x$h.2$h.2,x$h.3$h.3,
                                       x$h.1$UL.1, x$h.2$UL.2, x$h.3$UL.3,na.rm = TRUE), by = ygrid)
      }
    }
    if (plot.est == "Surv") {
      par(mfrow = c(1, 3))
      plot(range(T2seq), range(yLim), xlab = xlab[1], ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               S[1](t), "")), axes = FALSE)
      axis(1, at = T2seq)
      axis(2, at = yLim)
      lines(obj$S.1$time, obj$S.1$S.1, col = "blue", lwd = 3)
      lines(obj$S.1$time, obj$S.1$LL.1, col = "blue", lwd = 3,
            lty = 3)
      lines(obj$S.1$time, obj$S.1$UL.1, col = "blue", lwd = 3,
            lty = 3)
      plot(range(T2seq), range(yLim), xlab = xlab[2], ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               S[2](t), "")), axes = FALSE)
      axis(1, at = T2seq)
      axis(2, at = yLim)
      lines(obj$S.2$time, obj$S.2$S.2, col = "red", lwd = 3)
      lines(obj$S.2$time, obj$S.2$LL.2, col = "red", lwd = 3,
            lty = 3)
      lines(obj$S.2$time, obj$S.2$UL.2, col = "red", lwd = 3,
            lty = 3)
      plot(range(T2seq), range(yLim), xlab = xlab[3], ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               S[3](t), "")), axes = FALSE)
      axis(1, at = T2seq)
      axis(2, at = yLim)
      lines(obj$S.3$time, obj$S.3$S.3, col = "red", lwd = 3)
      lines(obj$S.3$time, obj$S.3$LL.3, col = "red", lwd = 3,
            lty = 3)
      lines(obj$S.3$time, obj$S.3$UL.3, col = "red", lwd = 3,
            lty = 3)
    }
    if (plot.est == "Haz") {
      par(mfrow = c(1, 3))
      plot(range(T2seq), range(yLim), xlab = xlab[1], ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               h[1](t), "")), axes = FALSE)
      axis(1, at = T2seq)
      axis(2, at = round(yLim, 4))
      lines(obj$h.1$time, obj$h.1$h.1, col = "blue", lwd = 3)
      lines(obj$h.1$time, obj$h.1$LL.1, col = "blue", lwd = 3,
            lty = 3)
      lines(obj$h.1$time, obj$h.1$UL.1, col = "blue", lwd = 3,
            lty = 3)
      plot(range(T2seq), range(yLim), xlab = xlab[2], ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               h[2](t), "")), axes = FALSE)
      axis(1, at = T2seq)
      axis(2, at = round(yLim, 4))
      lines(obj$h.2$time, obj$h.2$h.2, col = "red", lwd = 3)
      lines(obj$h.2$time, obj$h.2$LL.2, col = "red", lwd = 3,
            lty = 3)
      lines(obj$h.2$time, obj$h.2$UL.2, col = "red", lwd = 3,
            lty = 3)
      plot(range(T2seq), range(yLim), xlab = xlab[3], ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               h[3](t), "")), axes = FALSE)
      axis(1, at = T2seq)
      axis(2, at = round(yLim, 4))
      lines(obj$h.3$time, obj$h.3$h.3, col = "red", lwd = 3)
      lines(obj$h.3$time, obj$h.3$LL.3, col = "red", lwd = 3,
            lty = 3)
      lines(obj$h.3$time, obj$h.3$UL.3, col = "red", lwd = 3,
            lty = 3)
    }
  }
  invisible()
}
