####WRAPPER FUNCTION FOR REGULARIZED ESTIMATION####
##***********************************************##

#wrapper function for fitting a model for a single choice of lambdas, and taking the best of a number of different methods
FreqID_HReg_R <- function(Formula, data, na.action="na.fail", subset=NULL,
                          hazard=c("weibull"),frailty=TRUE,model, knots_list = NULL,
                          penalty=c("scad","mcp","lasso"), lambda, a=NULL, mm_epsilon=1e-8,
                          penalty_fusedcoef=c("none","fusedlasso"), lambda_fusedcoef=0,
                          penalty_fusedbaseline=c("none","fusedlasso"), lambda_fusedbaseline=0,
                          penweights_list=list(), verbose=FALSE, startVals=NULL,
                          # randomize_start=FALSE,randomize_seed=NULL,
                          est_var=FALSE,optimization_method="prox_grad",select_tol=1e-04,fusion_tol=1e-3,
                          maxit=300, step_size_init=1, step_size_min=1e-6, step_size_max=1e6,step_size_scale=0.5,
                          conv_crit="nll_pen_change", conv_tol=1e-7, num_restarts=2){

  ##INITIALIZE OPTIONS##
  ##******************##
  #checks on the penalties and on 'a' happen in the underlying functions to minimize duplication

  na.action <- match.arg(na.action) #technically na.fail I think is the only one currently implemented
  if (na.action != "na.fail" & na.action != "na.omit") {
    stop("na.action should be either na.fail or na.omit")
  }

  #construct modeling control list
  #randomize_start says whether to multiply starting values by a random Unif(0.8,1.2) value
  #est_var says whether to invert and report final hessian matrix as covariance matrix
  #optim_method says whether to use the newton raphson algorithm written, or to feed model into nleqslv
  #[method]_mle indicates to fit the unregularized MLE using provided start values, before beginning regularization
  #select_tol is the threshold for setting a small estimate to 0 for variable selection


  if(!all(tolower(optimization_method) %in% c('prox_grad','prox_grad_mle','prox_grad_nmaccel','prox_grad_nmaccel_mle',
                                              'newton_mle','newton'))){
    warning("unknown optimization method. Options are 'prox_grad','prox_grad_mle','prox_grad_nmaccel','prox_grad_nmaccel_mle',
                                              'newton_mle','newton'")
  }

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

  if(hazard %in% c("bspline","royston-parmar","piecewise")){
    if(is.null(knots_list)){
      p01 <- p02 <- p03 <- if(hazard=="bspline") 5 else 3
      knots_list <- get_default_knots_list(y1,y2,delta1,delta2,p01,p02,p03,hazard,model)
    }
    basis1 <- get_basis(x = y1,knots = knots_list[[1]],hazard = hazard)
    basis2 <- get_basis(x = y1,knots = knots_list[[2]],hazard = hazard)
    dbasis1 <- get_basis(x = y1,knots = knots_list[[1]],hazard = hazard,deriv = TRUE)
    dbasis2 <- get_basis(x = y1,knots = knots_list[[2]],hazard = hazard,deriv = TRUE)
    if(model=="semi-markov"){
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

  ##OPTIMIZATION ROUTINES##
  ##*********************##
  # rand_seed <- NA
  # if(!is.null(con[["randomize_start"]]) && con[["randomize_start"]]){
  #   if(is.null(con[["randomize_seed"]])){
  #     rand_seed <- sample(x=10000,size=1)
  #   } else{
  #     rand_seed <- con[["randomize_seed"]]
  #   }
  #   set.seed(rand_seed)
  #   startVals <- startVals*runif(n=length(startVals),min=0.8,max=1.2)
  # }

  ##GET INITIAL MLE##
  ##**************##
  if(any(tolower(optimization_method) %in% c("prox_grad_mle","prox_grad_nmaccel_mle","newton_mle"))){
    # parahat_mle <- proximal_gradient_descent_nmaccel(para=startVals,
    #                                                  y1=y1, y2=y2, delta1=delta1, delta2=delta2,
    #                                                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
    #                                                  hazard=hazard,frailty=frailty,model=model,
    #                                                  basis1 = basis1, basis2 = basis2, basis3 = basis3, basis3_y1 = basis3_y1,
    #                                                  dbasis1 = dbasis1, dbasis2 = dbasis2, dbasis3 = basis3,
    #                                                  penalty="lasso",lambda=0, a=1,
    #                                                  penalty_fusedcoef="none", lambda_fusedcoef=0,
    #                                                  penalty_fusedbaseline="none", lambda_fusedbaseline=0,
    #                                                  penweights_list = list(), mu_smooth_fused = NULL,
    #                                                  control=con_prox, verbose=FALSE)$estimate

    parahat_mle <- stats::optim(par = startVals, fn = nll_func, gr = ngrad_func,
                         y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                         Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                         hazard=hazard,frailty=frailty,model=model,
                         basis1 = basis1, basis2 = basis2, basis3 = basis3, basis3_y1 = basis3_y1,
                         dbasis1 = dbasis1, dbasis2 = dbasis2, dbasis3 = basis3,
                         control=list(maxit=300),method = if(hazard=="royston-parmar") "BFGS" else "L-BFGS")$par

  } else{
    parahat_mle <- NULL
  }

  ##And now, the multiple starting point game begins!
  ##I allow my newton algorithm, or nleqslv algorithm,
  ##each of which can start at the provided/default start vals or the MLE of the ID,
  ##and then I choose and report the winner
  final_nll_pen <- Inf

  if(any(tolower(optimization_method) %in% c("prox_grad"))){
    if(verbose){print(paste("prox algorithm starting:"))}
    optim_out <- proximal_gradient_descent(para=startVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                           Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                           hazard=hazard, frailty=frailty, model=model,
                                           basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                           dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                           penalty=penalty, lambda=lambda, a=a,
                                           penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                           penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                           penweights_list=penweights_list, mu_smooth_fused=0,
                                           step_size_init=step_size_init,step_size_min=step_size_min,step_size_max=step_size_max,step_size_scale=step_size_scale,
                                           maxit=maxit, conv_crit=conv_crit, conv_tol=conv_tol,
                                           ball_R=Inf, verbose=verbose)
    if(verbose){print(optim_out$final_nll_pen)}
    if(optim_out$final_nll_pen < final_nll_pen){
      fit_method <- "prox_grad"
      if(verbose){print("prox_grad is best yet!")}
      finalVals <- optim_out$estimate
      final_nll <- optim_out$final_nll
      final_nll_pen <- optim_out$final_nll_pen
      iter_vec <- optim_out$niter
      fit_code <- optim_out$fit_code
      trace_mat <- optim_out$trace_mat
      final_ngrad <- optim_out$final_ngrad
      final_nhess <- NULL
    }
  }

  if(any(tolower(optimization_method) %in% c("prox_grad_mle"))){
    if(verbose){print(paste("prox algorithm starting:"))}
    optim_out <- proximal_gradient_descent(para=parahat_mle, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                           Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                           hazard=hazard, frailty=frailty, model=model,
                                           basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                           dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                           penalty=penalty, lambda=lambda, a=a,
                                           penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                           penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                           penweights_list=penweights_list, mu_smooth_fused=0,
                                           step_size_init=step_size_init,step_size_min=step_size_min,step_size_max=step_size_max,step_size_scale=step_size_scale,
                                           maxit=maxit, conv_crit=conv_crit, conv_tol=conv_tol,
                                           ball_R=Inf, verbose=verbose)
    if(verbose){print(optim_out$final_nll_pen)}
    if(optim_out$final_nll_pen < final_nll_pen){
      fit_method <- "prox_grad_mle"
      if(verbose){print("prox_grad_mle is best yet!")}
      finalVals <- optim_out$estimate
      final_nll <- optim_out$final_nll
      final_nll_pen <- optim_out$final_nll_pen
      iter_vec <- optim_out$niter
      fit_code <- optim_out$fit_code
      trace_mat <- optim_out$trace_mat
      final_ngrad <- optim_out$final_ngrad
      final_nhess <- NULL
    }
  }


  if(any(tolower(optimization_method) %in% c("prox_grad_nmaccel"))){
    if(verbose){print(paste("prox_grad_nmaccel algorithm starting:"))}
    optim_out <- proximal_gradient_descent_nmaccel(para=startVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                                   Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                                   hazard=hazard, frailty=frailty, model=model,
                                                   basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                                   dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                                   penalty=penalty, lambda=lambda, a=a,
                                                   penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                                   penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                                   penweights_list=penweights_list, mu_smooth_fused=0, #for now, only consider inexact prox method
                                                   step_size_init_x=step_size_init,step_size_init_y=step_size_init,
                                                   step_size_min=step_size_min,step_size_max=step_size_max,step_size_scale=step_size_scale,
                                                   maxit=maxit, conv_crit=conv_crit, conv_tol=conv_tol,
                                                   ball_R=Inf, verbose=verbose)
    if(verbose){print(optim_out$final_nll_pen)}
    if(optim_out$final_nll_pen < final_nll_pen){
      fit_method <- "prox_grad_nmaccel"
      if(verbose){print("prox_grad_nmaccel is best yet!")}
      finalVals <- optim_out$estimate
      final_nll <- optim_out$final_nll
      final_nll_pen <- optim_out$final_nll_pen
      iter_vec <- optim_out$niter
      fit_code <- optim_out$fit_code
      trace_mat <- optim_out$trace_mat
      final_ngrad <- optim_out$final_ngrad
      final_nhess <- NULL
    }
  }

  if(any(tolower(optimization_method) %in% c("prox_grad_nmaccel_mle"))){
    if(verbose){print(paste("prox_grad_nmaccel_mle algorithm starting:"))}
    optim_out <- proximal_gradient_descent_nmaccel(para=parahat_mle, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                                   Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                                   hazard=hazard, frailty=frailty, model=model,
                                                   basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                                   dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                                   penalty=penalty, lambda=lambda, a=a,
                                                   penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                                   penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                                   penweights_list=penweights_list, mu_smooth_fused=0, #for now, only consider inexact prox method
                                                   step_size_init_x=step_size_init,step_size_init_y=step_size_init,
                                                   step_size_min=step_size_min,step_size_max=step_size_max,step_size_scale=step_size_scale,
                                                   maxit=maxit, conv_crit=conv_crit, conv_tol=conv_tol,
                                                   ball_R=Inf, verbose=verbose)
    if(verbose){print(optim_out$final_nll_pen)}
    if(optim_out$final_nll_pen < final_nll_pen){
      fit_method <- "prox_grad_nmaccel_mle"
      if(verbose){print("prox_grad_nmaccel_mle is best yet!")}
      finalVals <- optim_out$estimate
      final_nll <- optim_out$final_nll
      final_nll_pen <- optim_out$final_nll_pen
      iter_vec <- optim_out$niter
      fit_code <- optim_out$fit_code
      trace_mat <- optim_out$trace_mat
      final_ngrad <- optim_out$final_ngrad
      final_nhess <- NULL
    }
  }



  if(any(tolower(optimization_method) %in% c("newton"))){
    if(verbose){print("newton")}
    optim_out <- newton_raphson_mm(para=startVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                   Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                   hazard=hazard, frailty=frailty, model=model,
                                   penalty=penalty, lambda=lambda, a=a,
                                   penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                   penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                   penweights_list=penweights_list, mm_epsilon=mm_epsilon, select_tol=select_tol,
                                   verbose=verbose, maxit=maxit,step_size_min=step_size_min,
                                   conv_crit=conv_crit,conv_tol=conv_tol,num_restarts=num_restarts)
    if(verbose){print(optim_out$final_nll_pen)}
    if(optim_out$final_nll_pen < final_nll_pen){
      fit_method <- "newton"
      if(verbose){print("newton is best yet!")}
      finalVals <- optim_out$estimate
      final_nll <- optim_out$final_nll
      final_nll_pen <- optim_out$final_nll_pen
      iter_vec <- optim_out$niter
      fit_code <- optim_out$fit_code
      trace_mat <- optim_out$trace_mat
      final_ngrad <- optim_out$final_ngrad
      final_nhess <- optim_out$final_nhess
    }
  }

  if(any(tolower(optimization_method) %in% c("newton_mle"))){
    if(verbose){print("newton_mle")}
    optim_out <- newton_raphson_mm(para=parahat_mle, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                   Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                   hazard=hazard, frailty=frailty, model=model,
                                   penalty=penalty, lambda=lambda, a=a,
                                   penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                   penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                   penweights_list=penweights_list, mm_epsilon=mm_epsilon, select_tol=select_tol,
                                   verbose=verbose, maxit=maxit,step_size_min=step_size_min,
                                   conv_crit=conv_crit,conv_tol=conv_tol,num_restarts=num_restarts)
    if(verbose){print(optim_out$final_nll_pen)}
    if(optim_out$final_nll_pen < final_nll_pen){
      if(verbose){print("newton_mle is best yet!")}
      fit_method <- "newton_mle"
      finalVals <- optim_out$estimate
      final_nll <- optim_out$final_nll
      final_nll_pen <- optim_out$final_nll_pen
      iter_vec <- optim_out$niter
      fit_code <- optim_out$fit_code
      trace_mat <- optim_out$trace_mat
      final_ngrad <- optim_out$final_ngrad
      final_nhess <- optim_out$final_nhess
    }
  }

  names(final_ngrad) <- names(finalVals) <- names(startVals)

  #break out the beta vectors from the larger parameter vector, using nP0 to correctly pad out the baseline hazard and theta parameters
  # beta1 <- para[(1+nP0):(nP0+nP1)]
  # beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
  # beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]

  nPtot <- length(startVals)
  nP1 <- ncol(Xmat1)
  nP2 <- ncol(Xmat2)
  nP3 <- ncol(Xmat3)
  nP0 <- nPtot - nP1 - nP2 - nP3
  nP_vec <- c(nP0=nP0,nP1=nP1,nP2=nP2,nP3=nP3)
  n <- length(y1)


  beta1_selected <- if(nP1 != 0) finalVals[(1+nP0):(nP0+nP1)] else numeric(0)
  nP1_selected <- sum(beta1_selected != 0)

  beta2_selected <- if(nP2 != 0) finalVals[(1+nP0+nP1):(nP0+nP1+nP2)] else numeric(0)
  nP2_selected <- sum(beta2_selected != 0)

  beta3_selected <- if(nP3 != 0) finalVals[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)] else numeric(0)
  nP3_selected <- sum(beta3_selected != 0)

  nPtot_selected <- nP0 + nP1_selected + nP2_selected + nP3_selected
  nP_selected_vec <- c(nP0=nP0,nP1=nP1_selected,nP2=nP2_selected,nP3=nP3_selected)

  tempaic <- 2*final_nll + 2* nPtot_selected
  tempbic <- 2*final_nll + log(n) * nPtot_selected
  tempgcv <- final_nll/(n^2*(1-nPtot_selected/n)^2)

  # The provisional idea for degrees of freedom is to take the number of unique nonzero estimates, e.g.,
  # nPtot_selected_unique <- sum(unique(finalVals_selected) != 0)
  #but here we use a tolerance because we can't get exact equality
  #NOTE THIS DOES NOT INCORPORATE POTENTIAL EQUALITY OF BASELINE PARAMETERS
  #ALSO IT ASSUMES THE SAME COVARIATES IN THE SAME ORDER ACROSS ALL THREE HAZARDSSSS!!!
  if(nP1==nP2 & nP2==nP3){
    nPtot_unique <- nP0 + sum(beta1_selected != 0 &
                                abs(beta1_selected - beta2_selected) > fusion_tol &
                                abs(beta1_selected - beta3_selected) > fusion_tol) +
      sum(beta2_selected != 0 &
            abs(beta2_selected - beta3_selected) > fusion_tol) +
      sum(beta3_selected != 0)

    #currently, these account for fusion in the 'degrees of freedom' by counting unique parameters
    tempaic_unique <- 2*final_nll + 2* nPtot_unique
    tempbic_unique <- 2*final_nll + log(n) * nPtot_unique
    tempgcv_unique <- final_nll/(n^2*(1-nPtot_unique/n)^2)
  } else{
    nPtot_unique <- tempaic_unique <- tempbic_unique <- tempgcv_unique <- NA
  }

  #alternative definition of degrees of freedom from Sennhenn-Reulen & Kneib via Gray (1993)
  # browser()
  if(!is.null(final_nhess)){
    final_nhess_nopen <- nhess_func(para=finalVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                    hazard=hazard, frailty=frailty, model=model)
    final_cov <- tryCatch(solve(final_nhess),
                          error=function(cnd){
                            message(cnd)
                            cat("\n")
                            return(NULL)
                          })
    if(!is.null(final_cov)){
      df_est <- sum(diag( final_nhess_nopen %*% final_cov ))
    } else{
      df_est <- NA
    }

    tempaic_estdf <- 2*final_nll + 2* df_est
    tempbic_estdf <- 2*final_nll + log(n) * df_est
    tempgcv_estdf <- final_nll/(n^2*(1-df_est/n)^2)

  } else{
    df_est <- tempaic_estdf <- tempbic_estdf <- tempgcv_estdf <- NA
  }


  value <- list(estimate=finalVals,# estimate_selected=finalVals_selected,
                nll=final_nll, nll_pen=final_nll_pen,
                nP_vec=nP_vec,
                nP_selected_vec=nP_selected_vec,
                nPtot=nPtot, nPtot_selected=nPtot_selected,
                nPtot_unique=nPtot_unique, df_est=df_est,
                iter_vec=iter_vec,fit_code=fit_code,
                AIC=tempaic, BIC=tempbic, GCV=tempgcv,
                AIC_unique=tempaic_unique,BIC_unique=tempbic_unique,GCV_unique=tempgcv_unique,
                AIC_estdf=tempaic_estdf,BIC_estdf=tempbic_estdf,GCV_estdf=tempgcv_estdf,
                hazard=hazard, frailty=frailty, model=model,
                penalty=penalty, lambda=lambda, a=a, parahat=parahat_mle, mm_epsilon=mm_epsilon, select_tol=select_tol, fusion_tol=fusion_tol,
                penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                startVals=startVals, fit_method=fit_method, knots_list=knots_list,
                ngrad=final_ngrad)

  # if(con[["est_var"]]==TRUE){
  #   #compute naive variance
  #   final_ngrad_mat <- ngrad_mat_func(para=finalVals,y1=y1,y2=y2,delta1=delta1,delta2=delta2,
  #                             Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,frailty=frailty,hazard=hazard,model=model)
  #   #compute cheese (there is maybe a small discrepancy with Hunter & Li about how many times to divide by n, but this makes sense to me)
  #   cheese <- t(final_ngrad_mat) %*% final_ngrad_mat -
  #     (final_ngrad %*% t(final_ngrad)) / n
  #   #compute sandwich estimate
  #   sandwichvar <- final_cov %*% cheese %*% final_cov
  #   value[["naivevar"]]=final_cov
  #   value[["sandwichvar"]]=sandwichvar
  # }

  return(value)
  invisible()

}
