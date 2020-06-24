FreqID_HReg_Rpath <- function(Formula, data, na.action="na.fail", subset=NULL,
                              hazard=c("weibull"),frailty=TRUE,model, knots_list = NULL,
                              penalty=c("scad","mcp","lasso"),
                              lambda_path=NULL, lambda_target=0, N_path_steps = 40, #passing in lambda_path overrides automatic calculation of path
                              a=NULL, mm_epsilon=1e-8,  select_tol=1e-4, fusion_tol=1e-3,
                              penalty_fusedcoef=c("none","fusedlasso"), lambda_fusedcoef_path=0,
                              penalty_fusedbaseline=c("none","fusedlasso"), lambda_fusedbaseline=0,
                              penweights_list=list(), mu_smooth_path=0, fit_method="prox_grad",
                              startVals=NULL, ball_R=Inf,
                              warm_start=TRUE, step_size_min=1e-6, step_size_max=1e6, step_size_init=1,
                              step_size_scale=0.5, #no checks implemented on these values!!
                              step_delta=0.5, maxit=300,
                              conv_crit = "omega", conv_tol=1e-6,
                              verbose=FALSE){

  # To start, I'm going to implement the PISTA algorithm of Wang et al (2014).  I know it is somewhat deficient compared to
  # the APISTA and PICASSO algorithms subsequently proposed by the same authors, but it's a good start with what I have previously implemented
  # I also think that relatively speaking, the cost of the nll and ngrad functions doesn't nicely decompose for coordinate methods in the same way? idk.
  # Gotta start somewhere, that's all.

  ##INITIALIZE OPTIONS##
  ##******************##
  #checks on the penalties and on 'a' happen in the underlying functions to minimize duplication

  na.action <- match.arg(na.action) #technically na.fail I think is the only one currently implemented
  if (na.action != "na.fail" & na.action != "na.omit") {
    stop("na.action should be either na.fail or na.omit")
  }

  # browser()

  ##DATA PREPROCESSING##
  ##******************##
  ##MAKE THIS MORE EFFICIENT BY GOING DIRECTLY INTO NUMERICS

  # check_pen_params(penalty=penalty,penalty_fusedcoef=penalty_fusedcoef,penalty_fusedbaseline=penalty_fusedbaseline,
  #                  lambda=lambda_target,lambda_fusedcoef=lambda_fusedcoef,lambda_fusedbaseline=lambda_fusedbaseline,
  #                  a=a)

  #rearrange input Formula object (which stores the different pieces of the input formula)
  #this seems to ensure it's in the correct form, doesn't seem to make a diff
  form2 <- Formula::as.Formula(paste(Formula[2], Formula[1], Formula[3],
                                     sep=""))
  #arrange data into correct form, extracting subset of variables
  #and observations needed in model
  data <- model.frame(form2, data=data, na.action=na.action,
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
  Xmat1 <- as.matrix(model.frame(formula(form2, lhs=0, rhs=1),
                                 data=data))
  Xmat2 <- as.matrix(model.frame(formula(form2, lhs=0, rhs=2),
                                 data=data))
  Xmat3 <- as.matrix(model.frame(formula(form2, lhs=0, rhs=3),
                                 data=data))
  n <- length(y1)

  ##PREPARE KNOTS AND BASIS FUNCTIONS FOR FLEXIBLE MODELS##
  ##*****************************************************##

  if(hazard %in% c("bspline","royston-parmar","piecewise")){
    if(is.null(knots_list)){
      p01 <- p02 <- p03 <- if(hazard=="bspline") 5 else 4
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


  ##Set up algorithmic parameters##
  ##*****************************##

  #by setting 'sparse_start==TRUE' we just find MLE solutions for the baseline and theta params,
  #and set the beta params to 0.
  if(is.null(startVals)){
    startVals <- get_start(y1=y1,y2=y2,delta1=delta1,delta2=delta2,Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,
                           basis1=basis1,basis2=basis2,basis3=basis3,basis3_y1=basis3_y1,
                           dbasis1=dbasis1,dbasis2=dbasis2,dbasis3=dbasis3,
                           hazard=hazard,frailty=frailty,model=model,sparse_start=TRUE)
  }

  nPtot <- length(startVals)
  nP1 <- ncol(Xmat1)
  nP2 <- ncol(Xmat2)
  nP3 <- ncol(Xmat3)
  nP0 <- nPtot - nP1 - nP2 - nP3
  nP_vec <- c(nP0,nP1,nP2,nP3)

  #if lambda_path is null then compute it according to Wang (2014), otherwise go with what was provided
  if(is.null(lambda_path)){
    startVals_grad <- ngrad_func(para=startVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                 hazard=hazard, frailty=frailty, model=model,
                                 basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                 dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)/n #notice the normalization by n

    # N_path_steps <- 40
    lambda0 <- max(abs(startVals_grad[-(1:(1+p01+p02+p03))])) #largest gradient of the betas
    #using Wang (2014) 3.3, solving for 'eta' given the other three pieces
    #add the -1 to account for starting at 0, so that final results are of length N_path_steps
    step_eta <- exp((log(lambda_target) - log(lambda0)) / (N_path_steps-1))
    if(lambda_target > lambda0 | step_eta < 0.9 | step_eta > 1){
      print(c(N_path_steps=N_path_steps,step_eta=step_eta,lambda_target=lambda_target,lambda0=lambda0))
      stop("provided lambda_target is too big, or provided N_path_steps is too small. require step_eta > 0.9")
    }
    #for now, we're just gonna let the same lambda govern all the betas, we'll work on the rest later.
    #note that we already 'have the solution' for the max value, because it is the startVal, FIX
    lambda_path <- as.matrix(lambda0 * step_eta^(0:(N_path_steps-1))) #
  } else{
    step_eta <- NULL
    lambda_path <- as.matrix(lambda_path)
  }

  colnames(lambda_path) <- paste0("lambda",1:ncol(lambda_path))

  ##ACTUALLY BEGIN THE PATH##
  ##***********************##


  solution_path_out <- solution_path_function(para=startVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                              Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                              hazard=hazard, frailty=frailty, model=model,
                                              basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                              dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                              penalty=penalty, lambda_path=lambda_path, a=a,
                                              penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef_path=lambda_fusedcoef_path,
                                              penalty_fusedbaseline="none", lambda_fusedbaseline=0, #idk what to do about the baseline fusing, for now I skip
                                              penweights_list=penweights_list, mu_smooth_path=mu_smooth_path, ball_R=ball_R,
                                              fit_method=fit_method, warm_start=warm_start,
                                              step_size_min=step_size_min, step_size_max=step_size_max, step_size_init=step_size_init,
                                              step_size_scale=step_size_scale, #no checks implemented on these values!!
                                              step_delta=step_delta,
                                              maxit=maxit,
                                              conv_crit = conv_crit,
                                              conv_tol=conv_tol,
                                              verbose=verbose)

  # if(ncol(lambda_path) == 1){ #as long as all the fused lasso lambdas are the same, make the plots
  #   plot_out_ests <- reshape2::melt(as.data.frame(solution_path_out$out_ests), id="lambda1") %>%
  #     dplyr::mutate(vartype=ifelse(stringr::str_detect(variable,"_1"),"Beta1",
  #                                  ifelse(stringr::str_detect(variable,"_2"),"Beta2",
  #                                         ifelse(stringr::str_detect(variable,"_3"),"Beta3",
  #                                                "Baseline Params"))))
  #
  #
  #   plot_out <- ggplot2::qplot(lambda1, value, group=variable, geom="line", data=plot_out_ests) +
  #     ggplot2::facet_wrap(~ vartype,scales="free_y") +
  #     # geom_vline(xintercept=lambda[bic_best_index,]) +
  #     ggplot2::theme_classic() +
  #     ggplot2::theme(legend.position="none") +
  #     ggplot2::ggtitle(paste("Regularization Path for",penalty))
  #   # print(plot_out)
  # }
  #
  #

  solution_path_out$step_eta <- step_eta

  return(solution_path_out)

}


