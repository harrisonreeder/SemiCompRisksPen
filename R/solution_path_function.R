##**********************************##
####MANUAL NEWTON RAPHSON FUNCTION####
##**********************************##

#para is vector of starting values
#nP1-3 are integers for how many covariates are included in each arm (e.g., length of corresponding beta)
#penalty is character string indicating which type of penalty to apply
#y1-2 are first and second event times
#delta1-2 are indicators of observation of first and second events
#Xmat1-3 are design matrices corresponding to each transition
#lambda is either a scalar value for the penalty parameter shared by all three transitions,
#or a three-vector with the values for each of the three transitions
#a is a scalar value for the second penalty parameter in SCAD, MCP, and adaptive lasso
#Default is 3.7 for SCAD, 3 for MCP, and 1 for adaptive lasso
#frailty indicates whether to include a gamma frailty or not, currently everything is designed assuming a frailty
#mm_epsilon is the perturbation included in the local quadratic approximation of the penalty term
#verbose is a boolean indicating whether to print the running fitting details
#control is a list with options, outlined below
#parahat is a vector of weights for adaptive lasso, typically arising as the MLE fit of the full model
newton_raphson_mm <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                              hazard, frailty, model,
                              penalty, lambda, a,
                              penalty_fusedcoef, lambda_fusedcoef,
                              penalty_fusedbaseline, lambda_fusedbaseline,
                              penweights_list=penweights_list, mm_epsilon, verbose,
                              maxit=300,step_size_min=0.00001,
                              conv_crit="nll_pen_change",conv_tol=1e-05,num_restarts=2,
                              select_tol=1e-4){

  ##Set up control pieces##
  ##*********************##
  # browser()
  #construct control list
  #maxit is maximum number of newton-raphson iterations allowed
  #min_lr is the minimum learning rate allowed through step-halving process.
  #conv_crit determines criterion used to judge convergence
  #'est_change_norm' looks at l1 norm of change in parameter values (scaled by l1 norm of prior values): sum(abs(finalVals-prevVals))/sum(abs(prevVals))
  #'max_est_change' looks at largest absolute change in a parameter value
  #'ll_maj_change' looks at change in majorized log-likelihood (expanded around previous location)
  #'ll_reg_change' looks at change in regularized log-likelihood
  #'maj_grad_norm' looks at l1 norm of majorized gradient
  # if (length(noNms <- namc[!namc %in% nmsC])){
  #   warning("unknown names in control: ", paste(noNms, collapse=", "))
  # }
  if(!all(conv_crit %in% c("est_change_norm","max_est_change","nll_maj_change","nll_pen_change","maj_grad_norm"))){
    stop("unknown convergence criterion.")
  }
  if(maxit < 0){
    stop("maxit must be nonnegative.")
  }
  if(num_restarts <0){
    stop("number of algorithm restarts must be nonnegative.")
  }

  n <- length(y1)
  Xmat1 <- if(!is.null(Xmat1)) as.matrix(Xmat1) else matrix(nrow=n,ncol=0)
  Xmat2 <- if(!is.null(Xmat2)) as.matrix(Xmat2) else matrix(nrow=n,ncol=0)
  Xmat3 <- if(!is.null(Xmat3)) as.matrix(Xmat3) else matrix(nrow=n,ncol=0)

  nPtot <- length(para)
  nP1 <- ncol(Xmat1)
  nP2 <- ncol(Xmat2)
  nP3 <- ncol(Xmat3)
  nP0 <- nPtot - nP1 - nP2 - nP3

  #taking in any weights that have already been input,
  pen_mat_list_temp <- contrast_mat_list(nP0=nP0,nP1 = nP1,nP2 = nP2,nP3 = nP3,
                                         penalty_fusedcoef = penalty_fusedcoef, lambda_fusedcoef = lambda_fusedcoef,
                                         penalty_fusedbaseline = penalty_fusedbaseline,lambda_fusedbaseline = lambda_fusedbaseline,
                                         hazard = hazard, penweights_list = penweights_list)

  #append them into single weighted matrix
  #if there is no fusion, this should return NULL
  pen_mat_w <- do.call(what = rbind,args = pen_mat_list_temp[["pen_mat_list"]])

  #vector with length equal to the total number of rows of pen_mat_w, with each entry lambda corresponding to that contrast
  #if there is no fusion, this should return NULL
  lambda_f_vec <- pen_mat_list_temp[["lambda_f_vec"]]

  #now, the matrix with the penalty baked in for use in smooth approximation of the fused penalty
  if(is.null(pen_mat_w)){
    pen_mat_w_lambda <- NULL
  } else{
    pen_mat_w_lambda <- diag(lambda_f_vec) %*% pen_mat_w #consider shifting to Matrix package version, because it never goes to cpp
  }

  #THERE HAS GOT TO BE A BETTER WAY TO INTEGRATE THIS INTO THE MM ALGORITHM, BECAUSE AT PRESENT THE FUNCTIONS REALLY DUPLICATE A LOT OF EFFORT, BUT IT'S SOMETHING!


  ##Run algorithm##
  ##*************##

  nll_pen_trace <- rep(NA,maxit)
  trace_mat <- prevVals <- finalVals <- startVals <- para
  fit_code <- 4 #follows convention from 'nleqslv' L-BFGS that is 1 if converged, 4 if maxit reached, 3 if stalled
  bad_step_count <- restart_count <- 0
  i <- 1
  while(i <= maxit){
    if(verbose)print(i)
    lr <- 1

    Ek <- pen_maj_mat_func(para=finalVals,nP1=nP1,nP2=nP2,nP3=nP3,penalty=penalty,lambda=lambda, a=a,
                           penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                           penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                           penweights_list=penweights_list, mm_epsilon=mm_epsilon, hazard=hazard)

    curr_nll <- nll_func(para=finalVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                         Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                         hazard=hazard, frailty=frailty, model=model,
                         basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                         dbasis1=NULL, dbasis2=NULL, dbasis3=NULL)

    curr_pen <- pen_func(para=finalVals,nP1=nP1,nP2=nP2,nP3=nP3,penalty=penalty,lambda=lambda, a=a,
                         penweights_list=penweights_list, pen_mat_w_lambda = pen_mat_w_lambda)

    curr_nll_pen <- curr_nll + (n * curr_pen)

    ##Compute next newton step##
    ##************************##

    temp_ngrad <- ngrad_func(para=finalVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                             Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                             hazard=hazard, frailty=frailty, model=model,
                             basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                             dbasis1=NULL, dbasis2=NULL, dbasis3=NULL)

    temp_nhess <- nhess_func(para=finalVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                             Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                             hazard=hazard, frailty=frailty, model=model)

    grad_step <- tryCatch(solve(temp_nhess + n*Ek, temp_ngrad + n*Ek %*% finalVals),
                          error=function(cnd){
                            message(cnd)
                            cat("\n")
                            return(NULL)
                          }) #chol2inv(chol(temp_nhess_deriv + n*Ek))


    ##Check if step is good, otherwise either restart or exit algorithm##
    ##*****************************************************************##

    if(is.null(grad_step) & restart_count < num_restarts){
      print("Problem trying to invert Hessian matrix, restarting at perturbed start values.")
      trace_mat <- prevVals <- finalVals <- startVals <- runif(n=1,min=0.8,max=1.2)*para
      fit_code <- 4 #follows convention from 'nleqslv' that is 1 if converged, 4 if maxit reached, 3 if stalled
      bad_step_count <- 0
      restart_count <- restart_count + 1
      i <- 1
      next
    } else if(is.null(grad_step) & restart_count == num_restarts){
      warning("Problem trying to invert Hessian matrix, exiting with fit_code 5, likely not converged.")
      finalVals <- as.numeric(finalVals)
      names(finalVals) <- names(para)
      nll_pen_trace[i] <- next_nll_pen

      #replace any of the betas that are under the selection tolerance with 0, ignoring the baseline variables
      if(nP1+nP2+nP3 > 0){
        if(any(abs(finalVals[(1+nP0):nPtot]) < select_tol)) {
          finalVals[(1+nP0):nPtot][ abs(finalVals[(1+nP0):nPtot]) < select_tol] <- 0
        }
      }


      return(list(estimate=finalVals, niter=i,
                  final_ngrad_pen=as.numeric(temp_ngrad + n*Ek %*% finalVals),
                  final_nll=curr_nll,
                  final_nll_pen=curr_nll_pen,
                  nll_pen_trace=nll_pen_trace,
                  fit_code=5,
                  control=list(step_size_init=1,
                               step_size_min=step_size_min,
                               conv_crit=conv_crit,conv_tol=conv_tol,
                               maxit=maxit,mm_epsilon=mm_epsilon,
                               num_restarts=num_restarts),
                  # trace_mat=as.matrix(trace_mat),
                  # control=con,
                  final_step_size=lr,
                  restart_count=restart_count))

    }

    next_nll <- nll_func(para=finalVals-lr*grad_step, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                         Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                         frailty=frailty, hazard=hazard, model=model,
                         basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                         dbasis1=NULL, dbasis2=NULL, dbasis3=NULL)

    next_pen <- pen_func(para=finalVals-lr*grad_step,nP1=nP1,nP2=nP2,nP3=nP3,penalty=penalty,lambda=lambda, a=a,
                         penweights_list=penweights_list, pen_mat_w_lambda = pen_mat_w_lambda)

    next_maj <- maj_func(para=finalVals-lr*grad_step, para0=finalVals,
                         nP1=nP1,nP2=nP2,nP3=nP3,penalty=penalty,lambda=lambda, a=a,
                         penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                         penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline, pen_mat_w_lambda = pen_mat_w_lambda,
                         penweights_list=penweights_list, mm_epsilon=mm_epsilon, hazard=hazard)

    next_nll_pen <- next_nll + (n * next_pen)
    next_nll_maj <- next_nll + (n * next_maj)

    if(is.nan(next_nll_maj) & restart_count < num_restarts){
      print("Candidate iteration step yields infinite function value, restarting at perturbed start values.")
      trace_mat <- prevVals <- finalVals <- startVals <- runif(n=1,min=0.8,max=1.2)*para
      fit_code <- 4 #follows convention from 'nleqslv' that is 1 if converged, 4 if maxit reached, 3 if stalled
      bad_step_count <- 0
      restart_count <- restart_count + 1
      i <- 1
      next
    } else if(is.nan(next_nll_maj) & restart_count == num_restarts){
      warning("Candidate iteration step yields infinite function value. Exiting with fit_code 6, likely not converged.")

      finalVals <- as.numeric(finalVals)
      names(finalVals) <- names(para)
      nll_pen_trace[i] <- next_nll_pen

      #replace any of the betas that are under the selection tolerance with 0, ignoring the baseline variables
      if(nP1+nP2+nP3 > 0){
        if(any(abs(finalVals[(1+nP0):nPtot]) < select_tol)) {
          finalVals[(1+nP0):nPtot][ abs(finalVals[(1+nP0):nPtot]) < select_tol] <- 0
        }
      }


      return(list(estimate=finalVals, niter=i,
                  final_ngrad_pen=as.numeric(temp_ngrad + n*Ek %*% finalVals),
                  final_nll=curr_nll,
                  final_nll_pen=curr_nll_pen,
                  nll_pen_trace=nll_pen_trace,
                  fit_code=6,
                  control=list(step_size_init=1,
                               step_size_min=step_size_min,
                               conv_crit=conv_crit,conv_tol=conv_tol,
                               maxit=maxit,mm_epsilon=mm_epsilon,
                               num_restarts=num_restarts),
                  # trace_mat=as.matrix(trace_mat),
                  # control=con,
                  final_step_size=lr,
                  restart_count=restart_count))


    }

    #this step looks at majorized expansion around current point, and checks that the next step gets to a majorized value better than the current.
    #because current regularized value is also the point of majorization and majorized function is tangent to regularized function at current point
    #arguably, this could look just at regularized nll rather than majorized...because majorized has a tendency to be so pointy in some dimensions and leaves little wiggle room.
    while(lr > step_size_min & next_nll_maj > curr_nll_pen){
      if(verbose)print(lr)
      lr <- lr/2
      # browser()

      next_nll_pen <- next_nll_maj <- nll_maj_func(para=finalVals-lr*grad_step, para0=finalVals,
                                                   y1=y1, delta1=delta1, y2=y2, delta2=delta2,
                                                   Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                                   hazard=hazard, frailty=frailty, model=model,
                                                   penalty=penalty, lambda=lambda, a=a,
                                                   penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                                   penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline, pen_mat_w_lambda = pen_mat_w_lambda,
                                                   penweights_list=penweights_list, mm_epsilon=mm_epsilon)
    }

    if(lr < step_size_min){ #keep track of how many times the learning rate drops below its minimum but we still step
      bad_step_count <- bad_step_count+1
    }

    if(bad_step_count == 3 & restart_count < num_restarts){ #restart if a few bad steps in a row
      print("Algorithm cannot find good next step. Restarting at perturbed start values.")
      trace_mat <- prevVals <- finalVals <- startVals <- runif(n=1,min=0.8,max=1.2)*para
      fit_code <- 4 #follows convention from 'nleqslv' that is 1 if converged, 4 if maxit reached, 3 if stalled
      bad_step_count <- 0
      restart_count <- restart_count + 1
      i <- 1
      next
    } else if(bad_step_count == 3 & restart_count == num_restarts){
      warning("Algorithm cannot find good next step. Exiting with fit_code 7, likely not converged.")

      finalVals <- as.numeric(finalVals)
      names(finalVals) <- names(para)
      nll_pen_trace[i] <- next_nll_pen

      #replace any of the betas that are under the selection tolerance with 0, ignoring the baseline variables
      if(nP1+nP2+nP3 > 0){
        if(any(abs(finalVals[(1+nP0):nPtot]) < select_tol)) {
          finalVals[(1+nP0):nPtot][ abs(finalVals[(1+nP0):nPtot]) < select_tol] <- 0
        }
      }


      return(list(estimate=finalVals, niter=i,
                  final_ngrad_pen=as.numeric(temp_ngrad + n*Ek %*% finalVals),
                  final_nll=curr_nll,
                  final_nll_pen=curr_nll_pen,
                  nll_pen_trace=nll_pen_trace,
                  fit_code=7,
                  control=list(step_size_init=1,
                               step_size_min=step_size_min,
                               conv_crit=conv_crit,conv_tol=conv_tol,
                               maxit=maxit,mm_epsilon=mm_epsilon,
                               num_restarts=num_restarts),
                  # trace_mat=as.matrix(trace_mat),
                  # control=con,
                  final_step_size=lr,
                  restart_count=restart_count))

    }


    ##Take newton step and update values##
    ##**********************************##

    prevVals <- finalVals
    finalVals <- finalVals - lr*grad_step
    trace_mat <- cbind(trace_mat,finalVals)
    nll_pen_trace[i] <- next_nll_pen

    ##Check for convergence##
    ##*********************##

    max_change <- max(abs(finalVals-prevVals))
    scaled_est_change_norm <- sum(abs(finalVals-prevVals))/sum(abs(prevVals))
    nll_maj_change <- next_nll_maj-curr_nll_pen
    nll_pen_change <- next_nll_pen-curr_nll_pen
    prev_grad_mean_norm <- n^(-1)*sum(abs(temp_ngrad + n*Ek %*% prevVals))

    if(verbose){
      print(paste("max change in ests",max_change))
      print(paste("estimate with max change",names(para)[abs(finalVals-prevVals) == max(abs(finalVals-prevVals))]))
      print(paste("scaled change in l1 norm of ests", scaled_est_change_norm))
      print(paste("l1 mean norm of last gradient", prev_grad_mean_norm))
      print(paste("change in nll_maj", nll_maj_change))
      print(paste("new nll_maj", next_nll_maj))
      print(paste("new nll_pen", next_nll_pen))
    }

    if("maj_grad_norm" %in% conv_crit){
      if(prev_grad_mean_norm < conv_tol){
        fit_code <- 1
        break
      }
    }
    if("est_change_norm" %in% conv_crit){
      if(scaled_est_change_norm < conv_tol){
        fit_code <- 2
        break
      }
    }
    if("max_est_change" %in% conv_crit){
      if(max_change < conv_tol){
        fit_code <- 2
        break
      }
    }
    if("nll_pen_change" %in% conv_crit){
      if(abs(nll_pen_change) < conv_tol){
        fit_code <- 2
        break
      }
    }
    if("nll_maj_change" %in% conv_crit){
      if(abs(nll_maj_change) < conv_tol){
        fit_code <- 2
        break
      }
    }

    #iterate NR algorithm
    i <- i + 1
  }

  ##END ALGORITHM##
  ##*************##

  #if algorithm stopped because learning rate dropped too low, that is worth noting
  if(lr < step_size_min){
    fit_code <- 3
  }

  ##Compute traits of final estimates##
  ##*********************************##

  finalVals <- as.numeric(finalVals)
  names(finalVals) <- names(para)


  final_nll <- nll_func(para=finalVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                        hazard=hazard, frailty=frailty, model=model,
                        basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                        dbasis1=NULL, dbasis2=NULL, dbasis3=NULL)
  final_pen <- pen_func(para=finalVals,nP1=nP1,nP2=nP2,nP3=nP3,
                        penalty=penalty,lambda=lambda, a=a,
                        penweights_list=penweights_list, pen_mat_w_lambda = pen_mat_w_lambda)

  final_nll_pen <- final_nll + (n * final_pen)

  Ek <- pen_maj_mat_func(para=finalVals, nP1=nP1, nP2=nP2, nP3=nP3,
                         penalty=penalty,lambda=lambda, a=a,
                         penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                         penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                         penweights_list=penweights_list, mm_epsilon=mm_epsilon, hazard=hazard)

  #slight misnomer, but named this way for consistency with everything else
  final_ngrad_pen <- ngrad_func(para=finalVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                hazard=hazard, frailty=frailty, model=model,
                                basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                                dbasis1=NULL, dbasis2=NULL, dbasis3=NULL) +
    n*Ek %*% finalVals

  # final_nhess_pen <- nhess_func(para=finalVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
  #                               Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
  #                               hazard=hazard, frailty=frailty, model=model) +
  #   n*Ek

  ##Return list of final objects##
  ##****************************##

  #replace any of the betas that are under the selection tolerance with 0, ignoring the baseline variables
  if(nP1+nP2+nP3 > 0){
    if(any(abs(finalVals[(1+nP0):nPtot]) < select_tol)) {
      finalVals[(1+nP0):nPtot][ abs(finalVals[(1+nP0):nPtot]) < select_tol] <- 0
    }
  }

  return(list(estimate=finalVals, niter=i,
              final_ngrad_pen=as.numeric(final_ngrad_pen),
              final_nll=final_nll,
              final_nll_pen=final_nll_pen,
              nll_pen_trace=nll_pen_trace,
              fit_code=fit_code,
              control=list(step_size_init=1,
                           step_size_min=step_size_min,
                           conv_crit=conv_crit,conv_tol=conv_tol,
                           maxit=maxit,mm_epsilon=mm_epsilon,
                           num_restarts=num_restarts),
              # trace_mat=as.matrix(trace_mat),
              # control=con,
              final_step_size=lr,
              restart_count=restart_count))
}


##***********************************************##
####FUNCTION TO GENERATE SEQUENCES OF ESTIMATES####
##***********************************************##


#a function to loop through a path of lambda, lambda_fusedcoef, and mu_smooth values in a somewhat
#thoughtful way to maximize the pathwise connections between starting values and step sizes
#with some adjustments tailored to each approach to optimization
#ideally, the results of this function form the basis for the simulation outputs of interest
solution_path_function <- function(para, y1, y2, delta1, delta2,
                                   Xmat1 = matrix(nrow(length(y1)),ncol=0), #the default is a 'empty' matrix
                                   Xmat2 = matrix(nrow(length(y1)),ncol=0),
                                   Xmat3 = matrix(nrow(length(y1)),ncol=0),
                                   hazard, frailty, model,
                                   basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                                   dbasis1=NULL, dbasis2=NULL, dbasis3=NULL,
                                   penalty, lambda_path, a,
                                   penalty_fusedcoef, lambda_fusedcoef_path,
                                   penalty_fusedbaseline="none", lambda_fusedbaseline=0, #idk what to do about the baseline fusing, for now I skip
                                   penweights_list, mu_smooth_path, ball_R=Inf,
                                   fit_method, warm_start=TRUE,
                                   select_tol=1e-4, fusion_tol = 1e-3,
                                   step_size_min=1e-6, step_size_max=1e6,
                                   step_size_init=1,
                                   step_size_scale=0.5, #no checks implemented on these values!!
                                   step_delta=0.5,
                                   maxit=300,
                                   conv_crit = "omega",
                                   conv_tol=1e-6,
                                   verbose){

  # browser()

  lambda_path <- if(is.null(lambda_path)) as.matrix(0) else as.matrix(lambda_path)
  lambda_length <- nrow(lambda_path)
  colnames(lambda_path) <- paste0("lambda",1:ncol(lambda_path))

  lambda_fusedcoef_path <- if(is.null(lambda_fusedcoef_path)) as.matrix(0) else as.matrix(lambda_fusedcoef_path)
  lambda_fusedcoef_length <- nrow(lambda_fusedcoef_path)
  colnames(lambda_fusedcoef_path) <- paste0("lambda_fusedcoef",1:ncol(lambda_fusedcoef_path))

  mu_smooth_path <- if(is.null(mu_smooth_path)) as.matrix(0) else as.matrix(mu_smooth_path)
  mu_smooth_length <- nrow(mu_smooth_path)
  colnames(mu_smooth_path) <- "mu_smooth"


  n <- length(y1)
  Xmat1 <- if(!is.null(Xmat1)) as.matrix(Xmat1) else matrix(nrow=n,ncol=0)
  Xmat2 <- if(!is.null(Xmat2)) as.matrix(Xmat2) else matrix(nrow=n,ncol=0)
  Xmat3 <- if(!is.null(Xmat3)) as.matrix(Xmat3) else matrix(nrow=n,ncol=0)

  nPtot <- length(para)
  nP1 <- ncol(Xmat1)
  nP2 <- ncol(Xmat2)
  nP3 <- ncol(Xmat3)
  nP0 <- nPtot - nP1 - nP2 - nP3
  nP_vec <- c(nP0,nP1,nP2,nP3)


  #to correctly count number of iterations, we sum up number of iterations with no
  if(fit_method == "smooth_grad"){
    total_length <- lambda_length * lambda_fusedcoef_length * mu_smooth_length
  } else {
    total_length <- lambda_length*sum(lambda_fusedcoef_path %in% 0) +
      lambda_length*sum(!(lambda_fusedcoef_path %in% 0)) * mu_smooth_length
  }

  out_starts <- out_pen_ngrads <- out_ests <- matrix(nrow=total_length,
                                                     ncol=nPtot,dimnames=list(NULL,names(para)))
  out_ics <- matrix(nrow=total_length, ncol=15,
                    dimnames=list(NULL,c("nll", "nll_pen", "nPtot", "nPtot_selected", "nPtot_unique", "df_est",
                                         "AIC", "AIC_unique", "AIC_estdf",
                                         "BIC_unique", "BIC", "BIC_estdf",
                                         "GCV", "GCV_unique", "GCV_estdf")))
  out_info <- matrix(nrow=total_length,
                     ncol=ncol(lambda_path) + ncol(lambda_fusedcoef_path) + ncol(mu_smooth_path) + 1,
                     dimnames=list(NULL,c(colnames(lambda_path),
                                          colnames(lambda_fusedcoef_path),
                                          colnames(mu_smooth_path),"ball_R")))
  out_control <- matrix(nrow=total_length,
                        ncol=3,dimnames=list(NULL,c("niter","fit_code","final_step_size")))
  nll_pen_trace_mat <-  matrix(nrow=total_length,
                               ncol=maxit)

  # out_starts <- out_info <- out_pen_ngrads <- out_control <- out_ics <- out_ests <- NULL
  iter <- 1
  ##ACTUALLY BEGIN THE PATH##
  ##***********************##

  # browser()

  #keep running pathwise starting values and step sizes for the lambda, lambda_fused, and mu_smooth loops
  startVals_outer <- startVals_middle <- startVals_inner <- para
  step_size_init_outer <- step_size_init_middle <- step_size_init_inner <- step_size_init

  #OUTER LOOP: LOOPING THROUGH THE LAMBDA VALUES, GOVERNING PARAMETERWISE SPARSITY
  for(lambda_iter in 1:lambda_length){
    lambda <- lambda_path[lambda_iter,]
    # conv_tol <- if(lambda_iter==lambda_length) lambda/8 else lambda/4 #gotta be less than lambda_target/4, so we just halve it again.

    startVals_middle <- startVals_outer
    step_size_init_middle <- step_size_init_outer

    #MIDDLE LOOP: LOOPING THROUGH THE LAMBDA_FUSEDCOEF VALUES, GOVERNING FUSION SPARSITY
    for(lambda_fusedcoef_iter in 1:lambda_fusedcoef_length){
      lambda_fusedcoef <- lambda_fusedcoef_path[lambda_fusedcoef_iter,]

      startVals_inner <- startVals_middle
      step_size_init_inner <- step_size_init_middle

      if(verbose){
        print(paste0("step: ",lambda_iter, " lambda: ",lambda," lambda_fusedcoef: ",lambda_fusedcoef))#, " mu_smooth_fused: ",mu_smooth_fused))
      }

      #INNER LOOP: LOOPING THROUGH THE MU_SMOOTH VALUES, GOVERNING SMOOTH APPROXIMATION OF FUSION SPARSITY
      for(mu_smooth_iter in 1:mu_smooth_length){

        #if there is no fusion, there should be no smoothing unless we're explicitly smoothing the parameterwise penalties
        mu_smooth_fused <- if(fit_method != "smooth_grad" & lambda_fusedcoef==0) 0 else mu_smooth_path[mu_smooth_iter,]

        if(fit_method=="prox_grad"){
          tempfit <- proximal_gradient_descent(para=startVals_inner, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                               Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                               hazard=hazard, frailty=frailty, model=model,
                                               basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                               dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                               penalty=penalty, lambda=lambda, a=a,
                                               penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                               penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                               penweights_list=penweights_list, mu_smooth=0, mu_smooth_fused=mu_smooth_fused,
                                               step_size_init=step_size_init_inner,step_size_min=step_size_min,step_size_max=step_size_max,step_size_scale=step_size_scale,
                                               maxit=maxit, conv_crit=conv_crit, conv_tol=conv_tol,
                                               ball_R=Inf, verbose=verbose)
        }
        if(fit_method=="prox_grad_nmaccel"){
          tempfit <- proximal_gradient_descent_nmaccel(para=startVals_inner, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                                       Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                                       hazard=hazard, frailty=frailty, model=model,
                                                       basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                                       dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                                       penalty=penalty, lambda=lambda, a=a,
                                                       penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                                       penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                                       penweights_list=penweights_list, mu_smooth=0, mu_smooth_fused=mu_smooth_fused,
                                                       step_size_init_x=step_size_init_inner,step_size_init_y=step_size_init_inner,
                                                       step_size_min=step_size_min,step_size_max=step_size_max,step_size_scale=step_size_scale,
                                                       maxit=maxit, conv_crit=conv_crit, conv_tol=conv_tol,
                                                       ball_R=Inf, verbose=verbose)
        }
        if(fit_method=="smooth_grad"){
          tempfit <- smooth_gradient_descent(para=startVals_inner, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                             Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                             hazard=hazard, frailty=frailty, model=model,
                                             basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                             dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                             penalty=penalty, lambda=lambda, a=a,
                                             penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                             penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                             penweights_list=penweights_list, mu_smooth=mu_smooth_fused, mu_smooth_fused=mu_smooth_fused,
                                             step_size_init=step_size_init_inner,
                                             step_size_min=step_size_min,step_size_max=step_size_max,step_size_scale=step_size_scale,
                                             maxit=maxit, conv_crit=conv_crit, conv_tol=conv_tol, method="grad",#method="BFGS", #method is currently assumed to be standard gradient
                                             # ball_R=Inf, #currently not supported because it would require the prox function again.
                                             verbose=verbose)
        }
        if(fit_method=="newton"){
          tempfit <- newton_raphson_mm(para=startVals_inner, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                       Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                       hazard=hazard, frailty=frailty, model=model,
                                       penalty=penalty, lambda=lambda, a=a,
                                       penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                       penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                       penweights_list=penweights_list,
                                       mm_epsilon=1e-6,
                                       verbose=verbose, maxit=maxit,step_size_min=step_size_min,
                                       conv_crit=conv_crit,conv_tol=conv_tol,num_restarts=2)
        }

        #To prepare for the next iteration, let's update the inner loop inputs based on these results
        startVals_inner <- tempfit$estimate
        if(fit_method %in%  c("prox","smooth_grad")){
          step_size_init_inner <- tempfit$final_step_size
        }

        ##Now, let's COMPUTE SOME USEFUL FIT STATISTICS based on this fit##
        ##***************************************************************##
        final_nll <- tempfit$final_nll

        beta1_selected <- if(nP1 != 0) startVals_inner[(1+nP0):(nP0+nP1)] else numeric(0)
        nP1_selected <- sum(beta1_selected != 0)

        beta2_selected <- if(nP2 != 0) startVals_inner[(1+nP0+nP1):(nP0+nP1+nP2)] else numeric(0)
        nP2_selected <- sum(beta2_selected != 0)

        beta3_selected <- if(nP3 != 0) startVals_inner[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)] else numeric(0)
        nP3_selected <- sum(beta3_selected != 0)

        nPtot_selected <- nP0 + nP1_selected + nP2_selected + nP3_selected

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
        if(!is.null(tempfit$final_nhess)){
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

        out_info[iter,] <- c(lambda,lambda_fusedcoef, mu_smooth_fused, ball_R)
        out_ests[iter,] <- tempfit$estimate
        out_pen_ngrads[iter,] <- as.vector(tempfit$final_ngrad_pen)
        out_ics[iter,] <- c(tempfit$final_nll, tempfit$final_nll_pen,
                            nPtot, nPtot_selected, nPtot_unique, df_est,
                            tempaic, tempaic_unique,tempaic_estdf,
                            tempbic_unique,tempbic,tempbic_estdf,
                            tempgcv, tempgcv_unique,tempgcv_estdf)
        out_control[iter,] <- c(tempfit$niter,tempfit$fit_code,tempfit$final_step_size)
        out_starts[iter,] <- startVals_inner
        nll_pen_trace_mat[iter, ] <- tempfit$nll_pen_trace
        # out_info <- rbind(out_info,c(lambda,lambda_fusedcoef, mu_smooth_fused=mu_smooth_fused, ball_R=ball_R))
        # out_ests <- rbind(out_ests,tempfit$estimate)
        # out_pen_ngrads <- rbind(out_pen_ngrads,as.vector(tempfit$final_ngrad_pen))
        # out_ics <- rbind(out_ics,c(nll=tempfit$final_nll, nll_pen=tempfit$final_nll_pen,
        #                            nPtot_selected=nPtot_selected,
        #                            nPtot_unique=nPtot_unique, df_est=df_est,
        #                            AIC=tempaic, AIC_unique=tempaic_unique,AIC_estdf=tempaic_estdf,
        #                            BIC_unique=tempbic_unique,BIC=tempbic,BIC_estdf=tempbic_estdf,
        #                            GCV=tempgcv, GCV_unique=tempgcv_unique,GCV_estdf=tempgcv_estdf))
        # out_control <- rbind(out_control,c(niter=tempfit$niter,fit_code=tempfit$fit_code,final_step_size=tempfit$final_step_size))
        # out_starts <- rbind(out_starts,startVals_inner)


        #if this is the first time through this middle/inner combo,
        #then the results are relevant for the next iteration of the outer loop
        if(warm_start & lambda_fusedcoef_iter==1 & mu_smooth_iter==1){
          startVals_outer <- tempfit$estimate
          if(fit_method %in% c("prox","smooth_grad")){
            step_size_init_outer <- tempfit$final_step_size
          }
        }

        #if this is the first time through this inner loop,
        #then the results are relevant for the next iteration of the middle loop
        if(warm_start & mu_smooth_iter==1){
          if(fit_method %in% c("prox","smooth_grad")){
            step_size_init_middle <- tempfit$final_step_size
          }
        }


        iter <- iter + 1
        #if there is no smoothing, then break out of the inner smoothing loop
        if(mu_smooth_fused==0){
          break
        }

      }

      #now that we've reached the end of the inner loop, let's use those starting values for the next middle loop step
      #we have already above determined the correct starting step size for the next iteration of the middle loop
      if(warm_start){
        startVals_middle <- tempfit$estimate
      }

    }
  }

  # browser()

  ##MAKE SOME SUMMARY PLOTS AND FINALIZE THE OUTPUT TABLES##
  ##******************************************************##

  rownames(out_info) <- NULL
  rownames(out_ics) <- NULL
  rownames(out_ests) <- NULL
  rownames(out_pen_ngrads) <- NULL
  rownames(out_starts) <- NULL
  rownames(out_control) <- NULL
  rownames(nll_pen_trace_mat) <- NULL
  # colnames(out_control) <- c("niter","fitcode","final_step_size")
  # out_ests <- as.data.frame(cbind(out_info,out_ests))
  # out_ics <- as.data.frame(cbind(out_info,out_ics))
  # out_control <- as.data.frame(cbind(out_info,out_control))
  # out_pen_ngrads <- as.data.frame(cbind(out_info,out_pen_ngrads))
  # out_starts <- as.data.frame(cbind(out_info,out_starts))
  # bic_best_index <- which(out_ics$BIC == min(out_ics$BIC))
  # print(out_ests)
  # if(ncol(lambda_path) == 1){ #as long as all the fused lasso lambdas are the same, make the plots
  #   # colnames(out_ests) <- c("lambda",names(tempfit$estimate))
  #   plot_out_ests <- reshape2::melt(out_ests, id="lambda1") %>%
  #     dplyr::mutate(vartype=ifelse(stringr::str_detect(variable,"_1"),"Beta1",
  #                                  ifelse(stringr::str_detect(variable,"_2"),"Beta2",
  #                                         ifelse(stringr::str_detect(variable,"_3"),"Beta3",
  #                                                "Weibull Params"))))
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

  return(list(out_info=out_info,
              out_ests=out_ests,
              out_pen_ngrads=out_pen_ngrads,
              out_ics=out_ics,
              out_starts=out_starts,
              nll_pen_trace_mat=nll_pen_trace_mat,
              # step_eta=step_eta,
              # plot_out=plot_out,
              # plot_out_ests=plot_out_ests,
              out_control=out_control,
              hazard=hazard, frailty=frailty, model=model,
              penalty=penalty, a=a, mm_epsilon=1e-6,
              select_tol=select_tol, fusion_tol=fusion_tol,
              penalty_fusedcoef=penalty_fusedcoef,
              penalty_fusedbaseline=penalty_fusedbaseline,
              lambda_fusedbaseline=lambda_fusedbaseline,
              fit_method=fit_method))


}

