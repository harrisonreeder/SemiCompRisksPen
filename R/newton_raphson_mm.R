#' Newton-Type MM Algorithm
#'
#' This function runs a Newton-type MM (majorization-minimization) algorithm
#'   following the algorithm of Oelker & Tutz (2017).
#'
#' @inheritParams proximal_gradient_descent
#' @param mm_epsilon Positive numeric tolerance parameter for smooth approximation
#'   of absolute value function at 0.
#' @param num_restarts Number of times to allow algorithm to restart if it reaches
#'   a point where it can make no further progress.
#'
#' @return A list.
#' @export
newton_raphson_mm <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                              hazard, frailty, model,
                              penalty, lambda, a,
                              penalty_fusedcoef, lambda_fusedcoef,
                              penalty_fusedbaseline, lambda_fusedbaseline,
                              penweights_list=penweights_list, mm_epsilon, verbose,
                              maxit=300,step_size_min=0.00001,
                              conv_crit="nll_pen_change",conv_tol=1e-05,num_restarts=2,
                              select_tol=1e-4){


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



  ##Set up control pieces##
  ##*********************##
  #browser()
  #construct control list
  #maxit is maximum number of newton-raphson iterations allowed
  #min_lr is the minimum learning rate allowed through step-halving process.
  #conv_crit determines criterion used to judge convergence
  #est_change_norm looks at l1 norm of change in parameter values (scaled by l1 norm of prior values): sum(abs(finalVals-prevVals))/sum(abs(prevVals))
  #est_change_max looks at largest absolute change in a parameter value
  #ll_maj_change looks at change in majorized log-likelihood (expanded around previous location)
  #ll_reg_change looks at change in regularized log-likelihood
  #maj_grad_norm looks at l1 norm of majorized gradient
  #if (length(noNms <- namc[!namc %in% nmsC])){
  #warning("unknown names in control: ", paste(noNms, collapse=", "))
  #}
  #conv_crit determines criterion used to judge convergence
  #est_change_2norm' looks at l1 norm of change in parameter values (scaled by l1 norm of prior values): sum(abs(finalVals-prevVals))/sum(abs(prevVals))
  #est_change_max' looks at largest absolute change in a parameter value
  #nll_pen_change' looks at change in regularized log-likelihood
  #maj_grad_norm' looks at l1 norm of majorized gradient
  if(!all(conv_crit %in% c("est_change_2norm","est_change_max","nll_pen_change","nll_maj_change","suboptimality"))){
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

  #THERE HAS GOT TO BE A BETTER WAY TO INTEGRATE THIS INTO THE MM ALGORITHM, BECAUSE AT PRESENT THE FUNCTIONS REALLY DUPLICATE A LOT OF EFFORT, BUT ITS SOMETHING


  ##Run algorithm##
  ##*************##

  nll_pen_trace <- rep(NA,maxit)
  trace_mat <- prevVals <- finalVals <- startVals <- para
  fit_code <- 4 #follows convention from nleqslv L-BFGS that is 1 if converged, 4 if maxit reached, 3 if stalled
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
      trace_mat <- prevVals <- finalVals <- startVals <- stats::runif(n=1,min=0.8,max=1.2)*para
      fit_code <- 4 #follows convention from nleqslv that is 1 if converged, 4 if maxit reached, 3 if stalled
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
      trace_mat <- prevVals <- finalVals <- startVals <- stats::runif(n=1,min=0.8,max=1.2)*para
      fit_code <- 4 #follows convention from nleqslv that is 1 if converged, 4 if maxit reached, 3 if stalled
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
      trace_mat <- prevVals <- finalVals <- startVals <- stats::runif(n=1,min=0.8,max=1.2)*para
      fit_code <- 4 #follows convention from nleqslv that is 1 if converged, 4 if maxit reached, 3 if stalled
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

    est_change_max <- max(abs(finalVals-prevVals))
    est_change_2norm <- sqrt(sum((finalVals-prevVals)^2))
    # scaled_est_change_norm <- sum(abs(finalVals-prevVals))/sum(abs(prevVals))
    nll_maj_change <- next_nll_maj-curr_nll_pen
    nll_pen_change <- next_nll_pen-curr_nll_pen
    # prev_grad_mean_norm <- n^(-1)*sum(abs(temp_ngrad + n*Ek %*% prevVals))
    subopt_t <- n^(-1)*max(abs(temp_ngrad + n*Ek %*% prevVals))

    if(verbose){
      print(paste("max change in ests",est_change_max))
      print(paste("estimate with max change",names(para)[abs(finalVals-prevVals) == est_change_max]))
      print(paste("l2 norm of estimate change",est_change_2norm))
      # print(paste("scaled change in l1 norm of ests", scaled_est_change_norm))
      # print(paste("l1 mean norm of last gradient", prev_grad_mean_norm))
      print(paste("change in nll_maj", nll_maj_change))
      print(paste("new nll_maj", next_nll_maj))
      print(paste("new nll_pen", next_nll_pen))
    }

    if("suboptimality" %in% conv_crit){
      if(subopt_t < conv_tol){
        fit_code <- 2
        break
      }
    }
    if("est_change_2norm" %in% conv_crit){
      if(est_change_2norm < conv_tol){
        fit_code <- 2
        break
      }
    }
    if("est_change_max" %in% conv_crit){
      if(est_change_max < conv_tol){
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
