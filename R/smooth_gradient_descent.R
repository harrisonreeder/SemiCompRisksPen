smooth_gradient_descent <- function(para, y1, y2, delta1, delta2,
                                    Xmat1 = matrix(nrow(length(y1)),ncol=0), #the default is a 'empty' matrix
                                    Xmat2 = matrix(nrow(length(y1)),ncol=0),
                                    Xmat3 = matrix(nrow(length(y1)),ncol=0),
                                    hazard, frailty, model,
                                    basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                                    dbasis1=NULL, dbasis2=NULL, dbasis3=NULL,
                                    penalty, lambda, a,
                                    penalty_fusedcoef, lambda_fusedcoef,
                                    penalty_fusedbaseline, lambda_fusedbaseline,
                                    penweights_list, mu_smooth, mu_smooth_fused,
                                    step_size_init=1, step_size_min = 1e-6, step_size_max = 1e6,
                                    step_size_scale=1/2,  maxit=300, #currently this doesn't support ball_R constraint
                                    conv_crit = "omega", conv_tol=if(lambda>0) lambda/4 else 1e-6,
                                    method="BFGS", verbose){
  #START HERE AND UPDATE THE REST OF THE FUNCTION WITH THE NEW SYNTAX

  # browser()
  ##Set up control pieces##
  ##*********************##

  #construct control list
  #conv_crit determines criterion used to judge convergence
  #'est_change_norm' looks at l1 norm of change in parameter values (scaled by l1 norm of prior values): sum(abs(finalVals-prevVals))/sum(abs(prevVals))
  #'max_est_change' looks at largest absolute change in a parameter value
  #'nll_pen_change' looks at change in regularized log-likelihood
  #'maj_grad_norm' looks at l1 norm of majorized gradient
  if("maxit" <0){
    stop("maxit must be nonnegative.")
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

  ##CREATE CONTRAST MATRICES AND DECOMPOSITIONS FOR LATER USE IN PROXIMAL ALGORITHMS##
  ##I HAD PREVIOUSLY DEVISED A TWO-STAGE CONSTRUCTION TO MAKE UNWEIGHTED VERSION, COMPUTE WEIGHTS, AND THEN MAKE WEIGHTED VERSION
  ##BUT THEN I DECIDED THAT FOR NOW, COMPUTING WEIGHTS MIGHT HAPPEN OUTSIDE OF THIS FUNCTION
  #create list of unweighted contrast matrices
  # pen_mat_list <- contrast_mat_list(nP0=nP0,nP1 = nP1,nP2 = nP2,nP3 = nP3,
  #                                      penalty_fusedcoef = penalty_fusedcoef, lambda_fusedcoef = lambda_fusedcoef,
  #                                      penalty_fusedbaseline = penalty_fusedbaseline,lambda_fusedbaseline = lambda_fusedbaseline,
  #                                      hazard = hazard, penweights_list = NULL)
  # for(pen_name in names(D_list_noweight)){
  #create list of weights based on adaptive lasso inverse contrasts
  #   penweights_list[[pen_name]] <- penweights_internal(parahat=mle_optim,
  #                                                                D=pen_mat_list[[var]],
  #                                                                penweight_type="adaptive",addl_penweight=NULL)
  # }

  #create list of (adaptive) weighted contrast matrices
  #if there is no fusion, this should just be an empty list
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

  #compute eigenvector decomposition of this weighted contrast matrix
  #if there is no fusion, this should return a list with Q = as.matrix(0) and eigval = 0
  # pen_mat_w_eig <- pen_mat_decomp(pen_mat_w)
  pen_mat_w_eig <- NULL

  ##Run algorithm##
  ##*************##
  if(method %in% c("L-BFGS-B","BFGS")){
    out <- optim(par = para, fn=smooth_obj_func, gr=smooth_obj_grad_func,
                 y1=y1, y2=y2, delta1=delta1, delta2=delta2, Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                 hazard=hazard, frailty=frailty, model=model,
                 basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                 dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                 penalty=penalty, lambda=lambda, a=a,
                 penweights_list=penweights_list, mu_smooth=mu_smooth,
                 pen_mat_w_lambda=pen_mat_w_lambda, mu_smooth_fused=mu_smooth_fused,
                 method=method, control=list(maxit=maxit))
    finalVals <- out$par
    final_nll_smooth <- out$value
    niter <- out$count[2] #records the number of gradients computed, which is basically niter
    fit_code <- out$convergence

    #Here, report final nll on the SUM scale! No division by n
    final_nll <- nll_func(para=finalVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                          Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                          hazard=hazard, frailty=frailty, model=model,
                          basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                          dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
    final_pen <- pen_func(para=finalVals,nP1=nP1,nP2=nP2,nP3=nP3,
                          penalty=penalty,lambda=lambda, a=a,penweights_list=penweights_list,
                          pen_mat_w_lambda = pen_mat_w_lambda) #Here, we're reporting the final results under the un-smoothed fused lasso

    final_nll_pen <- final_nll + (n*final_pen)

    #note the division by n to put it on the mean scale--gradient descent works better then!
    final_ngrad <- smooth_obj_grad_func(para=finalVals,
                                        y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                        hazard=hazard, frailty=frailty, model=model,
                                        basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                        dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                        penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                        mu_smooth=mu_smooth, pen_mat_w_lambda = pen_mat_w_lambda,
                                        mu_smooth_fused = mu_smooth_fused) #here, we do still report the nesterov-smoothed results


    return(list(estimate=finalVals, niter=niter, fit_code=fit_code,
                final_nll=final_nll, final_nll_pen=final_nll_pen,
                final_nll_smooth=final_nll_smooth,
                final_ngrad=final_ngrad,
                final_ngrad_pen=final_ngrad, #because we don't have a prox in this case, just use regular grad
                startVals=para, hazard=hazard, frailty=frailty, model=model,
                penalty=penalty, lambda=lambda, a=a,
                penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                mu_smooth=mu_smooth,
                mu_smooth_fused=mu_smooth_fused,
                # trace_mat=as.matrix(trace_mat),
                nll_pen_trace=final_nll_pen,
                method=method,
                ball_R=Inf, #for this one we don't have the projection operator set up yet.
                control=NULL,
                final_step_size=NULL,
                bad_step_count=NULL))

  }




  nll_pen_trace <- rep(NA,maxit)
  trace_mat <- xcurr <- para
  fit_code <- 4 #follows convention from 'nleqslv' L-BFGS that is 1 if converged, 4 if maxit reached, 3 if stalled
  bad_step_count <- restart_count <- 0


  #note the division by n to put it on the mean scale--gradient descent works better then!
  obj_val_xcurr <- smooth_obj_func(para=xcurr,
                                   y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                   Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                   hazard=hazard, frailty=frailty, model=model,
                                   basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                   dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                   penalty=penalty,lambda=lambda, a=a,
                                   penweights_list=penweights_list,
                                   pen_mat_w_lambda = pen_mat_w_lambda,
                                   mu_smooth=mu_smooth,
                                   mu_smooth_fused = mu_smooth_fused)/n
  if(is.na(obj_val_xcurr)){stop("Initial values chosen yield infinite likelihood values. Consider other initial values.")}

  #note the division by n to put it on the mean scale--gradient descent works better then!
  ngrad_xcurr <- smooth_obj_grad_func(para=xcurr,
                                      y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                      hazard=hazard, frailty=frailty, model=model,
                                      basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                      dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                      penalty=penalty,lambda=lambda, a=a,
                                      penweights_list=penweights_list, mu_smooth=mu_smooth,
                                      pen_mat_w_lambda = pen_mat_w_lambda,
                                      mu_smooth_fused = mu_smooth_fused)/n


  #constants
  #define things using notation from Wang (2014), which also covers Zhao (2016) and Zhao (2018)
  #outputs
  #monitors

  #I'M TESTING SOMETHING HERE
  # step_L_ti <- step_L
  step_size_ti <- step_size_init #this is 1/L in the notation of Wang (2014)


  i <- 1 #this is instead of 'k' in Wang (2014)
  while(i <= maxit){
    if(verbose)print(i)
    # if(i==15){browser()}
    ##RUN ACTUAL STEP##
    ##***************##

    #slightly increase step size at each step (though, with line search this might get knocked back down.)
    step_size_ti <- min(step_size_max,step_size_ti/step_size_scale)

    #run prox function:
    #ADMM prox subroutine for fused lasso only if mu_smooth_fused=0
    #soft thresholding according to lambda*step_size*weight
    xnext <- xcurr-ngrad_xcurr * step_size_ti

    #note the division by n to put it on the mean scale--gradient descent works better then!
    obj_val_xnext <- smooth_obj_func(para=xnext,
                                     y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                     Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                     hazard=hazard, frailty=frailty, model=model,
                                     basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                     dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                     penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                     pen_mat_w_lambda = pen_mat_w_lambda,
                                     mu_smooth = mu_smooth,
                                     mu_smooth_fused = mu_smooth_fused)/n

    #in this case, this is just a smooth lqa , there is no additional penalty because everything is just smoothed
    # smooth_obj_lqa_pen_xnext <- smooth_obj_lqa_pen_func(para=xnext, prev_para=xcurr,
    #                                                     y1=y1, y2=y2, delta1=delta1, delta2=delta2,
    #                                                     Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
    #                                                     hazard=hazard, frailty=frailty, model=model,
    #                                                     basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
    #                                                     dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
    #                                                     penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
    #                                                     pen_mat_w_lambda = pen_mat_w_lambda,
    #                                                     mu_smooth=mu_smooth, mu_smooth_fused = mu_smooth_fused,
    #                                                     step_size=step_size_ti)/n

    if(is.nan(obj_val_xnext)){
      if(verbose)print("whoops, proposed step yielded infinite likelihood value, just fyi")
      obj_val_xnext <- smooth_obj_lqa_pen_xnext <- Inf
    }

    while(is.infinite(obj_val_xnext) ||
          obj_val_xnext >
          # smooth_obj_lqa_pen_xnext
          obj_val_xcurr - 0.5*step_size_ti*sum(ngrad_xcurr^2)
    ){
      if(step_size_ti < step_size_min){
        if(verbose)print("step size got too small, accepting result anyways")
        bad_step_count <- bad_step_count + 1
        break
      }
      step_size_ti <- step_size_ti * step_size_scale
      # if(verbose)print(paste0("nll_pen_xnext:",obj_val_xnext," bigger than smooth_obj_lqa_pen_xnext:",smooth_obj_lqa_pen_xnext))
      if(verbose)print(paste0("nll_pen_xnext:",obj_val_xnext," bigger than backtrack boundary:",obj_val_xcurr - 0.5*step_size_ti*sum(ngrad_xcurr^2)))
      if(verbose)print(paste0("effective step size reduced to: ",step_size_ti))

      xnext <- xcurr - ngrad_xcurr * step_size_ti

      #note the division by n to put it on the mean scale--gradient descent works better then!
      obj_val_xnext <-  smooth_obj_func(para=xnext,
                                        y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                        hazard=hazard, frailty=frailty, model=model,
                                        basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                        dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                        penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                        pen_mat_w_lambda = pen_mat_w_lambda,
                                        mu_smooth = mu_smooth,
                                        mu_smooth_fused = mu_smooth_fused)/n

      # smooth_obj_lqa_pen_xnext <- smooth_obj_lqa_pen_func(para=xnext, prev_para=xcurr,
      #                                                     y1=y1, y2=y2, delta1=delta1, delta2=delta2,
      #                                                     Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
      #                                                     hazard=hazard, frailty=frailty, model=model,
      #                                                     basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
      #                                                     dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
      #                                                     penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
      #                                                     mu_smooth=mu_smooth, pen_mat_w_lambda = pen_mat_w_lambda,
      #                                                     mu_smooth_fused = mu_smooth_fused,step_size=step_size_ti)/n

      if(is.nan(obj_val_xnext)){
        if(verbose)print("whoops, proposed step yielded infinite likelihood value, just fyi")
        obj_val_xnext <- smooth_obj_lqa_pen_xnext <- Inf
      }
    }

    ##UPDATE MONITORS##
    ##***************##

    xprev <- xcurr
    xcurr <- xnext

    # trace_mat <- cbind(trace_mat,xnext)

    # nll_pen_trace[i] <- nll_pen_xnext
    obj_val_xnext <-  smooth_obj_func(para=xnext,
                                      y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                      hazard=hazard, frailty=frailty, model=model,
                                      basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                      dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                      penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                      pen_mat_w_lambda = pen_mat_w_lambda,
                                      mu_smooth = mu_smooth,
                                      mu_smooth_fused = mu_smooth_fused)/n

    nll_pen_trace[i] <- nll_pen_func(para=xnext,
                                     y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                     Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                     basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                     dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                     hazard=hazard, frailty=frailty, model=model,
                                     penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                     pen_mat_w_lambda = pen_mat_w_lambda,mu_smooth_fused = 0)/n #RECORD RESULT WITHOUT SMOOTHING, TO PUT US ON A COMMON SCALE!

    ##Check for convergence##
    ##*********************##

    #note the division by n to put it on the mean scale--gradient descent works better then!
    ngrad_xcurr <- smooth_obj_grad_func(para=xcurr,
                                        y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                        hazard=hazard, frailty=frailty, model=model,
                                        basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                        dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                        penalty=penalty,lambda=lambda, a=a,
                                        penweights_list=penweights_list, mu_smooth=mu_smooth,
                                        pen_mat_w_lambda = pen_mat_w_lambda,mu_smooth_fused = mu_smooth_fused)/n

    #convergence criterion given in Wang (2014) adapted without prox operator
    omega_t <- max(abs(ngrad_xcurr))

    max_change <- max(abs(xcurr-xprev))
    norm_change <- sqrt(sum((xcurr-xprev)^2))
    obj_val_change <- obj_val_xnext-obj_val_xcurr

    if(verbose){
      print(paste("max change in ests",max_change))
      print(paste("omega_t (max norm of prox grad)", omega_t))
      print(paste("estimate with max change",names(para)[abs(xcurr-xprev) == max_change]))
      print(paste("max norm of change", norm_change)) #essentially a change in estimates norm
      print(paste("change in obj_val", obj_val_change))
      print(paste("new obj_val", obj_val_xnext))
    }

    if("omega" %in% conv_crit){
      if(is.null(omega_t) || !is.finite(omega_t)){
        warning("maximum absolute gradient value is not finite. Exiting with fit_code 5, likely not converged.")
        fit_code <- 5
        break
      }
      if(omega_t <= conv_tol){
        fit_code <- 2
        break
      } #conv_tol is called 'epsilon' in the Wang (2014) paper
    }
    if("est_change_norm" %in% conv_crit){
      if(norm_change < conv_tol){
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
      if(abs(obj_val_change) < conv_tol){
        fit_code <- 2
        break
      }
    }

    #after the updates, now update the nll monitor values
    obj_val_xcurr <- obj_val_xnext
    #iterate algorithm
    i <- i + 1
  }

  ##END ALGORITHM##
  ##*************##

  #if algorithm stopped because learning rate dropped too low, that is worth noting
  # if(lr < con[["min_lr"]]){
  #   fit_code <- 3
  # }

  ##Compute traits of final estimates##
  ##*********************************##

  finalVals <- as.numeric(xnext)
  names(finalVals) <- names(para)

  #Here, report final nll on the SUM scale! No division by n
  final_nll <- nll_func(para=finalVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                        hazard=hazard, frailty=frailty, model=model,
                        basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                        dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
  final_pen <- pen_func(para=finalVals,nP1=nP1,nP2=nP2,nP3=nP3,
                        penalty=penalty,lambda=lambda, a=a,penweights_list=penweights_list,
                        pen_mat_w_lambda = pen_mat_w_lambda) #Here, we're reporting the final results under the un-smoothed fused lasso

  final_nll_pen <- final_nll + (n*final_pen)

  final_nll_smoothed <- smooth_obj_func(para=finalVals,
                                        y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                        hazard=hazard, frailty=frailty, model=model,
                                        basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                        dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                        penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                        mu_smooth=mu_smooth, pen_mat_w_lambda = pen_mat_w_lambda,
                                        mu_smooth_fused = mu_smooth_fused)


  #note the division by n to put it on the mean scale--gradient descent works better then!
  final_ngrad <- smooth_obj_grad_func(para=finalVals,
                                      y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                      hazard=hazard, frailty=frailty, model=model,
                                      basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                      dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                      penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                      mu_smooth=mu_smooth, pen_mat_w_lambda = pen_mat_w_lambda,
                                      mu_smooth_fused = mu_smooth_fused) #here, we do still report the nesterov-smoothed results

  final_ngrad_pen <- prox_func(para=final_ngrad, prev_para = xprev,
                               nP1=nP1,nP2=nP2,nP3=nP3,step_size=1,
                               penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                               pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig,
                               lambda_f_vec=lambda_f_vec, mu_smooth=mu_smooth,
                               mu_smooth_fused = mu_smooth_fused,
                               ball_R=Inf)
  ##Return list of final objects##
  ##****************************##

  return(list(estimate=finalVals, niter=i, fit_code=fit_code,
              final_nll=final_nll, final_nll_pen=final_nll_pen,
              final_nll_smoothed=final_nll_smoothed,
              final_ngrad=final_ngrad, final_ngrad_pen=final_ngrad_pen,
              startVals=para, hazard=hazard, frailty=frailty, model=model,
              penalty=penalty, lambda=lambda, a=a,
              penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
              penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
              mu_smooth=mu_smooth,
              mu_smooth_fused=mu_smooth_fused,
              # trace_mat=as.matrix(trace_mat),
              nll_pen_trace = n*nll_pen_trace,
              method="grad",
              control=list(step_size_init=step_size_init,
                           step_size_min=step_size_min,
                           step_size_max=step_size_max,
                           step_size_scale=step_size_scale,
                           conv_crit=conv_crit,conv_tol=conv_tol),
              ball_R=Inf, #for this one we don't have the projection operator set up yet.
              final_step_size=step_size_ti,
              bad_step_count=bad_step_count))
}
