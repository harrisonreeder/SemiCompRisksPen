#' Proximal Gradient Descent Algorithm
#'
#' This function runs a proximal gradient descent algorithm with backtracking similar to that presented by
#'   Wang et al. (2014).
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
#' @param penalty A string value indicating the form of parameterwise penalty
#'   to apply. "lasso", "scad", and "mcp" are the options.
#' @param lambda The strength of the parameterwise penalty. Either a single non-negative numeric value
#'   for all three transitions, or a length 3 vector with elements corresponding to the three transitions.
#' @param a For two-parameter penalty functions (e.g., scad and mcp), the second parameter.
#' @param penalty_fusedcoef A string value indicating the form of the fusion penalty to apply
#'   to the regression parameters. "none" and "fusedlasso" are the options.
#' @param lambda_fusedcoef The strength of the fusion penalty on the regression parameters.
#'   Either a single non-negative numeric value
#'   for all three transitions, or a length 3 vector with elements corresponding to the three transitions.
#' @param penalty_fusedbaseline A string value indicating the form of the fusion penalty to apply
#'   to the baseline hazard parameters. "none" and "fusedlasso" are the options.
#' @param lambda_fusedbaseline The strength of the fusion penalty on the regression parameters.
#'   Either a single non-negative numeric value
#'   for all three transitions, or a length 3 vector with elements corresponding to the three transitions.
#' @param penweights_list A list of numeric vectors representing weights for each
#'   penalty term (e.g., for adaptive lasso.) Elements of the list should be indexed by the
#'   names "coef1", "coef2", "coef3", "fusedcoef12", "fusedcoef13", "fusedcoef23", "fusedbaseline12", "fusedbaseline13", and "fusedbaseline23"
#' @param mu_smooth A non-negative numeric value for the Nesterov smoothing parameter applied to the parameterwise penalty
#' @param mu_smooth_fused A non-negative numeric value for the Nesterov smoothing parameter applied to the fusion penalty.
#' @param step_size_init Positive numeric value for the initial step size.
#' @param step_size_min Positive numeric value for the minimum allowable step size to allow during backtracking.
#' @param step_size_max Positive numeric value for the maximum allowable step size to allow by size increase at each iteration.
#' @param step_size_scale Positive numeric value for the multiplicative change in step size at each step of backtracking.
#' @param ball_R Positive numeric value for \eqn{l_2} ball constraint around the origin for the regression parameters.
#'   Typically set to \code{Inf} indicating no constraint, otherwise equivalent to an extra \eqn{l_2} penalty.
#' @param maxit Positive integer maximum number of iterations.
#' @param conv_crit String (possibly vector) giving the convergence criterion.
#' @param conv_tol Positive numeric value giving the convergence tolerance for the chosen criterion.
#' @param verbose Boolean indicating whether information about each iteration should be printed.
#'
#' @return A list.
#' @export
proximal_gradient_descent <- function(para, y1, y2, delta1, delta2,
                                      Xmat1 = matrix(nrow(length(y1)),ncol=0), #the default is a 'empty' matrix
                                      Xmat2 = matrix(nrow(length(y1)),ncol=0),
                                      Xmat3 = matrix(nrow(length(y1)),ncol=0),
                                      hazard, frailty, model,
                                      basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                                      dbasis1=NULL, dbasis2=NULL, dbasis3=NULL,
                                      penalty, lambda, a,
                                      penalty_fusedcoef, lambda_fusedcoef,
                                      penalty_fusedbaseline, lambda_fusedbaseline,
                                      penweights_list, mu_smooth=0, mu_smooth_fused,
                                      step_size_init=1, step_size_min = 1e-6, step_size_max = 1e6,
                                      step_size_scale=1/2, ball_R=Inf, maxit=300,
                                      conv_crit = "omega", conv_tol=if(lambda>0) lambda/4 else 1e-6,
                                      verbose){


  #THIS IS STANDARD PROXIMAL GRADIENT DESCENT WITH A LINE SEARCH
  #THE BUILDING BLOCK OF THE PISTA ALGORITHM OF WANG ET AL. (2014) THAT I'M USING FOR A START
  #para is vector of starting values
  #y1-2 are first and second event times
  #delta1-2 are indicators of observation of first and second events
  #Xmat1-3 are design matrices corresponding to each transition (should be 'matrix' objects)
  #hazard is string saying form of hazard. "weibull"
  #frailty is boolean for whether to include a gamma frailty or not, currently everything is designed assuming a frailty
  #model is string saying "markov" or "semi-markov"
  #penalty is string for what form of coordinatewise penalty there should be
  #lambda is either a scalar value for the penalty parameter shared by all three transitions,
  #or a three-vector with the values for each of the three transitions
  #a is a scalar value for the second penalty parameter in SCAD, MCP, and adaptive lasso
  #Default is 3.7 for SCAD, 3 for MCP, and 1 for adaptive lasso
  #verbose is a boolean indicating whether to print the running fitting details
  #control is a list with options, outlined below
  #parahat is a vector of weights for adaptive lasso, typically arising as the MLE fit of the full model





  # browser()
  ##Set up control pieces##
  ##*********************##

  #construct control list
  #conv_crit determines criterion used to judge convergence
  #est_change_norm' looks at l1 norm of change in parameter values (scaled by l1 norm of prior values): sum(abs(finalVals-prevVals))/sum(abs(prevVals))
  #max_est_change' looks at largest absolute change in a parameter value
  #nll_pen_change' looks at change in regularized log-likelihood
  #maj_grad_norm' looks at l1 norm of majorized gradient
  if("maxit" < 0){
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
  nll_pen_trace <- rep(NA,maxit)
  trace_mat <- xcurr <- para
  fit_code <- 4 #follows convention from 'nleqslv' L-BFGS that is 1 if converged, 4 if maxit reached, 3 if stalled
  bad_step_count <- restart_count <- 0


  #note the division by n to put it on the mean scale--gradient descent works better then!
  nll_pen_xcurr <- nll_pen_func(para=xcurr,
                                y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                hazard=hazard, frailty=frailty, model=model,
                                basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                penalty=penalty,lambda=lambda, a=a,
                                penweights_list=penweights_list,
                                pen_mat_w_lambda = pen_mat_w_lambda,
                                mu_smooth_fused = mu_smooth_fused)/n
  if(is.na(nll_pen_xcurr)){stop("Initial values chosen yield infinite likelihood values. Consider other initial values.")}

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
    xnext <- prox_func(para=xcurr-ngrad_xcurr * step_size_ti, prev_para = xcurr,
                       nP1=nP1,nP2=nP2,nP3=nP3,step_size=step_size_ti,
                       penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                       pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig,
                       lambda_f_vec=lambda_f_vec, mu_smooth=mu_smooth,
                       mu_smooth_fused=mu_smooth_fused, ball_R=ball_R)

    #note the division by n to put it on the mean scale--gradient descent works better then!
    nll_pen_xnext <- nll_pen_func(para=xnext,
                                  y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                  hazard=hazard, frailty=frailty, model=model,
                                  basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                  penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                  pen_mat_w_lambda = pen_mat_w_lambda,
                                  mu_smooth_fused = mu_smooth_fused)/n

    smooth_obj_lqa_pen_xnext <- smooth_obj_lqa_pen_func(para=xnext, prev_para=xcurr,
                                                        y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                                        hazard=hazard, frailty=frailty, model=model,
                                                        basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                                        dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                                        penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                                        pen_mat_w_lambda = pen_mat_w_lambda,
                                                        mu_smooth=mu_smooth, mu_smooth_fused = mu_smooth_fused,
                                                        step_size=step_size_ti)/n

    if(is.nan(nll_pen_xnext)){
      if(verbose)print("whoops, proposed step yielded infinite likelihood value, just fyi")
      nll_pen_xnext <- smooth_obj_lqa_pen_xnext <- Inf
    }

    while(is.infinite(nll_pen_xnext) || nll_pen_xnext > smooth_obj_lqa_pen_xnext){
      if(step_size_ti < step_size_min){
        if(verbose)print("step size got too small, accepting result anyways")
        bad_step_count <- bad_step_count + 1
        break
      }
      step_size_ti <- step_size_ti * step_size_scale
      if(verbose)print(paste0("nll_pen_xnext:",nll_pen_xnext," bigger than smooth_obj_lqa_pen_xnext:",smooth_obj_lqa_pen_xnext))
      if(verbose)print(paste0("effective step size reduced to: ",step_size_ti))

      xnext <- prox_func(para=xcurr-ngrad_xcurr*step_size_ti, prev_para = xcurr,
                         nP1=nP1,nP2=nP2,nP3=nP3,step_size=step_size_ti,
                         penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                         pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig,
                         lambda_f_vec=lambda_f_vec, mu_smooth=mu_smooth,
                         mu_smooth_fused=mu_smooth_fused, ball_R=ball_R)

      #note the division by n to put it on the mean scale--gradient descent works better then!
      nll_pen_xnext <- nll_pen_func(para=xnext,
                                    y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                    hazard=hazard, frailty=frailty, model=model,
                                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                    penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                    pen_mat_w_lambda = pen_mat_w_lambda,
                                    mu_smooth_fused = mu_smooth_fused)/n

      smooth_obj_lqa_pen_xnext <- smooth_obj_lqa_pen_func(para=xnext, prev_para=xcurr,
                                                          y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                                          Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                                          hazard=hazard, frailty=frailty, model=model,
                                                          basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                                          dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                                          penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                                          mu_smooth=mu_smooth, pen_mat_w_lambda = pen_mat_w_lambda,
                                                          mu_smooth_fused = mu_smooth_fused,step_size=step_size_ti)/n

      if(is.nan(nll_pen_xnext)){
        if(verbose)print("whoops, proposed step yielded infinite likelihood value, just fyi")
        nll_pen_xnext <- smooth_obj_lqa_pen_xnext <- Inf
      }
    }

    ##UPDATE MONITORS##
    ##***************##

    xprev <- xcurr
    xcurr <- xnext

    # trace_mat <- cbind(trace_mat,xnext)

    # nll_pen_trace[i] <- nll_pen_xnext
    nll_pen_xnext <- nll_pen_trace[i] <- nll_pen_func(para=xnext,
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

    #convergence criterion given in Wang (2014)
    omega_t <- max(abs(prox_func(para=ngrad_xcurr, prev_para = xprev,
                                 nP1=nP1,nP2=nP2,nP3=nP3,step_size=1,
                                 penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                                 pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig,
                                 lambda_f_vec=lambda_f_vec, mu_smooth=mu_smooth,
                                 mu_smooth_fused = mu_smooth_fused, ball_R=ball_R)))

    max_change <- max(abs(xcurr-xprev))
    norm_change <- sqrt(sum((xcurr-xprev)^2))
    nll_pen_change <- nll_pen_xnext-nll_pen_xcurr

    if(verbose){
      print(paste("max change in ests",max_change))
      print(paste("omega_t (max norm of prox grad)", omega_t))
      print(paste("estimate with max change",names(para)[abs(xcurr-xprev) == max_change]))
      print(paste("max norm of change", norm_change)) #essentially a change in estimates norm
      print(paste("change in nll_pen", nll_pen_change))
      print(paste("new nll_pen", nll_pen_xnext))
    }

    if("omega" %in% conv_crit){
      if(omega_t <= conv_tol){break} #conv_tol is called 'epsilon' in the Wang (2014) paper
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
      if(abs(nll_pen_change) < conv_tol){
        fit_code <- 2
        break
      }
    }

    #after the updates, now update the nll monitor values
    nll_pen_xcurr <- nll_pen_xnext
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
                               ball_R=ball_R)
  ##Return list of final objects##
  ##****************************##

  return(list(estimate=finalVals, niter=i, fit_code=fit_code,
              final_nll=final_nll, final_nll_pen=final_nll_pen,
              final_ngrad=final_ngrad, final_ngrad_pen=final_ngrad_pen,
              startVals=para, hazard=hazard, frailty=frailty, model=model,
              penalty=penalty, lambda=lambda, a=a,
              penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
              penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
              mu_smooth_fused=mu_smooth_fused,
              # trace_mat=as.matrix(trace_mat),
              nll_pen_trace = n*nll_pen_trace,
              control=list(step_size_init=step_size_init,
                           step_size_min=step_size_min,
                           step_size_max=step_size_max,
                           step_size_scale=step_size_scale,
                           conv_crit=conv_crit,
                           conv_tol=conv_tol,
                           maxit=maxit),
              ball_R=ball_R,
              final_step_size=step_size_ti,
              bad_step_count=bad_step_count))
}
