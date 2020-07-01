#' Non-Monotone Accelerated Proximal Gradient Descent
#'
#' This function runs a non-monotone proximal gradient descent algorithm with Barzilai-Borwein
#'   backtracking line search, following the algorithm of Li & Lin (2015).
#'
#' @inheritParams proximal_gradient_descent
#' @param step_size_init_x,step_size_init_y Positive numeric values for the initial step sizes.
#' @param step_delta Positive numeric value for the parameter governing sufficient descent criterion.
#'
#' @return A list.
#' @export
proximal_gradient_descent_nmaccel <- function(para, y1, y2, delta1, delta2,
                                              Xmat1 = matrix(nrow(length(y1)),ncol=0), #the default is a 'empty' matrix
                                              Xmat2 = matrix(nrow(length(y1)),ncol=0),
                                              Xmat3 = matrix(nrow(length(y1)),ncol=0),
                                              hazard, frailty, model,
                                              basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                                              dbasis1=NULL, dbasis2=NULL, dbasis3=NULL,
                                              penalty, lambda, a,
                                              penalty_fusedcoef, lambda_fusedcoef,
                                              penalty_fusedbaseline, lambda_fusedbaseline,
                                              penweights_list, mu_smooth_fused,
                                              step_size_min=1e-6, step_size_max=1e6,
                                              step_size_init_x=1,step_size_init_y=1,
                                              step_size_scale=0.5, #no checks implemented on these values!!
                                              step_delta=0.5,ball_R=Inf, maxit=300,
                                              conv_crit = "nll_pen_change",
                                              conv_tol=1e-6,
                                              verbose){

  #THIS IS THE ALGORITHM OF LI & LIN (2015), WHICH IS NOT REALLY DESIGNED FOR A SOLUTION PATH BUT FOR FINDING A SPECIFIC SOLUTION.
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


  #START HERE AND UPDATE THE REST OF THE FUNCTION WITH THE NEW SYNTAX

  # browser()
  ##Set up control pieces##
  ##*********************##

  #construct control list
  #maxit is maximum number of newton-raphson iterations allowed
  #min_lr is the minimum learning rate allowed through step-halving process.
  #conv_crit determines criterion used to judge convergence
  #est_change_norm' looks at l1 norm of change in parameter values (scaled by l1 norm of prior values): sum(abs(finalVals-prevVals))/sum(abs(prevVals))
  #max_est_change' looks at largest absolute change in a parameter value
  #nll_pen_change' looks at change in regularized log-likelihood
  #maj_grad_norm' looks at l1 norm of majorized gradient
  # if (length(noNms <- namc[!namc %in% nmsC])){
  #   warning("unknown names in control: ", paste(noNms, collapse=", "))
  # }
  if(!all(conv_crit %in% c("est_change_norm","max_est_change","nll_pen_change","suboptimality"))){
    stop("unknown convergence criterion.")
  }
  if(maxit <0){
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
  pen_mat_w_eig <- pen_mat_decomp(pen_mat_w)

  ##Run algorithm##
  ##*************##
  nll_pen_trace <- rep(NA,maxit)
  trace_mat <- zcurr <- xnext <- xcurr <- xprev <- startVals <- para
  fit_code <- 4 #follows convention from 'nleqslv' L-BFGS that is 1 if converged, 4 if maxit reached, 3 if stalled
  bad_step_count <- restart_count <- 0

  #define things using notation from Li & Lin (2015) appendix with backtracking
  #outputs
  #x_{k} = xcurr #the actual estimates
  #x_{k-1} = xprev #the estimate at the previous step
  #intermediate outputs
  #y_{k} = ycurr #the extrapolated current point based on current nesterov momentum
  #z_{k} = zcurr #the proposed step based on gradient from y_{k-1}
  #z_{k+1} = znext #the proposed step based on gradient from y_{k}
  #v_{k+1} = vnext #the proposed step based on gradient from x_{k}

  #monitors
  #t_{k} = tcurr #value capturing the amount of nesterov momentum at the current step
  #t_{k-1} = tprev #value capturing the amount of nesterov momentum at the previous step
  #s_{k} = scurr #value capturing relative change in estimates between steps k-1 and k
  #r_{k} = rcurr #value capturing relative change in gradient between steps k-1 and k
  tcurr <- 1
  tprev <- 0
  #q_{k} = qcurr
  #q_{k+1} = qnext
  qnext <- qcurr <-1
  #c_{k} = ccurr

  #note the division by n to put it on the mean scale--gradient descent works better then!
  nll_pen_xnext <- nll_pen_xcurr <- ccurr <- nll_pen_func(para=xcurr,
                                                          y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                                          Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                                          hazard=hazard, frailty=frailty, model=model,
                                                          basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                                          dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                                          penalty=penalty,lambda=lambda, a=a,
                                                          penweights_list=penweights_list,
                                                          pen_mat_w_lambda = pen_mat_w_lambda, mu_smooth_fused = mu_smooth_fused)/n
  if(is.na(nll_pen_xcurr)){stop("Initial values chosen yield infinite likelihood values. Consider other initial values.")}
  #constants
  #alpha_y = step_alphay (step size for line search on y)
  #now called step_size_y
  #alpha_x = step_alphax (step size for line search on x)
  #now called step_size_x
  # step_alphay <- step_alphax <- 1
  step_size_x <- step_size_init_x
  step_size_y <- step_size_init_y


  #rho =  step_rho < 1
  # we call this step_size_scale
  # step_rho <- 1/2 #setting default step size change to step halving

  #delta = step_delta > 0
  #eta = step_eta \in [0,1)
  # step_delta <- 0.25 #parameter controlling 'sufficient descent' property, meaning
  #step has to not only reduce objective, but reduce by a certain amount
  #the smaller delta, the smaller the reduction must be
  step_eta <- 0.5 #this one doesn't feel like it should be able to be changed by input,
  #there is also a conflict with the use of eta in the Wang  (2014) paper and the Li & Lin (2015) paper so
  #I'd like to reduce confusion

  i <- 1
  while(i <= maxit){
    if(verbose)print(i)
    # if(i==15){browser()}
    ##RUN ACTUAL STEP##
    ##***************##

    ycurr <- xcurr + tprev/tcurr * (zcurr - xcurr) +
      ((tprev-1)/tcurr) * (xcurr - xprev)

    #note the division by n to put it on the mean scale--gradient descent works better then!
    ngrad_ycurr <- smooth_obj_grad_func(para=ycurr,
                                        y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                        hazard=hazard, frailty=frailty, model=model,
                                        basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                        dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                        penalty=penalty,lambda=lambda, a=a,
                                        penweights_list=penweights_list,
                                        pen_mat_w_lambda = pen_mat_w_lambda,
                                        mu_smooth_fused = mu_smooth_fused)/n

    #Barzilai-Borwein step size estimation
    if(i>1){
      scurr <- ycurr - yprev
      rcurr <- ngrad_ycurr - ngrad_yprev
      step_size_y <- min(sum(scurr^2)/sum(scurr*rcurr),step_size_max)
    }

    #run prox function:
    #ADMM prox subroutine for fused lasso only if mu_smooth_fused=0
    #soft thresholding according to lambda*step*weight
    znext <- prox_func(para=ycurr-step_size_y*ngrad_ycurr, prev_para = ycurr,
                       nP1=nP1,nP2=nP2,nP3=nP3,step_size=step_size_y,
                       penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                       pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig, ball_R=ball_R,
                       lambda_f_vec=lambda_f_vec, mu_smooth_fused=mu_smooth_fused)

    #note the division by n to put it on the mean scale--gradient descent works better then!
    nll_pen_znext <- nll_pen_func(para=znext,
                                  y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                  hazard=hazard, frailty=frailty, model=model,
                                  basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                  penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                  pen_mat_w_lambda = pen_mat_w_lambda,mu_smooth_fused = mu_smooth_fused)/n

    if(is.nan(nll_pen_znext)){
      if(verbose)print("whoops, proposed z step yielded infinite likelihood value, just fyi")
      nll_pen_znext <- Inf
    }


    #note the division by n to put it on the mean scale--gradient descent works better then!
    nll_pen_ycurr <- nll_pen_func(para=ycurr,
                                  y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                  basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                  hazard=hazard, frailty=frailty, model=model,
                                  penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                  pen_mat_w_lambda = pen_mat_w_lambda, mu_smooth_fused = mu_smooth_fused)/n

    #threshold for 'sufficient descent' (though I should monitor for whether this needs to be scaled by n)
    step_thresh1 <- nll_pen_ycurr - step_delta*sum((znext-ycurr)^2)
    #threshold for 'sufficient descent' (though I should monitor for whether this needs to be scaled by n)
    step_thresh2 <- ccurr - step_delta*sum((znext-ycurr)^2)

    while(nll_pen_znext > step_thresh1 & nll_pen_znext > step_thresh2){
      if(abs(step_size_y) < step_size_min){
        if(verbose)print("y step size got too small, accepting result anyways")
        bad_step_count <- bad_step_count + 1
        break
      }

      step_size_y <- step_size_scale * step_size_y

      znext <- prox_func(para=ycurr-step_size_y*ngrad_ycurr, prev_para=ycurr,
                         nP1=nP1,nP2=nP2,nP3=nP3,step_size=step_size_y,
                         penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                         pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig, ball_R=ball_R,
                         lambda_f_vec=lambda_f_vec, mu_smooth_fused=mu_smooth_fused)

      #note the division by n to put it on the mean scale--gradient descent works better then!
      nll_pen_znext <- nll_pen_func(para=znext,
                                    y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                    hazard=hazard, frailty=frailty, model=model,
                                    penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                    pen_mat_w_lambda = pen_mat_w_lambda,mu_smooth_fused = mu_smooth_fused)/n

      if(is.nan(nll_pen_znext)){
        if(verbose)print("whoops, proposed z step yielded infinite likelihood value, just fyi")
        nll_pen_znext <- Inf
      }

      #threshold for 'sufficient descent' (though I should monitor for whether this needs to be scaled by n)
      step_thresh1 <- nll_pen_ycurr - step_delta*mean((znext-ycurr)^2)
      #threshold for 'sufficient descent' (though I should monitor for whether this needs to be scaled by n)
      step_thresh2 <- ccurr - step_delta*mean((znext-ycurr)^2)

      if(verbose)print(paste0("y step size reduced to: ",step_size_y))

    }

    if(nll_pen_znext <= step_thresh2){
      xnext <- znext
    } else {

      #note the division by n to put it on the mean scale--gradient descent works better then!
      ngrad_xcurr <- smooth_obj_grad_func(para=xcurr,
                                          y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                          Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                          basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                          dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                          hazard=hazard, frailty=frailty, model=model,
                                          penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                          pen_mat_w_lambda = pen_mat_w_lambda,mu_smooth_fused = mu_smooth_fused)/n

      if(i>1){
        scurr <- xcurr - yprev
        rcurr <- ngrad_xcurr - ngrad_yprev
        step_size_x <- min(sum(scurr^2)/sum(scurr*rcurr),step_size_max)
      }

      vnext <- prox_func(para=xcurr-step_size_x*ngrad_xcurr, prev_para=xcurr,
                         nP1=nP1,nP2=nP2,nP3=nP3,step_size=step_size_x,
                         penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                         pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig, ball_R=ball_R,
                         lambda_f_vec=lambda_f_vec, mu_smooth_fused=mu_smooth_fused)

      #note the division by n to put it on the mean scale--gradient descent works better then!
      nll_pen_vnext <- nll_pen_func(para=vnext,
                                    y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                    hazard=hazard, frailty=frailty, model=model,
                                    penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                    pen_mat_w_lambda = pen_mat_w_lambda,mu_smooth_fused = mu_smooth_fused)/n

      #threshold for 'sufficient descent' (though I should monitor for whether this needs to be scaled by n)
      step_thresh3 <- ccurr - step_delta*mean((vnext-xcurr)^2)

      if(is.nan(nll_pen_vnext)){
        if(verbose)print("whoops, proposed v step yielded infinite likelihood value, just fyi")
        nll_pen_vnext=Inf
      }

      while(nll_pen_vnext >= step_thresh3){
        if(abs(step_size_x) < step_size_min){
          if(verbose)print("x step size got too small, accepting result anyways")
          bad_step_count <- bad_step_count + 1
          break
        }



        step_size_x <- step_size_scale * step_size_x

        vnext <- prox_func(para=xcurr-step_size_x*ngrad_xcurr, prev_para = xcurr,
                           nP1=nP1,nP2=nP2,nP3=nP3,step_size=step_size_x,
                           penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                           pen_mat_w=pen_mat_w, pen_mat_w_eig=pen_mat_w_eig, ball_R=ball_R,
                           lambda_f_vec=lambda_f_vec, mu_smooth_fused=mu_smooth_fused)

        #note the division by n to put it on the mean scale--gradient descent works better then!
        nll_pen_vnext <- nll_pen_func(para=vnext,
                                      y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                      basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                      dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                      hazard=hazard, frailty=frailty, model=model,
                                      penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                      pen_mat_w_lambda = pen_mat_w_lambda, mu_smooth_fused = mu_smooth_fused)/n

        if(is.nan(nll_pen_vnext)){
          if(verbose)print("whoops, proposed v step yielded infinite likelihood value, just fyi")
          nll_pen_vnext=Inf
        }


        #threshold for 'sufficient descent' (though I should monitor for whether this needs to be scaled by n)
        step_thresh3 <- ccurr - step_delta*sum((vnext-xcurr)^2)
        if(verbose)print(paste0("x step size reduced to: ",step_size_x))
      }

      if(nll_pen_znext < nll_pen_vnext){
        xnext <- znext
      } else{
        xnext <- vnext
      }

    }

    ##UPDATE MONITORS##
    ##***************##

    nll_pen_xcurr <- nll_pen_xnext
    #note the division by n to put it on the mean scale--gradient descent works better then!
    nll_pen_xnext <- nll_pen_func(para=xnext,
                                  y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                  basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                  hazard=hazard, frailty=frailty, model=model,
                                  penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                  pen_mat_w_lambda = pen_mat_w_lambda, mu_smooth_fused = mu_smooth_fused)/n

    tprev <- tcurr
    tcurr <- (1/2)*(sqrt(4*tprev^2 + 1) + 1)
    qcurr <- qnext
    qnext <- step_eta*qcurr + 1
    ccurr <- (step_eta*qcurr*ccurr + nll_pen_xnext) / qnext #updated directly rather than via a cnext, because it's the last step!

    yprev <- ycurr
    ngrad_yprev <- ngrad_ycurr

    xprev <- xcurr
    xcurr <- xnext
    # trace_mat <- cbind(trace_mat,xnext)

    #RECORD RESULT WITHOUT SMOOTHING, TO PUT US ON A COMMON SCALE!
    nll_pen_trace[i] <- nll_pen_func(para=xnext,
                                     y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                     Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                     basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                     dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                     hazard=hazard, frailty=frailty, model=model,
                                     penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                     pen_mat_w_lambda = pen_mat_w_lambda,mu_smooth_fused = 0)/n

    ##Check for convergence##
    ##*********************##

    max_est_change <- max(abs(xcurr-xprev))
    norm_est_change <- sqrt(sum((xcurr-xprev)^2))
    nll_pen_change <- nll_pen_xnext-nll_pen_xcurr



    if(verbose){
      print(paste("max change in ests",max_est_change))
      print(paste("estimate with max change",names(para)[abs(xcurr-xprev) == max_est_change]))
      # print(paste("suboptimality (max norm of prox grad)", subopt_t))
      print(paste("l2 norm of change", norm_est_change)) #essentially a change in estimates norm
      print(paste("change in nll_pen", nll_pen_change))
      print(paste("new nll_pen", nll_pen_xnext))
    }





    if("suboptimality" %in% conv_crit){
      #convergence criterion given in Wang (2014)
      #note the division by n to put it on the mean scale--gradient descent works better then!
      ngrad_xcurr <- smooth_obj_grad_func(para=xcurr,
                                          y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                          Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                          hazard=hazard, frailty=frailty, model=model,
                                          basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                          dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                          penalty=penalty,lambda=lambda, a=a,
                                          penweights_list=penweights_list,
                                          pen_mat_w_lambda = pen_mat_w_lambda,mu_smooth_fused = mu_smooth_fused)/n

      #suboptimality convergence criterion given as omega in Wang (2014)
      subopt_t <- max(abs(prox_func(para=ngrad_xcurr, prev_para = xprev,
                                   nP1=nP1,nP2=nP2,nP3=nP3,step_size=1,
                                   penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                                   pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig,
                                   lambda_f_vec=lambda_f_vec,
                                   mu_smooth_fused = mu_smooth_fused, ball_R=ball_R)))
      if(verbose){print(paste("suboptimality (max norm of prox grad)", subopt_t))}
      if(subopt_t < conv_tol){
        fit_code <- 2
        break
      }
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

  #convergence criterion given in Wang (2014)
  #note the division by n to put it on the mean scale--gradient descent works better then!
  ngrad_xcurr <- smooth_obj_grad_func(para=xcurr,
                                      y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                      hazard=hazard, frailty=frailty, model=model,
                                      basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                      dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                      penalty=penalty,lambda=lambda, a=a,
                                      penweights_list=penweights_list,
                                      pen_mat_w_lambda = pen_mat_w_lambda,mu_smooth_fused = mu_smooth_fused)/n

  #suboptimality convergence criterion given as omega in Wang (2014)
  subopt_t <- max(abs(prox_func(para=ngrad_xcurr, prev_para = xprev,
                                nP1=nP1,nP2=nP2,nP3=nP3,step_size=1,
                                penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                                pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig,
                                lambda_f_vec=lambda_f_vec,
                                mu_smooth_fused = mu_smooth_fused, ball_R=ball_R)))


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
                                      pen_mat_w_lambda = pen_mat_w_lambda,mu_smooth_fused = mu_smooth_fused) #here, we do still report the nesterov-smoothed results


  final_ngrad_pen <- prox_func(para=final_ngrad, prev_para = xprev,
                               nP1=nP1,nP2=nP2,nP3=nP3,step_size=1,
                               penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                               pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig,
                               lambda_f_vec=lambda_f_vec, mu_smooth_fused = mu_smooth_fused,
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
              control=list(step_size_init_x=step_size_init_x,
                           step_size_init_y=step_size_init_y,
                           step_size_min=step_size_min,
                           step_size_max=step_size_max,
                           step_size_scale=step_size_scale,
                           conv_crit=conv_crit,conv_tol=conv_tol,
                           maxit=maxit),
              conv_stats=c(max_est_change=max_est_change,norm_est_change=norm_est_change,nll_pen_change=nll_pen_change,subopt=subopt_t),
              ball_R=ball_R,
              final_step_size=step_size_x,
              final_step_size_x=step_size_x,
              final_step_size_y=step_size_y,
              bad_step_count=bad_step_count))
}
