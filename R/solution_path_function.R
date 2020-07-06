#' Estimate Penalized Illness-Death Model Solution Path
#'
#' This function estimates penalized illness-death model results along a range of
#'   penalty, fused penalty, and smoothing parameters.
#'
#' This is a function to loop through a path of lambda, lambda_fusedcoef, and mu_smooth values in a somewhat
#'   thoughtful way to maximize the pathwise connections between starting values and step sizes,
#'   with some adjustments tailored to each approach to optimization.
#' @inheritParams proximal_gradient_descent
#' @inheritParams proximal_gradient_descent_nmaccel
#' @inheritParams newton_raphson_mm
#' @param lambda_path Numeric sequence of decreasing regularization parameters
#'   for the parameterwise penalties, along which the solution path runs.
#'   Assumes a single shared penalty across transitions.
#' @param lambda_fusedcoef_path Numeric sequence of increasing regularization parameters
#'   for the fusion penalties.
#' @param mu_smooth_path Numeric sequence of decreasing Nesterov smoothing parameters for the fusion penalties.
#' @param fit_method String indicating which optimization method should be used at each step.
#' @param warm_start Boolean indicating whether each step of the solution should start from
#'   ending of previous (\code{TRUE}) or from the original starting point (\code{FALSE}).
#' @param extra_starts numeric indicating how many additional optimization runs from random start values
#'   should be performed at each grid point.
#' @param fusion_tol Numeric value indicating when to consider fused parameters
#'   that are close to be considered the same value, for estimating degrees of freedom.
#'
#' @return A list.
#' @export
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
                                   fit_method, warm_start=TRUE, extra_starts=0,
                                   select_tol=1e-4, fusion_tol = 1e-3,
                                   step_size_min=1e-6, step_size_max=1e6,
                                   step_size_init=1,
                                   step_size_scale=0.5, #no checks implemented on these values!!
                                   step_delta=0.5,
                                   maxit=300,
                                   conv_crit = "nll_pen_change",
                                   conv_tol=1e-6,
                                   mm_epsilon=1e-6,
                                   verbose){

  # browser()

  #If lambda_path is NULL, set it equal to 0 (aka no parameterwise regularization)
  #Otherwise, put it into matrix format
  # (note, this admits a 3-column matrix to let different transitions have different values)
  lambda_path <- if(is.null(lambda_path)) as.matrix(0) else as.matrix(lambda_path)
  lambda_length <- nrow(lambda_path)
  colnames(lambda_path) <- paste0("lambda",1:ncol(lambda_path))

  #If lambda_fusedcoef_path is NULL, set it equal to 0 (aka no fused regularization)
  #Otherwise, put it into matrix format
  # (note, this admits a 3-column matrix to let different transitions have different values)
  lambda_fusedcoef_path <- if(is.null(lambda_fusedcoef_path)) as.matrix(0) else as.matrix(lambda_fusedcoef_path)
  lambda_fusedcoef_length <- nrow(lambda_fusedcoef_path)
  colnames(lambda_fusedcoef_path) <- paste0("lambda_fusedcoef",1:ncol(lambda_fusedcoef_path))

  #If mu_smooth_path is NULL, set it equal to 0 (aka no fused smoothing)
  #Otherwise, put it into matrix format
  mu_smooth_path <- if(is.null(mu_smooth_path)) as.matrix(0) else as.matrix(mu_smooth_path)
  mu_smooth_length <- nrow(mu_smooth_path)
  colnames(mu_smooth_path) <- "mu_smooth" #for now, implicit that there is only one column


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
  # total_length <- lambda_length*sum(lambda_fusedcoef_path %in% 0) + #total number of steps with no fusion
  #   lambda_length*sum(!(lambda_fusedcoef_path %in% 0)) * mu_smooth_length #total number of fusion steps

  #in this new version, we're not storing every single step of the mu_smoothing loop, just the final one
  #therefore, we don't need to play around with how many iterations are smoothed vs not smoothed.
  #we just get a single output, regardless of whether there is any fusion or not.
  total_length <- lambda_length*lambda_fusedcoef_length

  out_starts <- out_pen_ngrads <- out_ests <- matrix(nrow=total_length,
                                                     ncol=nPtot,dimnames=list(NULL,names(para)))
  out_ics <- matrix(nrow=total_length, ncol=15,
                    dimnames=list(NULL,c("nll", "nll_pen",
                                         "nPtot", "nPtot_selected",
                                         "df_unique", "df_est",
                                         "AIC", "AIC_unique", "AIC_estdf",
                                         "BIC", "BIC_unique", "BIC_estdf",
                                         "GCV", "GCV_unique", "GCV_estdf")))
  out_info <- matrix(nrow=total_length,
                     ncol=ncol(lambda_path) + ncol(lambda_fusedcoef_path) + ncol(mu_smooth_path) + 1,
                     dimnames=list(NULL,c(colnames(lambda_path),
                                          colnames(lambda_fusedcoef_path),
                                          colnames(mu_smooth_path),"ball_R")))
  out_control <- matrix(nrow=total_length,
                        ncol=4,dimnames=list(NULL,c("niter","fit_code","final_step_size","best_start_iter")))

  out_conv_stats <- matrix(nrow=total_length,
                           ncol=4,dimnames=list(NULL,c("est_change_max","est_change_2norm",
                                                       "nll_pen_change","subopt")))

  nll_pen_trace_mat <-  matrix(nrow=total_length,
                               ncol=maxit)

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

    #Wang (2014) paper advises specific weaker tolerance earlier in the path, according to lambda value:
    # conv_tol <- if(lambda_iter==lambda_length) lambda/8 else lambda/4 #gotta be less than lambda_target/4, so we just halve it again.

    startVals_middle <- startVals_outer
    step_size_init_middle <- step_size_init_outer

    #MIDDLE LOOP: LOOPING THROUGH THE LAMBDA_FUSEDCOEF VALUES, GOVERNING FUSION SPARSITY
    for(lambda_fusedcoef_iter in 1:lambda_fusedcoef_length){

      lambda_fusedcoef <- lambda_fusedcoef_path[lambda_fusedcoef_iter,]
      startVals_inner <- startVals_middle
      step_size_init_inner <- step_size_init_middle

      if(verbose >= 1){
        print(paste0("step: ",lambda_iter,
                     " lambda: ",lambda,
                     " lambda_fusedcoef: ",lambda_fusedcoef))#, " mu_smooth_fused: ",mu_smooth_fused, " extra_start: ", extra_start_iter))
      }

      #monitor to track which random start achieves the lowest objective value
      best_nll_pen <- Inf

      #INNER LOOP: LOOPING THROUGH EXTRA START VALUES
      #Unfortunately this sort of disrupts the point of having multiple mu_smooth values, but because we're
      #moving away from that approach anyways to a single mu_smooth value, then this is fine.
      for(extra_start_iter in 1:(extra_starts+1)){

        #if this is the first time through, admit the warm start approach, otherwise randomize the start value
        #and reset the step size to the global default
        if(extra_start_iter == 1){
          startVals_innermost <- startVals_inner
          step_size_init_innermost <- step_size_init_inner
        } else{
          #Here I'm basing my perturbed start values based off of the 'middle' starting values
          startVals_innermost <- (startVals_inner + stats::rnorm(n = length(startVals_inner),mean = 0,sd=0.7)) *
            stats::runif(n = length(startVals_inner),min = 0.9, max = 1.1) #adding random multiplicative and additive noise
          step_size_init_innermost <- step_size_init
        }

        #INNERMOST LOOP: LOOPING THROUGH MU_SMOOTH VALUES, GOVERNING SMOOTH APPROXIMATION OF FUSION SPARSITY
        for(mu_smooth_iter in 1:mu_smooth_length){

          #if there is no fusion, there should be no smoothing
          mu_smooth_fused <- if(lambda_fusedcoef==0) 0 else mu_smooth_path[mu_smooth_iter,]

          if(fit_method=="prox_grad"){
            tempfit <- proximal_gradient_descent(para=startVals_innermost, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                                 hazard=hazard, frailty=frailty, model=model,
                                                 basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                                 dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                                 penalty=penalty, lambda=lambda, a=a,
                                                 penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                                 penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                                 penweights_list=penweights_list, mu_smooth_fused=mu_smooth_fused,
                                                 step_size_init=step_size_init_innermost,
                                                 step_size_min=step_size_min,step_size_max=step_size_max,
                                                 step_size_scale=step_size_scale,
                                                 maxit=maxit, conv_crit=conv_crit, conv_tol=conv_tol,
                                                 ball_R=Inf, verbose=verbose)
          }
          if(fit_method=="prox_grad_nmaccel"){
            tempfit <- proximal_gradient_descent_nmaccel(para=startVals_innermost, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                                         Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                                         hazard=hazard, frailty=frailty, model=model,
                                                         basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                                         dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                                         penalty=penalty, lambda=lambda, a=a,
                                                         penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                                         penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                                         penweights_list=penweights_list, mu_smooth_fused=mu_smooth_fused,
                                                         step_size_init_x=step_size_init_innermost,
                                                         step_size_init_y=step_size_init_innermost,
                                                         step_size_min=step_size_min,step_size_max=step_size_max,
                                                         step_size_scale=step_size_scale,
                                                         maxit=maxit, conv_crit=conv_crit, conv_tol=conv_tol,
                                                         ball_R=Inf, verbose=verbose)
          }
          if(fit_method=="newton"){
            tempfit <- newton_raphson_mm(para=startVals_innermost, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                         Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                         hazard=hazard, frailty=frailty, model=model,
                                         penalty=penalty, lambda=lambda, a=a,
                                         penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                         penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                         penweights_list=penweights_list,
                                         mm_epsilon=mm_epsilon,
                                         verbose=verbose, maxit=maxit,step_size_min=step_size_min,
                                         conv_crit=conv_crit,conv_tol=conv_tol,num_restarts=2)
          }

          #now, in this innermost loop we really just step through the iterations,
          #without really tracking the intermediate smoothing results
          startVals_innermost <- tempfit$estimate

          # Wang (2014) suggest in their pathwise approach to carry over the step size from each iteration
          # So we could do that here, but it's currently commented out, and even then
          # only written for the proximal algorithm.
          # The accelerated method re-estimates its step size
          # at each step anyways so the initial step being too big would be less of an issue.
          # if(fit_method %in%  c("prox_grad")){
          #   step_size_init_innermost <- tempfit$final_step_size
          # }

          #if there is no smoothing, then break out of the smoothing (innermost) loop
          #recall that above we set this to 0 if theres no fusion at all.
          if(mu_smooth_fused==0){
            break
          }
        } #END OF INNERMOST (SMOOTHING) LOOP

        #Now, having run the innermost loop though all of the mu_smooth steps (if there's only one, just the one)
        #We assess whether it has reached a better spot than the previous start

        #If this random start has reached a better penalized value  (this is a completely unsmoothed value, fyi)
        #than the previous start, then update all of the resulting outputs
        if(tempfit$final_nll_pen < best_nll_pen){
          if(verbose >= 1 && extra_start_iter > 1){
            print(paste0("start ",extra_start_iter,
                         "yielded better obj value ", tempfit$final_nll_pen,
                         " over ",best_nll_pen))
          }
          best_nll_pen <- tempfit$final_nll_pen

          ##Now, let's COMPUTE SOME USEFUL FIT STATISTICS based on this fit##
          ##***************************************************************##
          final_nll <- tempfit$final_nll

          beta1_selected <- if(nP1 != 0) tempfit$estimate[(1+nP0):(nP0+nP1)] else numeric(0)
          nP1_selected <- sum(beta1_selected != 0)
          beta2_selected <- if(nP2 != 0) tempfit$estimate[(1+nP0+nP1):(nP0+nP1+nP2)] else numeric(0)
          nP2_selected <- sum(beta2_selected != 0)
          beta3_selected <- if(nP3 != 0) tempfit$estimate[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)] else numeric(0)
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
            df_unique <- nP0 + sum(beta1_selected != 0 &
                                     abs(beta1_selected - beta2_selected) > fusion_tol &
                                     abs(beta1_selected - beta3_selected) > fusion_tol) +
              sum(beta2_selected != 0 &
                    abs(beta2_selected - beta3_selected) > fusion_tol) +
              sum(beta3_selected != 0)

            #currently, these account for fusion in the 'degrees of freedom' by counting unique parameters
            tempaic_unique <- 2*final_nll + 2* df_unique
            tempbic_unique <- 2*final_nll + log(n) * df_unique
            tempgcv_unique <- final_nll/(n^2*(1-df_unique/n)^2)
          } else{
            df_unique <- tempaic_unique <- tempbic_unique <- tempgcv_unique <- NA
          }

          #alternative definition of degrees of freedom from Sennhenn-Reulen & Kneib via Gray (1993)
          # browser()
          if(!is.null(tempfit$final_nhess)){
            final_nhess_nopen <- nhess_func(para=tempfit$estimate, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                            Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                            hazard=hazard, frailty=frailty, model=model)
            final_cov <- tryCatch(solve(tempfit$final_nhess),
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
                              nPtot, nPtot_selected, df_unique, df_est,
                              tempaic, tempaic_unique,tempaic_estdf,
                              tempbic, tempbic_unique,tempbic_estdf,
                              tempgcv, tempgcv_unique,tempgcv_estdf)
          out_control[iter,] <- c(tempfit$niter,tempfit$fit_code,tempfit$final_step_size,extra_start_iter)
          out_starts[iter,] <- startVals_innermost
          nll_pen_trace_mat[iter, ] <- tempfit$nll_pen_trace

          if(fit_method != "newton"){ #haven't implemented convergence stats output for newton yet
            out_conv_stats[iter,] <- tempfit$conv_stats
          }


          #If this new best value is arising during the first time through the fused coef (middle) loop,
          #then the estimates are relevant for the next iteration of the outer loop
          if(warm_start & lambda_fusedcoef_iter==1){
            startVals_outer <- tempfit$estimate
            #again, Wang (2014) suggests tempering the step size along the path but I'm setting that aside for now.
            # if(fit_method %in% c("prox_grad")){
            #   step_size_init_outer <- tempfit$final_step_size
            # }
          }

          #Moreover, if this new best value is arising during the first time through the fused coef (middle) loop,
          #then the step size at the end is also possibly relevant for the next iteration of the middle loop
          #though again, we're gonna set that aside for now
          # if(warm_start & mu_smooth_iter==1){
          #   if(fit_method %in% c("prox_grad")){
          #     step_size_init_middle <- tempfit$final_step_size
          #   }
          # }

        } #END OF IF STATEMENT HOUSING UPDATES IF NEW BEST VALUE IS REACHED

      } #END OF INNER (EXTRA STARTS) LOOP

      iter <- iter + 1

    } #END OF MIDDLE (FUSED COEFFICIENT LAMBDA) LOOP

    #now that we've reached the end of the inner loop, let's use the final best values
    #over all of the random starts as starting values for the next middle loop step
    #note we have already above determined the correct starting step size for the next iteration of the middle loop
    if(warm_start){
      startVals_middle <- startVals_inner
    }

  } #END OF OUTER (PARAMETERWISE LAMBDA) LOOP

  # browser()

  ##MAKE SOME SUMMARY PLOTS AND FINALIZE THE OUTPUT TABLES##
  ##******************************************************##

  rownames(out_info) <- NULL
  rownames(out_ics) <- NULL
  rownames(out_ests) <- NULL
  rownames(out_pen_ngrads) <- NULL
  rownames(out_starts) <- NULL
  rownames(out_control) <- NULL
  rownames(out_conv_stats) <- NULL
  rownames(nll_pen_trace_mat) <- NULL

  return(list(info=out_info,
              ests=out_ests,
              pen_ngrads=out_pen_ngrads,
              ics=out_ics,
              starts=out_starts,
              pen_trace_mat=nll_pen_trace_mat,
              # step_eta=step_eta,
              # plot_out=plot_out,
              # plot_out_ests=plot_out_ests,
              control=out_control,
              conv_stats=out_conv_stats,
              hazard=hazard, frailty=frailty, model=model,
              penalty=penalty, a=a, mm_epsilon=mm_epsilon,
              select_tol=select_tol, fusion_tol=fusion_tol,
              penalty_fusedcoef=penalty_fusedcoef,
              penalty_fusedbaseline=penalty_fusedbaseline,
              lambda_fusedbaseline=lambda_fusedbaseline,
              fit_method=fit_method))

}
