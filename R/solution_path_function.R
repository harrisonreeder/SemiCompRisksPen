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
          final_nhess_nopen <- nhess_func(para=startVals_inner, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
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
