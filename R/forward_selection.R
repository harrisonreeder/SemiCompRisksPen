#' Forward Selection Procedure for Parametric Illness-Death Model
#'
#' This method is meant as a comparator for the penalized estimation procedure.
#'
#' @inheritParams FreqID_HReg2
#' @inheritParams proximal_gradient_descent
#' @param vars character string with all of the possible variable names to be searched through.
#' @param select_crit a string indicating what criterion should be used to determine whether a covariate should be added.
#' @param fixed1 a string vector for parameters you want guaranteed included in arm 1
#' @param fixed2 a string vector for parameters you want guaranteed included in arm 2
#' @param fixed3 a string vector for parameters you want guaranteed included in arm 3
#'
#' @return A list.
#' @export
forward_selection <- function(vars, data, na.action="na.fail", subset=NULL,
                              hazard=c("weibull"),frailty=TRUE, model, knots_list = NULL,
                              fixed1=character(0),fixed2=character(0),fixed3=character(0),
                              optim_method="L-BFGS",
                              select_crit="bic", verbose=0, control=NULL){

  n <- nrow(data)

  #empty vectors/lists used to track the strings that have been 'used' aka are active in the model.
  used1 <- used2 <- used3 <- character(0)
  used1_list <- used2_list <- used3_list <- list()
  used1_vec <- used2_vec <- used3_vec <- character(0)

  #monitors of the search process.
  modelcount <- 1
  best_crit_round <- Inf
  continue_flag <- TRUE
  model_crit_vec <- numeric(0)
  counter <- 1

  #loop that iterates the forward selection process
  while(continue_flag){
    best_crit_round <- Inf
    if(verbose>=1){
      print(counter)
    }

    ##loop that iterates through the covariates
    for(currvar in unique(vars)) {
      if(verbose >=3){
        print(currvar)
      }

      run1 <- run2 <- run3 <- FALSE

      #During earlier project, had more complicated logic here to accomodate
      #interactions and squared terms, but for now things are simplified
      if(!(currvar %in% c(fixed1,used1))){
        run1 <- TRUE
      }
      if(!(currvar %in% c(fixed2,used2))){
        run2 <- TRUE
      }
      if(!(currvar %in% c(fixed3,used3))){
        run3 <- TRUE
      }

      # print(paste0("run1: ",run1,", run2: ",run2,", run3: ",run3))
      if(run1){
        form <- Formula::as.Formula(paste0("y1 + delta1 | y2 + delta2 ~", #note hardcoded names for now.
                                           paste(c(1,fixed1,used1,currvar),collapse = "+")," |",
                                           paste(c(1,fixed2,used2),collapse = "+")," |",
                                           paste(c(1,fixed3,used3),collapse = "+")," "))
        fit_temp <- FreqID_HReg2(Formula = form,
                                         data, na.action=na.action, subset=subset,
                                         hazard=hazard,frailty=frailty, hessian=FALSE,
                                         model=model, knots_list=knots_list,
                                         optim_method = optim_method)
        #rewrite to allow different criteria. Maybe make a specific function that computes different ones?
        model_crit_vec[modelcount] <- get_ic(nll = -fit_temp$logLike,
                                             df = sum(fit_temp$estimate != 0),
                                             n = n,
                                             ic = select_crit)
        if(model_crit_vec[modelcount] < best_crit_round) {best_crit_round <- model_crit_vec[modelcount]}
        used1_list[[modelcount]] <- c(used1,currvar)
        used2_list[[modelcount]] <- c(used2)
        used3_list[[modelcount]] <- c(used3)
        used1_vec[modelcount] <- paste0(c(used1,currvar),collapse=", ")
        used2_vec[modelcount] <- paste0(c(used2),collapse=", ")
        used3_vec[modelcount] <- paste0(c(used3),collapse=", ")

        if(verbose>= 4){
          print(paste0("used1: ",used1_vec[modelcount]))
          print(paste0("used2: ",used2_vec[modelcount]))
          print(paste0("used3: ",used3_vec[modelcount]))
        }
        modelcount = modelcount + 1
      }

      if(run2){
        form <- Formula::as.Formula(paste0("y1 + delta1 | y2 + delta2 ~",
                                  paste(c(1,fixed1,used1),collapse = "+")," |",
                                  paste(c(1,fixed2,used2,currvar),collapse = "+")," |",
                                  paste(c(1,fixed3,used3),collapse = "+")," "))
        fit_temp <- FreqID_HReg2(Formula = form,
                                         data, na.action=na.action, subset=subset,
                                         hazard=hazard,frailty=frailty, hessian=FALSE,
                                         model=model, knots_list=knots_list,
                                         optim_method = optim_method)
        #rewrite to allow different criteria. Maybe make a specific function that computes different ones?
        model_crit_vec[modelcount] <- get_ic(nll = -fit_temp$logLike,
                                             df = sum(fit_temp$estimate != 0),
                                             n = n,
                                             ic = select_crit)
        if(model_crit_vec[modelcount] < best_crit_round) {best_crit_round <- model_crit_vec[modelcount]}
        used1_list[[modelcount]] <- c(used1)
        used2_list[[modelcount]] <- c(used2,currvar)
        used3_list[[modelcount]] <- c(used3)
        used1_vec[modelcount] <- paste0(c(used1),collapse=", ")
        used2_vec[modelcount] <- paste0(c(used2,currvar),collapse=", ")
        used3_vec[modelcount] <- paste0(c(used3),collapse=", ")

        if(verbose>= 4){
          print(paste0("used1: ",used1_vec[modelcount]))
          print(paste0("used2: ",used2_vec[modelcount]))
          print(paste0("used3: ",used3_vec[modelcount]))
        }
        modelcount = modelcount + 1
      }

      if(run3){
        form <- Formula::as.Formula(paste0("y1 + delta1 | y2 + delta2 ~",
                                  paste(c(1,fixed1,used1),collapse = "+")," |",
                                  paste(c(1,fixed2,used2),collapse = "+")," |",
                                  paste(c(1,fixed3,used3,currvar),collapse = "+")," "))
        fit_temp <- FreqID_HReg2(Formula = form,
                                         data, na.action=na.action, subset=subset,
                                         hazard=hazard,frailty=frailty, hessian=FALSE,
                                         model=model, knots_list=knots_list,
                                         optim_method = optim_method)
        model_crit_vec[modelcount] <- get_ic(nll = -fit_temp$logLike,
                                             df = sum(fit_temp$estimate != 0),
                                             n = n,
                                             ic = select_crit)
        if(model_crit_vec[modelcount] < best_crit_round) {best_crit_round <- model_crit_vec[modelcount]}
        used1_list[[modelcount]] <- c(used1)
        used2_list[[modelcount]] <- c(used2)
        used3_list[[modelcount]] <- c(used3,currvar)
        used1_vec[modelcount] <- paste0(c(used1),collapse=", ")
        used2_vec[modelcount] <- paste0(c(used2),collapse=", ")
        used3_vec[modelcount] <- paste0(c(used3,currvar),collapse=", ")

        if(verbose>= 4){
          print(paste0("used1: ",used1_vec[modelcount]))
          print(paste0("used2: ",used2_vec[modelcount]))
          print(paste0("used3: ",used3_vec[modelcount]))
        }
        modelcount = modelcount + 1
      }
    }


    #tail used because if there's a tie for some reason, just pick the last
    best_model <- tail(which(min(model_crit_vec) == model_crit_vec),n=1)
    best_crit <- model_crit_vec[best_model]
    used1 <- used1_list[[best_model]]
    used2 <- used2_list[[best_model]]
    used3 <- used3_list[[best_model]]

    if(verbose>= 2){
      print(paste("best model at step", counter,"is:"))
      print(paste0("h1: "))
      print(c(fixed1,used1))
      print(paste0("h2: "))
      print(c(fixed2,used2))
      print(paste0("h3: "))
      print(c(fixed3,used3))
      print(paste0("with criterion value: ",best_crit))
    }

    if(best_crit < best_crit_round){
      continue_flag <- FALSE
    }
    counter <- counter + 1

  }

  crit_tab <- data.frame(model_crit_vec,used1_vec,used2_vec,used3_vec)

  #FIT BEST MODEL
  form <- Formula::as.Formula(paste0("y1 + delta1 | y2 + delta2 ~",
                                     paste(c(1,used1,fixed1),collapse = "+")," |",
                                     paste(c(1,used2,fixed2),collapse = "+")," |",
                                     paste(c(1,used3,fixed3),collapse = "+")," "))
  fit_temp <- FreqID_HReg2(Formula = form,
                            data, na.action=na.action, subset=subset,
                            hazard=hazard,frailty=frailty,
                            model=model, knots_list=knots_list,
                            optim_method = optim_method)

  return(list(final_fit=fit_temp,
              crit_tab=crit_tab,
              fixed1=fixed1,fixed2=fixed2,fixed3=fixed3,
              used1=used1,used2=used2,used3=used3,
              best_crit=best_crit,
              used1_list=used1_list,used2_list=used2_list,used3_list=used3_list,
              model_crit_vec=model_crit_vec,
              select_crit=select_crit))
}

