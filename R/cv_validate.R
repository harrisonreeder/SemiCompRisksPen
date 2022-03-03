#' Calculate Cross-validated prediction performance metrics
#'
#' This function calculates cross-validated prediction performance metrics
#'
#' @inheritParams FreqID_HReg2
#' @inheritParams calc_risk
#' @param n_folds Integer value specifying the number of partitions to divide the data into for cross validation
#' @param verbose Numeric indicating the amount of intermediate information that should be printed during the cross-validation process.
#'   Larger numbers correspond to more printed information.
#'
#' @return if Xmat has only one row, and t_cutoff is a scalar, then returns a 4 element row matrix
#'   of probabilities. If Xmat has \code{n} rows, then returns an \code{n} by 4 matrix of probabilities.
#'   If Xmat has \code{n} rows and t_cutoff is a vector of length \code{s}, then returns an \code{s} by 4 by \code{n} array.
#' @export
cv_predict <- function(Formula, data, na.action="na.fail", subset=NULL,
                         hazard=c("weibull"),frailty=TRUE, model,
                         knots_list = NULL, p0_vec = rep(4,3), startVals=NULL,
                         optim_method = if(tolower(hazard) %in% c("royston-parmar","rp")) "BFGS" else "L-BFGS-B",
                         n_folds, t_cutoff, t_start=0, tol=1e-3,
                         type="marginal", gamma=1, h3_tv, tv_knots,verbose=0){
  # browser()
  #na and subsetting actions first
  na.action <- match.arg(na.action) #technically na.fail I think is the only one currently implemented
  if (na.action != "na.fail" & na.action != "na.omit") {
    stop("na.action should be either na.fail or na.omit")
  }
  #rearrange input Formula object (which stores the different pieces of the input formula)
  #this seems to ensure it's in the correct form, doesn't seem to make a diff
  form2 <- Formula::as.Formula(paste(Formula[2], Formula[1], Formula[3],sep=""))
  #arrange data into correct form, extracting subset of variables
  #and observations needed in model
  data <- stats::model.frame(form2, data=data, na.action=na.action,subset=subset)

  #create matrices storing two outcomes, and then component vectors
  time1 <- Formula::model.part(form2, data=data, lhs=1)
  time2 <- Formula::model.part(form2, data=data, lhs=2)
  # Y <- cbind(time1[1], time1[2], time2[1], time2[2])
  y1 <- time1[[1]]
  delta1 <- time1[[2]]
  y2 <- time2[[1]]
  delta2 <- time2[[2]]



  N <- nrow(data)
  cv_index <- sample(rep(1:n_folds,length.out=N),size = N,replace = F)
  cv_fit_list <- list()
  cv_pred_list <- list()
  cv_outcome_list <- list()
  cv_ipcw_list <- list()

  for(leave_out in 1:n_folds){
    if(verbose>0){
      print("CV iteration")
      print(leave_out)
    }

    #divide test and training data
    training_data <- data[!(cv_index %in% leave_out),]
    test_data <- data[cv_index %in% leave_out,]
    y1_test <- y1[cv_index %in% leave_out]
    y2_test <- y2[cv_index %in% leave_out]
    delta1_test <- delta1[cv_index %in% leave_out]
    delta2_test <- delta2[cv_index %in% leave_out]

    #fit on training data
    fit_cv <- FreqID_HReg2(Formula = form2, data = training_data,
                             model=model,hazard=hazard,frailty=frailty,
                             knots_list = knots_list, p0_vec = p0_vec,hessian = FALSE,
                             startVals=startVals, optim_method = optim_method)
    cv_fit_list[[leave_out]] <- fit_cv
    #Create covariate matrices for each of three transition hazards
    Xmat1_test <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=1),
                                          data=test_data))
    Xmat2_test <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=2),
                                          data=test_data))
    Xmat3_test <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=3),
                                          data=test_data))

    #assume any time-varying quantities are at the end of the formula
    if(h3_tv == "linear"){
      Xmat3_test <- Xmat3_test[,-ncol(Xmat3_test)]
    } else if(h3_tv == "piecewise"){
      stopifnot(!is.null(tv_knots))
      if(tv_knots[1] != 0){tv_knots <- c(0,tv_knots)}
      if(utils::tail(tv_knots, n=1) != Inf){tv_knots <- c(tv_knots,Inf)}
      Xmat3_test <- Xmat3_test[,-((ncol(Xmat3_test) - length(tv_knots) + 3):ncol(Xmat3_test))]
    }

    cv_pred_list[[leave_out]] <- calc_risk(para = fit_cv$estimate,
                                              Xmat1 = Xmat1_test, Xmat2 = Xmat2_test, Xmat3 = Xmat3_test,
                                              t_cutoff = t_cutoff, t_start=t_start, tol=tol, frailty=frailty,
                                              type=type, gamma=gamma,model=model,hazard=hazard,
                                              h3_tv=h3_tv,tv_knots=tv_knots)
    cv_outcome_list[[leave_out]] <- get_outcome_mat(y1=y1_test,y2 = y2_test,
                                                       delta1 = delta1_test,delta2 = delta2_test,
                                                       t_cutoff = t_cutoff)
    cv_ipcw_list[[leave_out]] <- get_ipcw_mat(y2 = y2_test,delta2 = delta2_test,t_cutoff = t_cutoff)

  }

  list(data=data,cv_index=cv_index,cv_fit_list=cv_fit_list,
       cv_pred_list=cv_pred_list,
       cv_outcome_list=cv_outcome_list,
       cv_ipcw_list=cv_ipcw_list,
       prediction_details=list(
         para = fit_cv$estimate,
         Xmat1 = Xmat1_test, Xmat2 = Xmat2_test, Xmat3 = Xmat3_test,
         t_cutoff = t_cutoff, t_start=t_start, tol=tol, frailty=frailty,
         type=type, gamma=gamma,model=model,
         h3_tv=h3_tv,tv_knots=tv_knots
       ))

}

#' Calculate Cross-validated prediction performance metrics
#'
#' This function calculates cross-validated prediction performance metrics
#'
#' @param cv_predict_out_list is an object output from \code{cv_predict} function
#' @param metrics is a vector of strings naming what metrics to compute. Options are "brier", "auc", "hum", "ccp",and "pdi".
#' @param verbose Numeric indicating the amount of intermediate information that should be printed during the cross-validation process.
#'   Larger numbers correspond to more printed information.
#'
#' @return if Xmat has only one row, and t_cutoff is a scalar, then returns a 4 element row matrix
#'   of probabilities. If Xmat has \code{n} rows, then returns an \code{n} by 4 matrix of probabilities.
#'   If Xmat has \code{n} rows and t_cutoff is a vector of length \code{s}, then returns an \code{s} by 4 by \code{n} array.
#' @export
cv_validate <- function(cv_predict_out_list,
                        metrics=c("brier","auc","hum","ccp","pdi"),verbose=0){
  # browser()
  n_folds <- max(cv_predict_out_list$cv_index)
  t_cutoff <- cv_predict_out_list$prediction_details$t_cutoff
  metrics <- tolower(metrics)
  outlist <- list()
  for(metric  in metrics){
    outlist[[metric]] <- list()
  }


  for(leave_out in 1:n_folds){
    if(verbose>0){
      print("CV iteration")
      print(leave_out)
    }

    outcome_mat <- cv_predict_out_list$cv_outcome_list[[leave_out]]
    pred_mat <- cv_predict_out_list$cv_pred_list[[leave_out]]
    ipcw_mat <- cv_predict_out_list$cv_ipcw_list[[leave_out]]

    if("brier" %in% metrics){
      outlist[["brier"]][[leave_out]] <- compute_score(outcome_mat = outcome_mat,pred_mat = pred_mat, ipcw_mat = ipcw_mat, score="brier")
    }
    if("hum" %in% metrics){
      outlist[["hum"]][[leave_out]] <- compute_hum(outcome_mat = outcome_mat,pred_mat = pred_mat, ipcw_mat = ipcw_mat)
    }
    if("ccp" %in% metrics){
      outlist[["ccp"]][[leave_out]] <- compute_ccp(outcome_mat = outcome_mat,pred_mat = pred_mat, ipcw_mat = ipcw_mat)
    }
    if("pdi" %in% metrics){
      outlist[["pdi"]][[leave_out]] <- compute_pdi(outcome_mat = outcome_mat,pred_mat = pred_mat, ipcw_mat = ipcw_mat)
    }
    if("auc" %in% metrics){
      t_cutoff <- cv_predict_out_list$prediction_details$t_cutoff
      stopifnot(length(t_cutoff)==1)
      test_data <- cv_predict_out_list$data[cv_predict_out_list$cv_index %in% leave_out,]
      form_temp <- cv_predict_out_list$cv_fit_list[[1]]$formula
      #create matrices storing two outcomes, and then component vectors
      time1 <- Formula::model.part(form_temp, data=test_data, lhs=1)
      time2 <- Formula::model.part(form_temp, data=test_data, lhs=2)
      dat <- data.frame(y1=time1[[1]], delta1=time1[[2]], y2=time2[[1]], delta2=time2[[2]])

      outlist[["auc"]][[leave_out]] <- compute_auc(dat = dat,t_cutoff = t_cutoff,pred_mat = pred_mat)
    }
  }


  if(length(t_cutoff)==1){
    cv_metrics <- NULL
    lab <- NULL
    for(i in 1:length(metrics)){
      cv_metrics <- c(cv_metrics,colMeans(do.call(rbind,outlist[[metrics[i]]])))
      if(metrics[i]=="auc"){
        lab <- c(lab,"AUC_nt","AUC_t")
      } else{
        lab <- c(lab,metrics[i])
      }
    }
    outlist[["cv_metrics"]] <- matrix(cv_metrics,nrow = length(t_cutoff),byrow = TRUE,
                                      dimnames = list(as.character(t_cutoff),lab))
  }

  outlist
}


