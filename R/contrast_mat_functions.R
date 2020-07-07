#' Create List of Contrast Matrices for Parameter Fusion
#'
#' Creates list of numeric matrices corresponding to the fusion of blocks of parameters.
#'   The list is indexed by pairs of parameter blocks being fused, i.e.,
#'   \itemize{
#'   \item \code{fusedcoef12},\code{fusedcoef13},\code{fusedcoef23} list elements correspond to
#'     matrices that, when multiplied to the parameter vector, yield a vector of fused differences
#'     between corresponding regression parameters between the stated transition hazards.
#'   \item \code{fusedbaseline12},\code{fusedbaseline13},\code{fusedbaseline23} list elements
#'     are as above, but fusing the baseline hazard parameters.
#'   }
#'
#' @inheritParams proximal_gradient_descent
#' @param nP0,nP1,nP2,nP3 Number of parameters.
#'
#' @return Returns a list with elements named as above, each containing a numeric contrast matrix.
#' @export
contrast_mat_list <- function(nP0,nP1,nP2,nP3, #NEEDS UPDATE TO WORK WITH WEIBULL
                              penalty_fusedcoef,lambda_fusedcoef,
                              penalty_fusedbaseline,lambda_fusedbaseline,
                              hazard, penweights_list){

  ####FUNCTIONS TO CREATE CONTRAST MATRICES FOR FUSED PENALIZATION####

  #create single (weighted) contrast matrix that applies all of the differences at once
  #some circularity is induced because we need matrix to get adaptive weights to get weighted matrix, but that's ok

  pen_mat_list <-list()
  lambda_f_vec <- NULL
  nPtot <- nP0 + nP1 + nP2 + nP3

  if(tolower(penalty_fusedcoef) %in% c("fusedlasso","adafusedlasso")){

    if(length(lambda_fusedcoef)==1){
      lambda_fusedcoef12 <- lambda_fusedcoef13 <- lambda_fusedcoef23 <- lambda_fusedcoef
    } else if(length(lambda_fusedcoef)==3){
      lambda_fusedcoef12 <- lambda_fusedcoef[1]; lambda_fusedcoef13 <- lambda_fusedcoef[2];lambda_fusedcoef23 <- lambda_fusedcoef[3]
    } else{ stop("lambda_fusedcoef is neither a single value or a 3-vector!!") }

    #create a difference matrix connecting first and second sets of betas, left padded by 7 columns of zeros for the hazards
    if(lambda_fusedcoef12 != 0){
      stopifnot(nP1==nP2)
      #in principle, any penalty can have weights added, e.g., to set some terms to 0
      if(!is.null(penweights_list[["fusedcoef12"]]) && length(penweights_list[["fusedcoef12"]])==nP1){
        pen_mat_list[["fusedcoef12"]] <- cbind(matrix(data=0,nrow=nP1,ncol=nP0),
                                               diag(as.vector(penweights_list[["fusedcoef12"]])),
                                               -diag(as.vector(penweights_list[["fusedcoef12"]])),
                                               matrix(data=0,nrow=nP1,ncol=nP1))
      } else{
        pen_mat_list[["fusedcoef12"]] <- cbind(matrix(data=0,nrow=nP1,ncol=nP0), diag(nP1), -diag(nP1), matrix(data=0,nrow=nP1,ncol=nP1))
      }
      lambda_f_vec <- c(lambda_f_vec, rep(lambda_fusedcoef12, nP1))
    }

    #create a difference matrix connecting first and third sets of betas, left padded by 7 columns of zeros for the hazards
    if(lambda_fusedcoef13 != 0){
      stopifnot(nP1==nP3)
      #in principle, any penalty can have weights added, e.g., to set some terms to 0
      if(!is.null(penweights_list[["fusedcoef13"]]) && length(penweights_list[["fusedcoef13"]])==nP1){
        pen_mat_list[["fusedcoef13"]] <- cbind(matrix(data=0,nrow=nP1,ncol=nP0),
                                               diag(as.vector(penweights_list[["fusedcoef13"]])),
                                               matrix(data=0,nrow=nP1,ncol=nP1),
                                               -diag(as.vector(penweights_list[["fusedcoef13"]])))
      } else{
        pen_mat_list[["fusedcoef13"]] <- cbind(matrix(data=0,nrow=nP1,ncol=nP0), diag(nP1),matrix(data=0,nrow=nP1,ncol=nP1),-diag(nP1))
      }
      lambda_f_vec <- c(lambda_f_vec, rep(lambda_fusedcoef13, nP1))
    }

    #create a difference matrix connecting second and third sets of betas, left padded by 7 columns of zeros for the hazards
    if(lambda_fusedcoef23 != 0){
      stopifnot(nP2==nP3)
      #in principle, any penalty can have weights added, e.g., to set some terms to 0
      if(!is.null(penweights_list[["fusedcoef23"]]) && length(penweights_list[["fusedcoef23"]])==nP2){
        pen_mat_list[["fusedcoef23"]] <- cbind(matrix(data=0,nrow=nP2,ncol=nP0),
                                               matrix(data=0,nrow=nP2,ncol=nP2),
                                               diag(as.vector(penweights_list[["fusedcoef23"]])),
                                               -diag(as.vector(penweights_list[["fusedcoef23"]])))
      } else{
        pen_mat_list[["fusedcoef23"]] <- cbind(matrix(data=0,nrow=nP2,ncol=nP0), matrix(data=0,nrow=nP2,ncol=nP2),diag(nP2),-diag(nP2))
      }
      lambda_f_vec <- c(lambda_f_vec, rep(lambda_fusedcoef23, nP2))
    }

  }

  #Finally, add penalty of fusion for baseline parameters
  if(tolower(penalty_fusedbaseline) != "none" & !(is.null(hazard))){
    if(!(tolower(hazard) %in% c("weibull","wb"))){
      stop("non-Weibull not yet implemented")
    }

    if(length(lambda_fusedbaseline)==1){
      lambda_fusedbaseline12 <- lambda_fusedbaseline13 <- lambda_fusedbaseline23 <- lambda_fusedbaseline
    } else if(length(lambda_fusedbaseline)==3){
      lambda_fusedbaseline12 <- lambda_fusedbaseline[1]
      lambda_fusedbaseline13 <- lambda_fusedbaseline[2]
      lambda_fusedbaseline23 <- lambda_fusedbaseline[3]
    } else{ stop("lambda_fusedbaseline is neither a single value or a 3-vector!!") }

    #create a difference matrix connecting first and second sets of baseline parameters
    if(lambda_fusedbaseline12 != 0){
      D_temp <- matrix(data=0,nrow=2,ncol=nPtot)

      #FIGURE OUT A BETTER WAY TO INCORPORATE THESE WEIGHTS
      if(!is.null(penweights_list[["fusedbaseline12"]]) && length(penweights_list[["fusedbaseline12"]])==2){
        D_temp[1,1] <- D_temp[2,2] <- penweights_list[["fusedbaseline12"]][1]
        D_temp[1,3] <- D_temp[2,4] <- -penweights_list[["fusedbaseline12"]][2]
      } else{
        D_temp[1,1] <- D_temp[2,2] <- 1
        D_temp[1,3] <- D_temp[2,4] <- -1
      }
      pen_mat_list[["fusedbaseline12"]] <- D_temp
      lambda_f_vec <- c(lambda_f_vec, rep(lambda_fusedbaseline12, 2))
    }
    #create a difference matrix connecting first and third sets of baseline parameters
    if(lambda_fusedbaseline13 != 0){
      D_temp <- matrix(data=0,nrow=2,ncol=nPtot)

      if(!is.null(penweights_list[["fusedbaseline13"]]) && length(penweights_list[["fusedbaseline13"]])==2){
        D_temp[1,1] <- D_temp[2,2] <- penweights_list[["fusedbaseline13"]][1]
        D_temp[1,5] <- D_temp[2,6] <- -penweights_list[["fusedbaseline13"]][2]

      } else{
        D_temp[1,1] <- D_temp[2,2] <- 1
        D_temp[1,5] <- D_temp[2,6] <- -1
      }
      pen_mat_list[["fusedbaseline13"]] <- D_temp
      lambda_f_vec <- c(lambda_f_vec, rep(lambda_fusedbaseline13, 2))
    }
    #create a difference matrix connecting second and third sets of baseline parameters
    if(lambda_fusedbaseline23 != 0){
      D_temp <- matrix(data=0,nrow=2,ncol=nPtot)

      if(!is.null(penweights_list[["fusedbaseline23"]]) && length(penweights_list[["fusedbaseline23"]])==2){
        D_temp[1,3] <- D_temp[2,4] <- penweights_list[["fusedbaseline23"]]
        D_temp[1,5] <- D_temp[2,6] <- -penweights_list[["fusedbaseline23"]]
      } else{
        D_temp[1,3] <- D_temp[2,4] <- 1
        D_temp[1,5] <- D_temp[2,6] <- -1
      }
      pen_mat_list[["fusedbaseline23"]] <- D_temp
      lambda_f_vec <- c(lambda_f_vec, rep(lambda_fusedbaseline23, 2))
    }
  }

  return(list(pen_mat_list=pen_mat_list,lambda_f_vec=lambda_f_vec))
}



pen_mat_decomp <- function(pen_mat) {

  #FROM SMURF CODE "PENALTY_MATRICES"
  # Compute eigenvalue decomposition of t(pen_mat[[j]]) %*% pen_mat[[j]] for fused penalties
  # except "none", "lasso" and "grouplasso"
  #
  # pen_mat: (weighted) penalty matrix


  pen_mat_aux <- list()
  # Return NULL if error
  tmp <- tryCatch(eigen(t(pen_mat) %*% pen_mat), error = function(e) NULL)

  if (!is.null(tmp)) {
    # Get eigenvectors and -values if eigen did not give an error
    pen_mat_aux$Q <- tmp$vectors
    pen_mat_aux$eigval <- tmp$values

  } else {
    # eigen gave an error, use slower ADMM version (check happens when calling C++ code)
    pen_mat_aux$Q <- as.matrix(0)
    pen_mat_aux$eigval <- 0
  }
  return(pen_mat_aux)
}
