#general majorized ll function
nll_maj_func <- function(para, para0, y1, y2, delta1, delta2,
                         Xmat1, Xmat2, Xmat3,
                         hazard, frailty, model,
                         penalty, lambda, a,
                         penalty_fusedcoef, lambda_fusedcoef,
                         penalty_fusedbaseline, lambda_fusedbaseline, pen_mat_w_lambda,
                         mm_epsilon, penweights_list){

  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  nll <- nll_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                  frailty=frailty,hazard=hazard,model=model,
                  basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                  dbasis1=NULL, dbasis2=NULL, dbasis3=NULL)

  maj <- maj_func(para=para,para0=para0,nP1=nP1,nP2=nP2,nP3=nP3,
                  penalty=penalty,lambda=lambda, a=a,
                  penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                  penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline, pen_mat_w_lambda = pen_mat_w_lambda,
                  penweights_list=penweights_list, mm_epsilon=mm_epsilon, hazard=hazard)

  return( nll + (n * maj) )

}




#general majorized penalty negative gradient function
ngrad_maj_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                           hazard, frailty, model,
                           penalty, lambda, a,
                           penalty_fusedcoef, lambda_fusedcoef,
                           penalty_fusedbaseline, lambda_fusedbaseline, pen_mat_w_lambda,
                           penweights_list, mm_epsilon){

  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)


  ngrad <- ngrad_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2, Xmat1=Xmat1, Xmat2=Xmat2,
                      Xmat3=Xmat3, frailty=frailty, hazard=hazard, model=model,
                      basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                      dbasis1=NULL, dbasis2=NULL, dbasis3=NULL)

  Ek <- pen_maj_mat_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,
                         penalty=penalty,lambda=lambda, a=a,
                         penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                         penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                         penweights_list=penweights_list, mm_epsilon=mm_epsilon, hazard=hazard)

  return( ngrad + ( n * (Ek %*% para) ) )

}


nhess_maj_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                           hazard, frailty, model,
                           penalty, lambda, a,
                           penalty_fusedcoef, lambda_fusedcoef,
                           penalty_fusedbaseline, lambda_fusedbaseline, pen_mat_w_lambda,
                           penweights_list, mm_epsilon){

  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)


  nhess <- nhess_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2, Xmat1=Xmat1, Xmat2=Xmat2,
                      Xmat3=Xmat3, frailty=frailty, hazard=hazard, model=model)

  Ek <- pen_maj_mat_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,
                         penalty=penalty,lambda=lambda, a=a,
                         penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                         penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                         penweights_list=penweights_list, mm_epsilon=mm_epsilon, hazard=hazard)

  return( nhess + ( n * Ek ) )
}
