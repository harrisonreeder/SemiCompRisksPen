// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

arma::vec getCommonVec(const arma::vec& delta1,const arma::vec& delta2,const arma::vec& AVec, double h){
	return (delta1 + delta2 + exp(-h)) / (1 + exp(h) * AVec);
}


/*******
Royston-Parmar
*******/

// [[Rcpp::export]] 
double nlogLikRP_uni(const arma::vec& para, const arma::vec& delta, const arma::vec& y,
					const arma::mat& basis, const arma::mat& dbasis, const arma::mat& X){

	//define constants
	int p0 = basis.n_cols; 
	int p1 = X.n_cols; 
	int n = X.n_rows;

	arma::vec phi = para(arma::span(0, p0-1)); 

	//define linear predictors
	arma::vec eta;
	if(p1 == 0){
		eta = arma::zeros(n);
	} else{
		arma::vec beta =  para(arma::span(p0, p0+p1-1));
		eta = X * beta;
	}

	//define splines
	arma::vec s = basis * phi;
	arma::vec sprime = dbasis * phi;

	//delta * (log lambda) + log(survival)
	double obj_val = arma::accu( delta % (arma::log(sprime) + s + eta - arma::log(y)) - arma::exp(s + eta)  ); 
	//note that in theory the '-delta*log(y)' term can be proportioned out, but I added it in to keep a properly comparable -2LL scale for AIC

	return(-obj_val);
}



// [[Rcpp::export]] 
arma::vec ngradRP_uni(const arma::vec& para, const arma::vec& y, const arma::vec& delta, const arma::mat& basis, const arma::mat& dbasis, const arma::mat& X){

	//define constants
	int p0 = basis.n_cols; 
	int p1 = X.n_cols; 
	int n = X.n_rows;

	arma::vec phi = para(arma::span(0, p0-1)); 

	//define linear predictors
	arma::vec eta;
	if(p1 == 0){
		eta = arma::zeros(n);
	} else{
		arma::vec beta =  para(arma::span(p0, p0+p1-1));
		eta = X * beta;
	}

	//define splines
	arma::vec s = basis * phi;
	arma::vec sprime = dbasis * phi;

	arma::vec temp_scorevec = arma::zeros<arma::vec>(p0 + p1);

	temp_scorevec(arma::span(0, p0-1)) = dbasis.t() * (delta % arma::pow(sprime,-1)) + basis.t() * (delta - arma::exp(s + eta) );

	if(p1 > 0){
		temp_scorevec(arma::span(p0,p0+p1-1)) = X.t() * ( delta - arma::exp(s + eta) );
	}	
	return(-temp_scorevec);

}

// [[Rcpp::export]] 
double nlogLikRP_ID_frail_SM(const arma::vec& para, const arma::vec& y1, const arma::vec& y2, 
							 const arma::vec& delta1, const arma::vec& delta2, 
						   	 const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, 
						   	 const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3, 
						   	 const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3){
	//define constants
	int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	
	int n = X1.n_rows;

	arma::vec phi1 = para(arma::span(0, p01-1));
	arma::vec phi2 = para(arma::span(p01, p01+p02-1)); 	
	arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1)); 	
	double h = para(p01+p02+p03);

	//define linear predictors
	arma::vec eta1, eta2, eta3;
	if(p1 == 0){
		eta1 = arma::zeros(n);
	} else{
		eta1 = X1 * para(arma::span(p01+p02+p03+1, p01+p02+p03+p1));
	}
	if(p2 == 0){
		eta2 = arma::zeros(n);
	} else{
		eta2 = X2 * para(arma::span(p01+p02+p03+1+p1, p01+p02+p03+p1+p2));
	}
	if(p3 == 0){
		eta3 = arma::zeros(n);
	} else{
		eta3 = X3 * para(arma::span(p01+p02+p03+1+p1+p2, p01+p02+p03+p1+p2+p3));
	}

	//you might think that all of the information about y1 and y2 is coming through the basis, but
	//to keep the likelihood on the same scale as other parametric models, we need to include the 'extraneous'
	//terms involving y1, y2, and y2-y1 in the final calculation
	arma::vec logdiff = arma::log(y2-y1);
	logdiff = logdiff.replace(-arma::datum::inf, 0); //log(y2-y1) with the negative infinity values replaced with 0's.


	//define splines
	//assumptions: basis 1 is coming from y_1, with placement of knots depending on quantiles of y[delta1==1]
	arma::vec s1 = basis1 * phi1;
	arma::vec s1prime = dbasis1 * phi1;
	//assumptions: basis 2 is coming from y_1, with placement of knots depending on quantiles of y[delta1==0 & delta2==0]
	arma::vec s2 = basis2 * phi2;
	arma::vec s2prime = dbasis2 * phi2;
	//assumptions: under semi-markov basis 3 is coming from y_2-y1, with placement of knots depending on quantiles of (y2-y1)[delta1==1 & delta2==1]
	arma::vec s3 = basis3 * phi3;
	arma::vec s3prime = dbasis3 * phi3;

	//my big concern is what to do in the semi-markov setting, when log(y2-y1)=-Infty is not well characterized here?
	//If we look at where this comes into play, it arises in the likelihood term that is already multiplied by delta1*delta2, 
	//and then in AVec so really we just need to ensure that observations with delta1=0 are zeroed out in the final sum.
	arma::vec AVec;
	AVec = arma::exp(s1 + eta1) + arma::exp(s2 + eta2) + delta1 % arma::exp(s3 + eta3);


	double obj_val = arma::accu(  delta1 % ( arma::log(s1prime) + s1 + eta1 - arma::log(y1)) 
	+ (1-delta1) % delta2 % ( arma::log(s2prime) + s2 + eta2 - arma::log(y1) )
	+ delta1 % delta2 % ( arma::log(s3prime) + s3 + eta3 - logdiff + log1p( exp(h) )  )
	- (exp(-h) + delta1 + delta2) % arma::log1p( exp(h) * AVec )  ) ;

  return -obj_val;    
}


// [[Rcpp::export]] 
double nlogLikRP_ID_frail_M(const arma::vec& para, const arma::vec& y1, const arma::vec& y2,
							const arma::vec& delta1, const arma::vec& delta2, 
						    const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, 
						    const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3, const arma::mat& basis3_y1, 
						    const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3){
	//define constants
	int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	
	int n = X1.n_rows;

	arma::vec phi1 = para(arma::span(0, p01-1));
	arma::vec phi2 = para(arma::span(p01, p01+p02-1)); 	
	arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1)); 	
	double h = para(p01+p02+p03);

	//define linear predictors
	arma::vec eta1, eta2, eta3;
	if(p1 == 0){
		eta1 = arma::zeros(n);
	} else{
		eta1 = X1 * para(arma::span(p01+p02+p03+1, p01+p02+p03+p1));
	}
	if(p2 == 0){
		eta2 = arma::zeros(n);
	} else{
		eta2 = X2 * para(arma::span(p01+p02+p03+1+p1, p01+p02+p03+p1+p2));
	}
	if(p3 == 0){
		eta3 = arma::zeros(n);
	} else{
		eta3 = X3 * para(arma::span(p01+p02+p03+1+p1+p2, p01+p02+p03+p1+p2+p3));
	}

	//define splines
	//assumptions: basis 1 is coming from y_1, with placement of knots depending on quantiles of y[delta1==1]
	arma::vec s1 = basis1 * phi1;
	arma::vec s1prime = dbasis1 * phi1;
	//assumptions: basis 2 is coming from y_1, with placement of knots depending on quantiles of y[delta1==0 & delta2==0]
	arma::vec s2 = basis2 * phi2;
	arma::vec s2prime = dbasis2 * phi2;
	//assumptions: under markov, basis 3 is coming from y_2, with placement of knots depending on quantiles of y2[delta1==1 & delta2==1]
	arma::vec s3 = basis3 * phi3;
	arma::vec s3prime = dbasis3 * phi3;

	//if we want a markov approach, we also need extra piece, which is y_1 set at same knots of y2[delta1==1 & delta2==1]
	arma::vec s3_y1 = basis3_y1 * phi3;

	//my big concern is what to do in the semi-markov setting, when log(y2-y1)=-Infty is not well characterized here?
	//If we look at where this comes into play, it arises in the likelihood term that is already multiplied by delta1*delta2, 
	//and then in AVec so really we just need to ensure that observations with delta1=0 are zeroed out in the final sum.
	arma::vec AVec;
	AVec = arma::exp(s1 + eta1) + arma::exp(s2 + eta2) + delta1 % (arma::exp(s3 + eta3) - arma::exp(s3_y1 + eta3));


	double obj_val = arma::accu(  delta1 % ( arma::log(s1prime) + s1 + eta1 - arma::log(y1) ) 
	+ (1-delta1) % delta2 % ( arma::log(s2prime) + s2 + eta2 - arma::log(y1) )
	+ delta1 % delta2 % ( arma::log(s3prime) + s3 + eta3 - arma::log(y2) + log1p( exp(h) )  )
	- (exp(-h) + delta1 + delta2) % arma::log1p( exp(h) * AVec )  ) ;

  return -obj_val;    
}

//This function implements both markov and semi-markov specifications, but I think I'll need to write separate functions for the grad and hess
// [[Rcpp::export]] 
arma::vec ngradRP_ID_frail_SM(const arma::vec& para, const arma::vec& y1, const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2, 
						   	  const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, 
						      const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3, 
						      const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3){
	//define constants
	int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	
	int n = X1.n_rows;

	arma::vec phi1 = para(arma::span(0, p01-1));
	arma::vec phi2 = para(arma::span(p01, p01+p02-1)); 	
	arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1)); 	
	double h = para(p01+p02+p03);

	//define linear predictors
	arma::vec eta1, eta2, eta3;
	if(p1 == 0){
		eta1 = arma::zeros(n);
	} else{
		eta1 = X1 * para(arma::span(p01+p02+p03+1, p01+p02+p03+p1));
	}
	if(p2 == 0){
		eta2 = arma::zeros(n);
	} else{
		eta2 = X2 * para(arma::span(p01+p02+p03+1+p1, p01+p02+p03+p1+p2));
	}
	if(p3 == 0){
		eta3 = arma::zeros(n);
	} else{
		eta3 = X3 * para(arma::span(p01+p02+p03+1+p1+p2, p01+p02+p03+p1+p2+p3));
	}

	//define splines
	//assumptions: basis 1 is coming from y_1, with placement of knots depending on quantiles of y[delta1==1]
	arma::vec s1 = basis1 * phi1;
	arma::vec s1prime = dbasis1 * phi1;
	//assumptions: basis 2 is coming from y_1, with placement of knots depending on quantiles of y[delta1==0 & delta2==0]
	arma::vec s2 = basis2 * phi2;
	arma::vec s2prime = dbasis2 * phi2;
	//assumptions: under semi-markov basis 3 is coming from y_2-y1, with placement of knots depending on quantiles of (y2-y1)[delta1==1 & delta2==1]
	arma::vec s3 = basis3 * phi3;
	arma::vec s3prime = dbasis3 * phi3;

	//my big concern is what to do in the semi-markov setting, when log(y2-y1)=-Infty is not well characterized here?
	//If we look at where this comes into play, it arises in the likelihood term that is already multiplied by delta1*delta2, 
	//and then in AVec so really we just need to ensure that observations with delta1=0 are zeroed out in the final sum.
	arma::vec AVec = arma::exp(s1 + eta1) + arma::exp(s2 + eta2) + delta1 % arma::exp(s3 + eta3);

  	arma::vec temp_scorevec = arma::zeros<arma::vec>(1+p01+p02+p03+p1+p2+p3);

	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);


	//phi1
	temp_scorevec(arma::span(0, p01 - 1)) = dbasis1.t() * (delta1 % arma::pow(s1prime,-1)) 
											+ basis1.t() * (delta1 - commonVec % arma::exp(h + s1 + eta1) );

	//phi2
	temp_scorevec(arma::span(p01, p01 + p02 - 1)) = dbasis2.t() * ((1-delta1) % delta2 % arma::pow(s2prime,-1)) 
													+ basis2.t() * ((1-delta1) % delta2 - commonVec % arma::exp(h + s2 + eta2) );

	//phi3
	temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = dbasis3.t() * (delta1 % delta2 % arma::pow(s3prime,-1)) 
													+ basis3.t() * (delta1 % delta2 - 
														delta1 % commonVec % arma::exp(h + s3 + eta3) );  //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2-y1). Helps make it more robust

	//h
	temp_scorevec(p01 + p02 + p03) = arma::accu( exp(h)*( delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec ));

	//beta1 (what ina calls u2)
	if(p1 > 0){
		temp_scorevec(arma::span(1 + p01 + p02 + p03, 1 + p01 + p02 + p03 + p1 - 1)) = X1.t() * (delta1 - commonVec % arma::exp(h + s1 + eta1));
	}

	//beta2 (what ina calls u3)
	if(p2 > 0){
		temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1, 1 + p01 + p02 + p03 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % arma::exp(h + s2 + eta2)) ;
	}

	//beta3 (what ina calls u4)
	if(p3 > 0){
		temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1 + p2,1 + p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - 
								 delta1 % commonVec % arma::exp(s3 + h + eta3)); //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2-y1)
	}

    return -temp_scorevec;

}



//This function implements both markov and semi-markov specifications, but I think I'll need to write separate functions for the grad and hess
// [[Rcpp::export]] 
arma::vec ngradRP_ID_frail_M(const arma::vec& para, const arma::vec& y1, const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2, 
						     const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, 
						     const arma::mat& basis1, const arma::mat& basis2, const arma::mat& basis3, const arma::mat& basis3_y1,
						     const arma::mat& dbasis1, const arma::mat& dbasis2, const arma::mat& dbasis3){
	//define constants
	int p01 = basis1.n_cols; int p02 = basis2.n_cols;	int p03 = basis3.n_cols;
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	
	int n = X1.n_rows;

	arma::vec phi1 = para(arma::span(0, p01-1));
	arma::vec phi2 = para(arma::span(p01, p01+p02-1)); 	
	arma::vec phi3 = para(arma::span(p01+p02, p01+p02+p03-1)); 	
	double h = para(p01+p02+p03);

	//define linear predictors
	arma::vec eta1, eta2, eta3;
	if(p1 == 0){
		eta1 = arma::zeros(n);
	} else{
		eta1 = X1 * para(arma::span(p01+p02+p03+1, p01+p02+p03+p1));
	}
	if(p2 == 0){
		eta2 = arma::zeros(n);
	} else{
		eta2 = X2 * para(arma::span(p01+p02+p03+1+p1, p01+p02+p03+p1+p2));
	}
	if(p3 == 0){
		eta3 = arma::zeros(n);
	} else{
		eta3 = X3 * para(arma::span(p01+p02+p03+1+p1+p2, p01+p02+p03+p1+p2+p3));
	}

	//define splines
	//assumptions: basis 1 is coming from y_1, with placement of knots depending on quantiles of y[delta1==1]
	arma::vec s1 = basis1 * phi1;
	arma::vec s1prime = dbasis1 * phi1;
	//assumptions: basis 2 is coming from y_1, with placement of knots depending on quantiles of y[delta1==0 & delta2==0]
	arma::vec s2 = basis2 * phi2;
	arma::vec s2prime = dbasis2 * phi2;
	//assumptions: under markov, basis 3 is coming from y_2, with placement of knots depending on quantiles of y2[delta1==1 & delta2==1]
	arma::vec s3 = basis3 * phi3;
	arma::vec s3prime = dbasis3 * phi3;

	//if we want a markov approach, we also need extra piece, which is y_1 set at same knots of y2[delta1==1 & delta2==1]
	arma::vec s3_y1 = basis3_y1 * phi3;

	//my big concern is what to do in the semi-markov setting, when log(y2-y1)=-Infty is not well characterized here?
	//If we look at where this comes into play, it arises in the likelihood term that is already multiplied by delta1*delta2, 
	//and then in AVec so really we just need to ensure that observations with delta1=0 are zeroed out in the final sum.
	arma::vec AVec = arma::exp(s1 + eta1) + arma::exp(s2 + eta2) + delta1 % (arma::exp(s3 + eta3) - arma::exp(s3_y1 + eta3));

  	arma::vec temp_scorevec = arma::zeros<arma::vec>(1 + p01 + p02 + p03 + p1 + p2 + p3);

	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);

	//phi1
	temp_scorevec(arma::span(0, p01 - 1)) = dbasis1.t() * (delta1 % arma::pow(s1prime,-1)) 
											+ basis1.t() * (delta1 - commonVec % arma::exp(h + s1 + eta1) );

	//phi2
	temp_scorevec(arma::span(p01, p01 + p02 - 1)) = dbasis2.t() * ((1-delta1) % delta2 % arma::pow(s2prime,-1)) 
													+ basis2.t() * ((1-delta1) % delta2 - commonVec % arma::exp(h + s2 + eta2) );

	//phi3
	temp_scorevec(arma::span(p01 + p02, p01 + p02 + p03 - 1)) = dbasis3.t() * (delta1 % delta2 % arma::pow(s3prime, -1)) 
													+ basis3.t() * (delta1 % delta2) 
													- (basis3.t() * (delta1 % commonVec % arma::exp(h + s3 + eta3))
													   - basis3_y1.t() * (delta1 % commonVec % arma::exp(h + s3_y1 + eta3)));  //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2)-H(y1). Helps make it more robust.


	//h
	temp_scorevec(p01 + p02 + p03) = arma::accu( exp(h)*( delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec ));

	//beta1 (what ina calls u2)
	if(p1 > 0){
		temp_scorevec(arma::span(1 + p01 + p02 + p03, 1 + p01 + p02 + p03 + p1 - 1)) = X1.t() * (delta1 - commonVec % arma::exp(h + s1 + eta1));
	}

	//beta2 (what ina calls u3)
	if(p2 > 0){
		temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1, 1 + p01 + p02 + p03 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % arma::exp(h + s2 + eta2)) ;
	}

	//beta3 (what ina calls u4)
	if(p3 > 0){
		temp_scorevec(arma::span(1 + p01 + p02 + p03 + p1 + p2,1 + p01 + p02 + p03 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - 
								 delta1 % commonVec % arma::exp(h + eta3) % (arma::exp(s3) - arma::exp(s3_y1))); //note, extra delta1 multiplied to last term because only those terms contribute nonzero H(y2)-H(y1). Helps make it more robust.
	}

    return -temp_scorevec;

}



/*******
WEIBULL
*******/


arma::vec getLambda0gWB(const arma::vec& yg, double ag, double kg){
	return exp(kg) * arma::pow(yg, exp(ag));
}

arma::vec getAVecWB_SM(const arma::vec& y1,const arma::vec& y2,
					   double k1,double a1, double k2,double a2, double k3, double a3,
					   const arma::vec& eta1, const arma::vec& eta2, const arma::vec& eta3){
	//this is for semi-markov specification
	return getLambda0gWB(y1,a1,k1) % arma::exp(eta1) + getLambda0gWB(y1,a2,k2) % arma::exp(eta2) + getLambda0gWB(y2-y1,a3,k3) % arma::exp(eta3);
}

arma::vec getAVecWB_M(const arma::vec& y1,const arma::vec& y2, 
					  double k1,double a1, double k2,double a2, double k3, double a3, 
					  const arma::vec& eta1, const arma::vec& eta2, const arma::vec& eta3){
	//this is for markov specification
	return getLambda0gWB(y1,a1,k1) % arma::exp(eta1) + getLambda0gWB(y1,a2,k2) % arma::exp(eta2) + getLambda0gWB(y2,a3,k3) % arma::exp(eta3) - getLambda0gWB(y1,a3,k3) % arma::exp(eta3);
}









//This function implements both markov and semi-markov specifications, but I think I'll need to write separate functions for the grad and hess
// [[Rcpp::export]] 
double nlogLikWB_ID_frail(const arma::vec& para,
						 const arma::vec& y1,const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2, 
						 const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const std::string model = "semi-markov"){
	//define constants
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

	double k1 = para(0);
	double a1 = para(1);
	double k2 = para(2);
	double a2 = para(3);
	double k3 = para(4);
	double a3 = para(5);
	double h = para(6);

	//define linear predictors
	arma::vec eta1, eta2, eta3, AVec, logdiff;
	if(p1 == 0){
		eta1 = arma::zeros(n);
	} else{
		eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
	}
	if(p2 == 0){
		eta2 = arma::zeros(n);
	} else{
		eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
	}
	if(p3 == 0){
		eta3 = arma::zeros(n);
	} else{
		eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
	}

	//define collected terms
	if (model.compare("markov") == 0){
		AVec = getAVecWB_M(y1,y2, k1,a1,k2,a2,k3,a3, eta1, eta2, eta3);
		logdiff = arma::log(y2); //awkward name reflects the fact that the markov assumption sets log(y2|y1) = log(y2), so there's no 'difference' here
	} else {
		AVec = getAVecWB_SM(y1,y2, k1,a1,k2,a2,k3,a3, eta1, eta2, eta3);
		logdiff = arma::log(y2-y1);
		logdiff = logdiff.replace(-arma::datum::inf, 0); //log(y2-y1) with the negative infinity values replaced with 0's.
	}

	//this is the slow way of writing it, I could simplify to make this a bit quicker
	arma::vec LL1 = delta1 % delta2 % (  a1 + k1 + (exp(a1) - 1) * arma::log(y1) + eta1 
									   + a3 + k3 + (exp(a3) - 1) * logdiff       + eta3
									   + log1p(exp(h)) //added because I think it's missing
									   - (exp(-h) + 2) * arma::log1p(exp(h) * AVec)
									  ); //LL1
	arma::vec LL3 = delta1 % (1-delta2) % (a1 + k1 + (exp(a1) - 1) * arma::log(y1) + eta1 
											- (exp(-h)+1) * arma::log1p(exp(h) * AVec)); //LL3	;
	arma::vec LL2 = (1-delta1) % delta2 % (a2 + k2 + (exp(a2) - 1) * arma::log(y1) + eta2 
											- (exp(-h)+1) * arma::log1p(exp(h) * AVec)); // LL 2
	arma::vec LL4 = (1-delta1) % (1-delta2) % (-exp(-h) * arma::log1p(exp(h) * AVec)); // LL4

    double obj_val = arma::accu( LL1 + LL2 + LL3 + LL4);
    return -obj_val;
}




// [[Rcpp::export]] 
arma::vec ngradWB_ID_frail_SM(const arma::vec& para,
							  const arma::vec& y1,const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2, 
							  const arma::mat& X1, const arma::mat& X2, const arma::mat& X3){
	//define constants
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

	double k1 = para(0);
	double a1 = para(1);
	double k2 = para(2);
	double a2 = para(3);
	double k3 = para(4);
	double a3 = para(5);
	double h = para(6);

	//define linear predictors
	arma::vec eta1, eta2, eta3;
	if(p1 == 0){
		eta1 = arma::zeros(n);
	} else{
		eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
	}
	if(p2 == 0){
		eta2 = arma::zeros(n);
	} else{
		eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
	}
	if(p3 == 0){
		eta3 = arma::zeros(n);
	} else{
		eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
	}

	arma::vec temp_scorevec = arma::zeros<arma::vec>(p1+p2+p3+7);

	arma::vec AVec = getAVecWB_SM(y1, y2, k1, a1, k2, a2, k3, a3, eta1, eta2, eta3);
	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);
	arma::vec Lambda01 = getLambda0gWB(y1, a1, k1);
	arma::vec Lambda02 = getLambda0gWB(y1, a2, k2);
	arma::vec Lambda03 = getLambda0gWB(y2-y1, a3, k3);

	arma::vec logdiff = arma::log(y2-y1);
	logdiff = logdiff.replace(-arma::datum::inf, 0); //log(y2-y1) with the negative infinity values replaced with 0's.


	//k1 (what ina calls u5)
	temp_scorevec(0) = arma::accu(delta1 - commonVec % Lambda01 % arma::exp(h + eta1));

	//a1 (what ina calls u6)
	temp_scorevec(1) = arma::accu( delta1 % (1 + exp(a1) * arma::log(y1)) -
									exp(a1) * commonVec % Lambda01 % arma::exp(h+eta1) % arma::log(y1));
	//k2 (what ina calls u7)
	temp_scorevec(2) = arma::accu( (1-delta1) % delta2 - 
									commonVec % Lambda02 % arma::exp(h+eta2));
	//a2 (what ina calls u8)
	temp_scorevec(3) = arma::accu( (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) - 
									commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1) );
	//k3 (what ina calls u9)
	temp_scorevec(4) = arma::accu( delta1 % delta2 - 
									commonVec % Lambda03 % arma::exp(h+eta3));

	//a3 (what ina calls u10)
	temp_scorevec(5) = arma::accu( delta1 % delta2 % (1 + exp(a3) * logdiff) - 
									commonVec % Lambda03 % arma::exp(h+a3+eta3) % logdiff );
	//h (what ina calls u1)
	temp_scorevec(6) = arma::accu(exp(h)*(delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec));

	//beta1 (what ina calls u2)
	if(p1 > 0){
		temp_scorevec(arma::span(7, 7 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	}

	//beta2 (what ina calls u3)
	if(p2 > 0){
		temp_scorevec(arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
	}

	//beta3 (what ina calls u4)
	if(p3 > 0){
		temp_scorevec(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % Lambda03 % arma::exp(h + eta3));
	}

    return -temp_scorevec;
}

//this is the gradient with the markov assumption

// [[Rcpp::export]] 
arma::vec ngradWB_ID_frail_M(const arma::vec& para,
							 const arma::vec& y1,const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2,
							 const arma::mat& X1, const arma::mat& X2, const arma::mat& X3){
	//define constants
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

	double k1 = para(0);
	double a1 = para(1);
	double k2 = para(2);
	double a2 = para(3);
	double k3 = para(4);
	double a3 = para(5);
	double h = para(6);

	//define linear predictors
	arma::vec eta1, eta2, eta3;
	if(p1 == 0){
		eta1 = arma::zeros(n);
	} else{
		eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
	}
	if(p2 == 0){
		eta2 = arma::zeros(n);
	} else{
		eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
	}
	if(p3 == 0){
		eta3 = arma::zeros(n);
	} else{
		eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
	}

	arma::vec temp_scorevec = arma::zeros<arma::vec>(p1+p2+p3+7);

	arma::vec AVec = getAVecWB_M(y1, y2, k1, a1, k2, a2, k3, a3, eta1, eta2, eta3);
	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);
	arma::vec Lambda01 = getLambda0gWB(y1, a1, k1);
	arma::vec Lambda02 = getLambda0gWB(y1, a2, k2);
	arma::vec Lambda03y2 = getLambda0gWB(y2, a3, k3);
	arma::vec Lambda03y1 = getLambda0gWB(y1, a3, k3);

	//k1 (what ina calls u5)
	temp_scorevec(0) = arma::accu(delta1 - commonVec % Lambda01 % arma::exp(h + eta1));

	//a1 (what ina calls u6)
	temp_scorevec(1) = arma::accu( delta1 % (1 + exp(a1) * arma::log(y1)) -
									exp(a1) * commonVec % Lambda01 % arma::exp(h+eta1) % arma::log(y1));
	//k2 (what ina calls u7)
	temp_scorevec(2) = arma::accu( (1-delta1) % delta2 - 
									commonVec % Lambda02 % arma::exp(h+eta2));
	//a2 (what ina calls u8)
	temp_scorevec(3) = arma::accu( (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) - 
									commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1) );

	//k3 (what ina calls u9)
	temp_scorevec(4) = arma::accu( delta1 % delta2 - 
									commonVec % (Lambda03y2 - Lambda03y1) % arma::exp(h+eta3));

	//a3 (what ina calls u10)
	temp_scorevec(5) = arma::accu( delta1 % delta2 % (1 + exp(a3) * arma::log(y2)) - 
									commonVec % arma::exp(h+a3+eta3) % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) );
	//h (what ina calls u1)
	temp_scorevec(6) = arma::accu(exp(h)*(delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec));

	//beta1 (what ina calls u2)
	if(p1 > 0){
		temp_scorevec(arma::span(7, 7 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	}

	//beta2 (what ina calls u3)
	if(p2 > 0){
		temp_scorevec(arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
	}

	//beta3 (what ina calls u4)
	if(p3 > 0){
		temp_scorevec(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % (Lambda03y2 - Lambda03y1)  % arma::exp(h + eta3));
	}

    return -temp_scorevec;
}


//this is the hessian with the semi-markov assumption

// [[Rcpp::export]] 
arma::mat nhessWB_ID_frail_SM(const arma::vec& para,
							  const arma::vec& y1,const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2,
							  const arma::mat& X1, const arma::mat& X2, const arma::mat& X3){
	//define constants
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

	double k1 = para(0);
	double a1 = para(1);
	double k2 = para(2);
	double a2 = para(3);
	double k3 = para(4);
	double a3 = para(5);
	double h = para(6);

	//define linear predictors
	arma::vec eta1, eta2, eta3;
	if(p1 == 0){
		eta1 = arma::zeros(n);
	} else{
		eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
	}
	if(p2 == 0){
		eta2 = arma::zeros(n);
	} else{
		eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
	}
	if(p3 == 0){
		eta3 = arma::zeros(n);
	} else{
		eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
	}

	arma::mat temp_hessmat = arma::zeros<arma::mat>(p1+p2+p3+7,p1+p2+p3+7);

	arma::vec AVec = getAVecWB_SM(y1, y2, k1, a1, k2, a2, k3, a3, eta1, eta2, eta3);
	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);
	arma::vec Lambda01 = getLambda0gWB(y1, a1, k1);
	arma::vec Lambda02 = getLambda0gWB(y1, a2, k2);
	arma::vec Lambda03 = getLambda0gWB(y2-y1, a3, k3);

	arma::vec logdiff = arma::log(y2-y1);
	logdiff = logdiff.replace(-arma::datum::inf, 0); //log(y2-y1) with the negative infinity values replaced with 0's.

	arma::vec dAVecdk1 = Lambda01 % arma::exp(eta1);
	arma::vec dAVecdk2 = Lambda02 % arma::exp(eta2);
	arma::vec dAVecdk3 = Lambda03 % arma::exp(eta3);

	arma::vec dAVecda1 = Lambda01 % arma::exp(a1 + eta1) % arma::log(y1);
	arma::vec dAVecda2 = Lambda02 % arma::exp(a2 + eta2) % arma::log(y1);
	arma::vec dAVecda3 = Lambda03 % arma::exp(a3 + eta3) % logdiff;
	
	arma::vec dcommonVecdk1 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecdk1;
	arma::vec dcommonVecdk2 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecdk2;
	arma::vec dcommonVecdk3 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecdk3;

	arma::vec dcommonVecda1 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecda1;
	arma::vec dcommonVecda2 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecda2;
	arma::vec dcommonVecda3 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecda3;

//	double dAVecdh = 0;

	//***********k1 section***********//
	//temp_scorevec(0) = arma::accu(delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	
	//k1k1
	temp_hessmat(0,0) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecdk1 % Lambda01 + commonVec % Lambda01 ) );
	//k1a1
	temp_hessmat(0,1) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecda1 % Lambda01 + exp(a1) * commonVec % Lambda01 % arma::log(y1) ) );
	temp_hessmat(1,0) = temp_hessmat(0,1);
	//k1k2
	temp_hessmat(0,2) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecdk2 % Lambda01 ) );
	temp_hessmat(2,0) = temp_hessmat(0,2);
	//k1a2
	temp_hessmat(0,3) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecda2 % Lambda01 ) );
	temp_hessmat(3,0) = temp_hessmat(0,3);
	//k1k3
	temp_hessmat(0,4) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecdk3 % Lambda01 ) );
	temp_hessmat(4,0) = temp_hessmat(0,4);
	//k1a3
	temp_hessmat(0,5) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecda3 % Lambda01 ) );
	temp_hessmat(5,0) = temp_hessmat(0,5);
	//k1h
	temp_hessmat(6,0) = arma::accu( dAVecdk1 / (1 + exp(h)*AVec) - 
							exp(h) * (AVec % dcommonVecdk1 + commonVec % dAVecdk1 ));
	temp_hessmat(0,6) = temp_hessmat(6,0);

	//***********k2 section***********//
	//temp_scorevec(2) = arma::accu( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h+eta2));

	//k2k2
	temp_hessmat(2,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecdk2 % Lambda02 + commonVec % Lambda02 ) );
	//k2a1
	temp_hessmat(1,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecda1 % Lambda02 ) );
	temp_hessmat(2,1) = temp_hessmat(1,2);
	//k2a2
	temp_hessmat(3,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecda2 % Lambda02 + exp(a2) * commonVec % Lambda02 % arma::log(y1) ) );
	temp_hessmat(2,3) = temp_hessmat(3,2);
	//k2k3
	temp_hessmat(4,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecdk3 % Lambda02 ) );
	temp_hessmat(2,4) = temp_hessmat(4,2);
	//k2a3
	temp_hessmat(5,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecda3 % Lambda02 ) );
	temp_hessmat(2,5) = temp_hessmat(5,2);

	//***********k3 section***********// 
	//temp_scorevec(4) = arma::accu( delta1 % delta2 - commonVec % Lambda03 % arma::exp(h+eta3));

	//k3k3
	temp_hessmat(4,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecdk3 % Lambda03 + commonVec % Lambda03 ) );
	//k3a1
	temp_hessmat(1,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecda1 % Lambda03 ) );
	temp_hessmat(4,1) = temp_hessmat(1,4);
	//k3a2
	temp_hessmat(3,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecda2 % Lambda03 ) );
	temp_hessmat(4,3) = temp_hessmat(3,4);
	//k3a3
	temp_hessmat(5,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecda3 % Lambda03 + exp(a3) * commonVec % Lambda03 % logdiff ) );
	temp_hessmat(4,5) = temp_hessmat(5,4);

	//***********a1 section**********// 
	//temp_scorevec(1) = arma::accu( delta1 % (1 + exp(a1) * arma::log(y1)) - exp(a1) * commonVec % Lambda01 % arma::exp(h+eta1) % arma::log(y1));
	//rewritten as 		 arma::accu( delta1 + exp(a1) * delta1 % arma::log(y1) - arma::exp(h+eta1) % arma::log(y1) % (exp(a1) * commonVec % Lambda01) ); 

	//a1a1
	temp_hessmat(1,1) = arma::accu( exp(a1) * delta1 % arma::log(y1) - 
									arma::exp(h+eta1) % arma::log(y1) % ( exp(a1) * commonVec % Lambda01 + 
																		  exp(a1) * dcommonVecda1 % Lambda01 + 
																		  exp(2 * a1) * commonVec % Lambda01 % arma::log(y1) ) ); 
	//a1a2
	temp_hessmat(1,3) = arma::accu( - arma::exp(h+eta1) % arma::log(y1) % (exp(a1) * dcommonVecda2 % Lambda01) ); 
	temp_hessmat(3,1) = temp_hessmat(1,3);
	//a1a3
	temp_hessmat(1,5) = arma::accu( - arma::exp(h+eta1) % arma::log(y1) % (exp(a1) * dcommonVecda3 % Lambda01) ); 
	temp_hessmat(5,1) = temp_hessmat(1,5);

	//***********a2 section**********// 
	//temp_scorevec(3) = arma::accu( (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) - commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1) );

	//rewritten as arma::accu( (1-delta1) % delta2 + exp(a2) * (1-delta1) % delta2 % arma::log(y1) - arma::exp(h+eta2) % arma::log(y1) % (exp(a2) * commonVec % Lambda02) );
	//a2a2
	temp_hessmat(3,3) = arma::accu( exp(a2) * (1-delta1) % delta2 % arma::log(y1) - arma::exp(h+eta2) % arma::log(y1) % (exp(a2) * commonVec % Lambda02 + exp(a2) * dcommonVecda2 % Lambda02 + exp(2 * a2) * commonVec % Lambda02 % arma::log(y1)) );
	//a2a3
	temp_hessmat(3,5) = arma::accu( - arma::exp(h+eta2) % arma::log(y1) % (exp(a2) * dcommonVecda3 % Lambda02) ); 
	temp_hessmat(5,3) = temp_hessmat(3,5);

	//***********a3 section**********// 
	//temp_scorevec(5) = arma::accu( delta1 % delta2 % (1 + exp(a3) * logdiff) - commonVec % Lambda03 % arma::exp(h+a3+eta3) % logdiff );

	//a3a3
	temp_hessmat(5,5) = arma::accu( exp(a3) * delta1 % delta2 % logdiff - arma::exp(h+eta3) % logdiff % (exp(a3) * commonVec % Lambda03 + exp(a3) * dcommonVecda3 % Lambda03 + exp(2 * a3) * commonVec % Lambda03 % logdiff) );

	//***********h section***********//
	//temp_scorevec(6) = arma::accu(exp(h)*(delta1 % delta2 / (1+exp(h)) + arma::log(1+exp(h) * AVec)/exp(2*h) - commonVec % AVec) );

	//hh
	// from wolfram alpha (as an alternate derivation, with delta1=a and delta2=b): 
	// https://www.wolframalpha.com/input/?i=derivative+of+e%5Eh+*+(+a*b+%2F+(1%2Be%5Eh)+%2B+log(1+%2B+e%5Eh+*+A)+%2F+(e%5E(2*h))+-+(e%5E(-h)+%2B+a+%2B+b)+*+A+%2F+(1+%2B+e%5Eh+*+A))+with+respect+to+h
	temp_hessmat(6,6) = arma::accu( (-delta1 + 2 * AVec - delta2) / (AVec * exp(h)+1) 
								  + (delta1 - AVec + delta2) / (arma::pow(AVec * exp(h) + 1 , 2))
								  + delta1 % delta2 / (exp(h) + 1)
								  - delta1 % delta2 / (pow(exp(h) + 1, 2))
								  - exp(-h) * arma::log1p(AVec * exp(h))
								  );	
	//ha1
	temp_hessmat(6,1) = arma::accu( dAVecda1 / (1 + exp(h)*AVec) - 
							exp(h) * (AVec % dcommonVecda1 + commonVec % dAVecda1 ));
	temp_hessmat(1,6) = temp_hessmat(6,1);
	//hk2
	temp_hessmat(6,2) = arma::accu( dAVecdk2 / (1 + exp(h)*AVec) - 
							exp(h) * (AVec % dcommonVecdk2 + commonVec % dAVecdk2 ));
	temp_hessmat(2,6) = temp_hessmat(6,2);
	//ha2
	temp_hessmat(6,3) = arma::accu( dAVecda2 / (1 + exp(h)*AVec) - 
							exp(h) * (AVec % dcommonVecda2 + commonVec % dAVecda2 ));
	temp_hessmat(3,6) = temp_hessmat(6,3);
	//hk3
	temp_hessmat(6,4) = arma::accu( dAVecdk3 / (1 + exp(h)*AVec) - 
							exp(h) * (AVec % dcommonVecdk3 + commonVec % dAVecdk3 ));
	temp_hessmat(4,6) = temp_hessmat(6,4);
	//ha3
	temp_hessmat(6,5) = arma::accu( dAVecda3 / (1 + exp(h)*AVec) - 
							exp(h) * (AVec % dcommonVecda3 + commonVec % dAVecda3 ));
	temp_hessmat(5,6) = temp_hessmat(6,5);

	//***********beta1 section******//
	//temp_scorevec(arma::span(7, 7 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	if(p1 > 0){
		//beta1k1
		temp_hessmat(arma::span(7,7 + p1 - 1),0) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecdk1 % Lambda01 + commonVec % Lambda01 ) );
	    temp_hessmat(0,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),0).t();
		//beta1a1
		temp_hessmat(arma::span(7,7 + p1 - 1),1) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecda1 % Lambda01 + exp(a1) * commonVec % Lambda01 % arma::log(y1) ) );
	    temp_hessmat(1,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),1).t();
		//beta1k2
		temp_hessmat(arma::span(7,7 + p1 - 1),2) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecdk2 % Lambda01 ) );
	    temp_hessmat(2,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),2).t();
		//beta1a2
		temp_hessmat(arma::span(7,7 + p1 - 1),3) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecda2 % Lambda01 ) );
	    temp_hessmat(3,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),3).t();
		//beta1k3
		temp_hessmat(arma::span(7,7 + p1 - 1),4) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecdk3 % Lambda01 ) );
	    temp_hessmat(4,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),4).t();
		//beta1a3
		temp_hessmat(arma::span(7,7 + p1 - 1),5) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecda3 % Lambda01 ) );
	    temp_hessmat(5,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),5).t();
	    //beta1h
	    temp_hessmat(arma::span(7,7 + p1 - 1),6) = X1.t() * ( - Lambda01 % arma::exp(h + eta1) % (commonVec - (commonVec % AVec * exp(h) + exp(-h) ) / (1+exp(h)*AVec) ) );
	    temp_hessmat(6,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),6).t();

	    //beta1beta1
	    temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7,7 + p1 - 1)) = X1.t() * (X1.each_col() % ( - (Lambda01 % arma::exp(h + eta1) % ( commonVec + dcommonVecdk1 )))); //computes sum of w_i * x_ix_i^T

		if(p2 > 0){
		    //beta1beta2
		    //from univariate case: temp_hess(arma::span(2,p+1),arma::span(2,p+1)) = -( X.t() * (X.each_col() % (exp(k)*arma::exp(eta) % arma::pow(y,exp(a))))); //computes sum of w_i * x_ix_i^T
		    temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)) = X1.t() * (X2.each_col() % ( - dcommonVecdk2 % Lambda01 % arma::exp(h + eta1))); //computes sum of w_i * x_ix_i^T
			temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)).t();
		}

		if(p3 > 0){
		    //beta1beta3
		    //from univariate case: temp_hess(arma::span(2,p+1),arma::span(2,p+1)) = -( X.t() * (X.each_col() % (exp(k)*arma::exp(eta) % arma::pow(y,exp(a))))); //computes sum of w_i * x_ix_i^T
		    temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = X1.t() * (X3.each_col() % ( - dcommonVecdk3 % Lambda01 % arma::exp(h + eta1))); //computes sum of w_i * x_ix_i^T
			temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)).t();
		}
	}

	//***********beta2 section******//
    //temp_scorevec(arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
	if(p2 > 0){
		//beta2k1
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),0) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecdk1 % Lambda02 ) );
	    temp_hessmat(0,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),0).t();
		//beta2a1
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),1) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecda1 % Lambda02 ) );
	    temp_hessmat(1,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),1).t();
		//beta2k2
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),2) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecdk2 % Lambda02 + commonVec % Lambda02) );
	    temp_hessmat(2,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),2).t();
		//beta2a2
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),3) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecda2 % Lambda02 + exp(a2) * commonVec % Lambda02 % arma::log(y1) ) );
	    temp_hessmat(3,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),3).t();
		//beta2k3
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),4) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecdk3 % Lambda02 ) );
	    temp_hessmat(4,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),4).t();
		//beta2a3
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),5) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecda3 % Lambda02 ) );
	    temp_hessmat(5,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),5).t();
	    //beta2h
	    temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),6) = X2.t() * ( - Lambda02 % arma::exp(h + eta2) % (commonVec - (commonVec % AVec * exp(h) + exp(-h) ) / (1+exp(h)*AVec) ) );
	    temp_hessmat(6,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),6).t();

	    //beta2beta2
	    temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)) = X2.t() * (X2.each_col() % ( - (Lambda02 % arma::exp(h + eta2) % ( commonVec + dcommonVecdk2 )))); //computes sum of w_i * x_ix_i^T
	    if(p3 > 0){
		    //beta2beta3
		    //from univariate case: temp_hess(arma::span(2,p+1),arma::span(2,p+1)) = -( X.t() * (X.each_col() % (exp(k)*arma::exp(eta) % arma::pow(y,exp(a))))); //computes sum of w_i * x_ix_i^T
		    temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = X2.t() * (X3.each_col() % ( - dcommonVecdk3 % Lambda02 % arma::exp(h + eta2))); //computes sum of w_i * x_ix_i^T
			temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)).t();
		}
	}

	//***********beta3 section******//
    //temp_scorevec(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % Lambda03 % arma::exp(h + eta3));
	if(p3 > 0){
		//beta3k1
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),0) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecdk1 % Lambda03 ) );
	    temp_hessmat(0,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),0).t();
		//beta3a1
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),1) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecda1 % Lambda03 ) );
	    temp_hessmat(1,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),1).t();
		//beta3k2
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),2) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecdk2 % Lambda03 ) );
	    temp_hessmat(2,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),2).t();
		//beta3a2
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),3) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecda2 % Lambda03 ) );
	    temp_hessmat(3,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),3).t();
		//beta3k3
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),4) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecdk3 % Lambda03  + commonVec % Lambda03 ) );
	    temp_hessmat(4,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),4).t();
		//beta3a3
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),5) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecda3 % Lambda03 + exp(a3) * commonVec % Lambda03 % logdiff ) );
	    temp_hessmat(5,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),5).t();
	    //beta3h
	    temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),6) = X3.t() * ( - Lambda03 % arma::exp(h + eta3) % (commonVec - (commonVec % AVec * exp(h) + exp(-h) ) / (1+exp(h)*AVec) ) );
	    temp_hessmat(6,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),6).t();    

	    //beta3beta3
	    temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = X3.t() * (X3.each_col() % ( - (Lambda03 % arma::exp(h + eta3) % ( commonVec + dcommonVecdk3 )))); //computes sum of w_i * x_ix_i^T
	}

	return -temp_hessmat;
}








//this is the hessian with the markov assumption

// [[Rcpp::export]] 
arma::mat nhessWB_ID_frail_M(const arma::vec& para,
							 const arma::vec& y1,const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2,
							 const arma::mat& X1, const arma::mat& X2, const arma::mat& X3){

	//define constants
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

	double k1 = para(0);
	double a1 = para(1);
	double k2 = para(2);
	double a2 = para(3);
	double k3 = para(4);
	double a3 = para(5);
	double h = para(6);

	//define linear predictors
	arma::vec eta1, eta2, eta3;
	if(p1 == 0){
		eta1 = arma::zeros(n);
	} else{
		eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
	}
	if(p2 == 0){
		eta2 = arma::zeros(n);
	} else{
		eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
	}
	if(p3 == 0){
		eta3 = arma::zeros(n);
	} else{
		eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
	}

	arma::mat temp_hessmat = arma::zeros<arma::mat>(p1+p2+p3+7,p1+p2+p3+7);

	arma::vec AVec = getAVecWB_M(y1, y2, k1, a1, k2, a2, k3, a3, eta1, eta2, eta3);
	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);
	arma::vec Lambda01 = getLambda0gWB(y1, a1, k1);
	arma::vec Lambda02 = getLambda0gWB(y1, a2, k2);
	arma::vec Lambda03y2 = getLambda0gWB(y2, a3, k3);
	arma::vec Lambda03y1 = getLambda0gWB(y1, a3, k3);
	arma::vec Lambda03diff = Lambda03y2 - Lambda03y1;

	arma::vec dAVecdk1 = Lambda01 % arma::exp(eta1);
	arma::vec dAVecdk2 = Lambda02 % arma::exp(eta2);
	arma::vec dAVecdk3 = Lambda03diff % arma::exp(eta3);

	arma::vec dAVecda1 = Lambda01 % arma::exp(a1 + eta1) % arma::log(y1);
	arma::vec dAVecda2 = Lambda02 % arma::exp(a2 + eta2) % arma::log(y1);
	arma::vec dAVecda3 = arma::exp(a3 + eta3) % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1));
	
	arma::vec dcommonVecdk1 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecdk1;
	arma::vec dcommonVecdk2 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecdk2;
	arma::vec dcommonVecdk3 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecdk3;

	arma::vec dcommonVecda1 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecda1;
	arma::vec dcommonVecda2 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecda2;
	arma::vec dcommonVecda3 = -(exp(h) * commonVec) / (1 + exp(h)*AVec) % dAVecda3;

//	double dAVecdh = 0;

	//***********k1 section***********//
	//temp_scorevec(0) = arma::accu(delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	
	//k1k1
	temp_hessmat(0,0) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecdk1 % Lambda01 + commonVec % Lambda01 ) );
	//k1a1
	temp_hessmat(0,1) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecda1 % Lambda01 + exp(a1) * commonVec % Lambda01 % arma::log(y1) ) );
	temp_hessmat(1,0) = temp_hessmat(0,1);
	//k1k2
	temp_hessmat(0,2) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecdk2 % Lambda01 ) );
	temp_hessmat(2,0) = temp_hessmat(0,2);
	//k1a2
	temp_hessmat(0,3) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecda2 % Lambda01 ) );
	temp_hessmat(3,0) = temp_hessmat(0,3);
	//k1k3
	temp_hessmat(0,4) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecdk3 % Lambda01 ) );
	temp_hessmat(4,0) = temp_hessmat(0,4);
	//k1a3
	temp_hessmat(0,5) = arma::accu( - arma::exp(h + eta1) % ( dcommonVecda3 % Lambda01 ) );
	temp_hessmat(5,0) = temp_hessmat(0,5);
	//k1h
	temp_hessmat(6,0) = arma::accu( dAVecdk1 / (1 + exp(h)*AVec) - 
							exp(h) * (AVec % dcommonVecdk1 + commonVec % dAVecdk1 ));
	temp_hessmat(0,6) = temp_hessmat(6,0);

	//***********k2 section***********//
	//temp_scorevec(2) = arma::accu( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h+eta2));

	//k2k2
	temp_hessmat(2,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecdk2 % Lambda02 + commonVec % Lambda02 ) );
	//k2a1
	temp_hessmat(1,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecda1 % Lambda02 ) );
	temp_hessmat(2,1) = temp_hessmat(1,2);
	//k2a2
	temp_hessmat(3,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecda2 % Lambda02 + exp(a2) * commonVec % Lambda02 % arma::log(y1) ) );
	temp_hessmat(2,3) = temp_hessmat(3,2);
	//k2k3
	temp_hessmat(4,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecdk3 % Lambda02 ) );
	temp_hessmat(2,4) = temp_hessmat(4,2);
	//k2a3
	temp_hessmat(5,2) = arma::accu( - arma::exp(h + eta2) % ( dcommonVecda3 % Lambda02 ) );
	temp_hessmat(2,5) = temp_hessmat(5,2);

	//***********k3 section***********// 
	//temp_scorevec(4) = arma::accu( delta1 % delta2 - commonVec % Lambda03 % arma::exp(h+eta3));

	//k3k3
	temp_hessmat(4,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecdk3 % Lambda03diff + commonVec % Lambda03diff ) );
	//k3a1
	temp_hessmat(1,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecda1 % Lambda03diff ) );
	temp_hessmat(4,1) = temp_hessmat(1,4);
	//k3a2
	temp_hessmat(3,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecda2 % Lambda03diff ) );
	temp_hessmat(4,3) = temp_hessmat(3,4);
	//k3a3
	temp_hessmat(5,4) = arma::accu( - arma::exp(h + eta3) % ( dcommonVecda3 % Lambda03diff + exp(a3) * commonVec % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) ) );
	temp_hessmat(4,5) = temp_hessmat(5,4);

	//***********a1 section**********// 
	//temp_scorevec(1) = arma::accu( delta1 % (1 + exp(a1) * arma::log(y1)) - exp(a1) * commonVec % Lambda01 % arma::exp(h+eta1) % arma::log(y1));
	//rewritten as 		 arma::accu( delta1 + exp(a1) * delta1 % arma::log(y1) - arma::exp(h+eta1) % arma::log(y1) % (exp(a1) * commonVec % Lambda01) ); 

	//a1a1
	temp_hessmat(1,1) = arma::accu( exp(a1) * delta1 % arma::log(y1) - 
									arma::exp(h+eta1) % arma::log(y1) % ( exp(a1) * commonVec % Lambda01 + 
																		  exp(a1) * dcommonVecda1 % Lambda01 + 
																		  exp(2 * a1) * commonVec % Lambda01 % arma::log(y1) ) ); 
	//a1a2
	temp_hessmat(1,3) = arma::accu( - arma::exp(h+eta1) % arma::log(y1) % (exp(a1) * dcommonVecda2 % Lambda01) ); 
	temp_hessmat(3,1) = temp_hessmat(1,3);
	//a1a3
	temp_hessmat(1,5) = arma::accu( - arma::exp(h+eta1) % arma::log(y1) % (exp(a1) * dcommonVecda3 % Lambda01) ); 
	temp_hessmat(5,1) = temp_hessmat(1,5);

	//***********a2 section**********// 
	//temp_scorevec(3) = arma::accu( (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) - commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1) );

	//rewritten as arma::accu( (1-delta1) % delta2 + exp(a2) * (1-delta1) % delta2 % arma::log(y1) - arma::exp(h+eta2) % arma::log(y1) % (exp(a2) * commonVec % Lambda02) );
	//a2a2
	temp_hessmat(3,3) = arma::accu( exp(a2) * (1-delta1) % delta2 % arma::log(y1) - arma::exp(h+eta2) % arma::log(y1) % (exp(a2) * commonVec % Lambda02 + exp(a2) * dcommonVecda2 % Lambda02 + exp(2 * a2) * commonVec % Lambda02 % arma::log(y1)) );
	//a2a3
	temp_hessmat(3,5) = arma::accu( - arma::exp(h+eta2) % arma::log(y1) % (exp(a2) * dcommonVecda3 % Lambda02) ); 
	temp_hessmat(5,3) = temp_hessmat(3,5);

	//***********a3 section**********// 
	//temp_scorevec(5) = arma::accu( delta1 % delta2 % (1 + exp(a3) * arma::log(y2)) - 
	//								commonVec % arma::exp(h+a3+eta3) % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) );

	//a3a3
	temp_hessmat(5,5) = arma::accu( exp(a3) * delta1 % delta2 % arma::log(y2) 
									- arma::exp(h+eta3) % ( exp(a3) * commonVec % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) 
														    + exp(a3) * dcommonVecda3 % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) 
														    + exp(2 * a3) * commonVec % (Lambda03y2 % arma::log(y2) % arma::log(y2) - Lambda03y1 % arma::log(y1) % arma::log(y1)) ) ) ;

	//***********h section***********//
	//temp_scorevec(6) = arma::accu(exp(h)*(delta1 % delta2 / (1+exp(h)) + arma::log(1+exp(h) * AVec)/exp(2*h) - commonVec % AVec) );

	//hh
	// from wolfram alpha (as an alternate derivation, with delta1=a and delta2=b): 
	// https://www.wolframalpha.com/input/?i=derivative+of+e%5Eh+*+(+a*b+%2F+(1%2Be%5Eh)+%2B+log(1+%2B+e%5Eh+*+A)+%2F+(e%5E(2*h))+-+(e%5E(-h)+%2B+a+%2B+b)+*+A+%2F+(1+%2B+e%5Eh+*+A))+with+respect+to+h
	temp_hessmat(6,6) = arma::accu( (-delta1 + 2 * AVec - delta2) / (AVec * exp(h)+1) 
								  + (delta1 - AVec + delta2) / (arma::pow(AVec * exp(h) + 1 , 2))
								  + delta1 % delta2 / (exp(h) + 1)
								  - delta1 % delta2 / (pow(exp(h) + 1, 2))
								  - exp(-h) * arma::log1p(AVec * exp(h))
								  );	
	//ha1
	temp_hessmat(6,1) = arma::accu( dAVecda1 / (1 + exp(h)*AVec) - 
							exp(h) * (AVec % dcommonVecda1 + commonVec % dAVecda1 ));
	temp_hessmat(1,6) = temp_hessmat(6,1);
	//hk2
	temp_hessmat(6,2) = arma::accu( dAVecdk2 / (1 + exp(h)*AVec) - 
							exp(h) * (AVec % dcommonVecdk2 + commonVec % dAVecdk2 ));
	temp_hessmat(2,6) = temp_hessmat(6,2);
	//ha2
	temp_hessmat(6,3) = arma::accu( dAVecda2 / (1 + exp(h)*AVec) - 
							exp(h) * (AVec % dcommonVecda2 + commonVec % dAVecda2 ));
	temp_hessmat(3,6) = temp_hessmat(6,3);
	//hk3
	temp_hessmat(6,4) = arma::accu( dAVecdk3 / (1 + exp(h)*AVec) - 
							exp(h) * (AVec % dcommonVecdk3 + commonVec % dAVecdk3 ));
	temp_hessmat(4,6) = temp_hessmat(6,4);
	//ha3
	temp_hessmat(6,5) = arma::accu( dAVecda3 / (1 + exp(h)*AVec) - 
							exp(h) * (AVec % dcommonVecda3 + commonVec % dAVecda3 ));
	temp_hessmat(5,6) = temp_hessmat(6,5);

	//***********beta1 section******//
	//temp_scorevec(arma::span(7, 7 + p1 - 1)) = X1.t() * (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	if(p1 > 0){
		//beta1k1
		temp_hessmat(arma::span(7,7 + p1 - 1),0) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecdk1 % Lambda01 + commonVec % Lambda01 ) );
	    temp_hessmat(0,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),0).t();
		//beta1a1
		temp_hessmat(arma::span(7,7 + p1 - 1),1) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecda1 % Lambda01 + exp(a1) * commonVec % Lambda01 % arma::log(y1) ) );
	    temp_hessmat(1,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),1).t();
		//beta1k2
		temp_hessmat(arma::span(7,7 + p1 - 1),2) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecdk2 % Lambda01 ) );
	    temp_hessmat(2,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),2).t();
		//beta1a2
		temp_hessmat(arma::span(7,7 + p1 - 1),3) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecda2 % Lambda01 ) );
	    temp_hessmat(3,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),3).t();
		//beta1k3
		temp_hessmat(arma::span(7,7 + p1 - 1),4) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecdk3 % Lambda01 ) );
	    temp_hessmat(4,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),4).t();
		//beta1a3
		temp_hessmat(arma::span(7,7 + p1 - 1),5) = X1.t() * ( - arma::exp(h + eta1) % ( dcommonVecda3 % Lambda01 ) );
	    temp_hessmat(5,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),5).t();
	    //beta1h
	    temp_hessmat(arma::span(7,7 + p1 - 1),6) = X1.t() * ( - Lambda01 % arma::exp(h + eta1) % (commonVec - (commonVec % AVec * exp(h) + exp(-h) ) / (1+exp(h)*AVec) ) );
	    temp_hessmat(6,arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),6).t();

	    //beta1beta1
	    temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7,7 + p1 - 1)) = X1.t() * (X1.each_col() % ( - (Lambda01 % arma::exp(h + eta1) % ( commonVec + dcommonVecdk1 )))); //computes sum of w_i * x_ix_i^T

	    if(p2 > 0){
		    //beta1beta2
		    //from univariate case: temp_hess(arma::span(2,p+1),arma::span(2,p+1)) = -( X.t() * (X.each_col() % (exp(k)*arma::exp(eta) % arma::pow(y,exp(a))))); //computes sum of w_i * x_ix_i^T
		    temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)) = X1.t() * (X2.each_col() % ( - dcommonVecdk2 % Lambda01 % arma::exp(h + eta1))); //computes sum of w_i * x_ix_i^T
			temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)).t();
		}

	    if(p3 > 0){
		//beta1beta3
	    //from univariate case: temp_hess(arma::span(2,p+1),arma::span(2,p+1)) = -( X.t() * (X.each_col() % (exp(k)*arma::exp(eta) % arma::pow(y,exp(a))))); //computes sum of w_i * x_ix_i^T
	    temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = X1.t() * (X3.each_col() % ( - dcommonVecdk3 % Lambda01 % arma::exp(h + eta1))); //computes sum of w_i * x_ix_i^T
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),arma::span(7,7 + p1 - 1)) = temp_hessmat(arma::span(7,7 + p1 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)).t();
		}
	}

	//***********beta2 section******//
    //temp_scorevec(arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2.t() * ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
	if(p2 > 0){
		//beta2k1
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),0) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecdk1 % Lambda02 ) );
	    temp_hessmat(0,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),0).t();
		//beta2a1
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),1) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecda1 % Lambda02 ) );
	    temp_hessmat(1,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),1).t();
		//beta2k2
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),2) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecdk2 % Lambda02 + commonVec % Lambda02) );
	    temp_hessmat(2,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),2).t();
		//beta2a2
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),3) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecda2 % Lambda02 + exp(a2) * commonVec % Lambda02 % arma::log(y1) ) );
	    temp_hessmat(3,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),3).t();
		//beta2k3
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),4) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecdk3 % Lambda02 ) );
	    temp_hessmat(4,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),4).t();
		//beta2a3
		temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),5) = X2.t() * ( - arma::exp(h + eta2) % ( dcommonVecda3 % Lambda02 ) );
	    temp_hessmat(5,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),5).t();
	    //beta2h
	    temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),6) = X2.t() * ( - Lambda02 % arma::exp(h + eta2) % (commonVec - (commonVec % AVec * exp(h) + exp(-h) ) / (1+exp(h)*AVec) ) );
	    temp_hessmat(6,arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),6).t();

	    //beta2beta2
	    temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)) = X2.t() * (X2.each_col() % ( - (Lambda02 % arma::exp(h + eta2) % ( commonVec + dcommonVecdk2 )))); //computes sum of w_i * x_ix_i^T

		if(p3 > 0){
		    //beta2beta3
		    //from univariate case: temp_hess(arma::span(2,p+1),arma::span(2,p+1)) = -( X.t() * (X.each_col() % (exp(k)*arma::exp(eta) % arma::pow(y,exp(a))))); //computes sum of w_i * x_ix_i^T
		    temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = X2.t() * (X3.each_col() % ( - dcommonVecdk3 % Lambda02 % arma::exp(h + eta2))); //computes sum of w_i * x_ix_i^T
			temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),arma::span(7 + p1,7 + p1 + p2 - 1)) = temp_hessmat(arma::span(7 + p1,7 + p1 + p2 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)).t();
		}
	}

	//***********beta3 section******//
	//temp_scorevec(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % (Lambda03y2 - Lambda03y1)  % arma::exp(h + eta3));
    //temp_scorevec(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)) = X3.t() * (delta1 % delta2 - commonVec % Lambda03 % arma::exp(h + eta3));
	if(p3 > 0){
		//beta3k1
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),0) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecdk1 % Lambda03diff ) );
	    temp_hessmat(0,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),0).t();
		//beta3a1
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),1) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecda1 % Lambda03diff ) );
	    temp_hessmat(1,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),1).t();
		//beta3k2
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),2) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecdk2 % Lambda03diff ) );
	    temp_hessmat(2,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),2).t();
		//beta3a2
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),3) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecda2 % Lambda03diff ) );
	    temp_hessmat(3,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),3).t();
		//beta3k3
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),4) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecdk3 % Lambda03diff  + commonVec % Lambda03diff ) );
	    temp_hessmat(4,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),4).t();
		//beta3a3
		temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),5) = X3.t() * ( - arma::exp(h + eta3) % ( dcommonVecda3 % Lambda03diff + exp(a3) * commonVec % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) ) );
	    temp_hessmat(5,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),5).t();
	    //beta3h
	    temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),6) = X3.t() * ( - Lambda03diff % arma::exp(h + eta3) % (commonVec - (commonVec % AVec * exp(h) + exp(-h) ) / (1+exp(h)*AVec) ) );
	    temp_hessmat(6,arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),6).t();    

	    //beta3beta3
	    temp_hessmat(arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1),arma::span(7 + p1 + p2,7 + p1 + p2 + p3 - 1)) = X3.t() * (X3.each_col() % ( - (Lambda03diff % arma::exp(h + eta3) % ( commonVec + dcommonVecdk3 )))); //computes sum of w_i * x_ix_i^T
	}
	return -temp_hessmat;
}





// [[Rcpp::export]] 
arma::mat ngradWB_ID_frail_mat_SM(const arma::vec& para,
								  const arma::vec& y1,const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2, 
								  const arma::mat& X1, const arma::mat& X2, const arma::mat& X3){
	//define constants
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

	double k1 = para(0);
	double a1 = para(1);
	double k2 = para(2);
	double a2 = para(3);
	double k3 = para(4);
	double a3 = para(5);
	double h = para(6);

	//define linear predictors
	arma::vec eta1, eta2, eta3;
	if(p1 == 0){
		eta1 = arma::zeros(n);
	} else{
		eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
	}
	if(p2 == 0){
		eta2 = arma::zeros(n);
	} else{
		eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
	}
	if(p3 == 0){
		eta3 = arma::zeros(n);
	} else{
		eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
	}

	arma::mat temp_scoremat = arma::zeros<arma::mat>(n,p1+p2+p3+7);

	arma::vec AVec = getAVecWB_SM(y1, y2, k1, a1, k2, a2, k3, a3, eta1, eta2, eta3);
	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);
	arma::vec Lambda01 = getLambda0gWB(y1, a1, k1);
	arma::vec Lambda02 = getLambda0gWB(y1, a2, k2);
	arma::vec Lambda03 = getLambda0gWB(y2-y1, a3, k3);

	arma::vec logdiff = arma::log(y2-y1);
	logdiff = logdiff.replace(-arma::datum::inf, 0); //log(y2-y1) with the negative infinity values replaced with 0's.


	//k1 (what ina calls u5)
	temp_scoremat(arma::span(0,n-1),0) = delta1 - commonVec % Lambda01 % arma::exp(h + eta1);

	//a1 (what ina calls u6)
	temp_scoremat(arma::span(0,n-1),1) = delta1 % (1 + exp(a1) * arma::log(y1)) -
									exp(a1) * commonVec % Lambda01 % arma::exp(h+eta1) % arma::log(y1);
	//k2 (what ina calls u7)
	temp_scoremat(arma::span(0,n-1),2) = (1-delta1) % delta2 - 
									commonVec % Lambda02 % arma::exp(h+eta2);
	//a2 (what ina calls u8)
	temp_scoremat(arma::span(0,n-1),3) = (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) - 
									commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1);
	//k3 (what ina calls u9)
	temp_scoremat(arma::span(0,n-1),4) = delta1 % delta2 - 
									commonVec % Lambda03 % arma::exp(h+eta3);

	//a3 (what ina calls u10)
	temp_scoremat(arma::span(0,n-1),5) = delta1 % delta2 % (1 + exp(a3) * logdiff) - 
									commonVec % Lambda03 % arma::exp(h+a3+eta3) % logdiff;
	//h (what ina calls u1)
	temp_scoremat(arma::span(0,n-1),6) = exp(h)*(delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec);

	//beta1 (what ina calls u2)
	if(p1 > 0){
		temp_scoremat(arma::span(0,n-1),arma::span(7, 7 + p1 - 1)) = X1;
		temp_scoremat(arma::span(0,n-1),arma::span(7, 7 + p1 - 1)).each_col() %= (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));
	}

	//beta2 (what ina calls u3)
	if(p2 > 0){
		temp_scoremat(arma::span(0,n-1),arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2;
		temp_scoremat(arma::span(0,n-1),arma::span(7 + p1, 7 + p1 + p2 - 1)).each_col() %= ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
	}

	//beta3 (what ina calls u4)
	if(p3 > 0){
		temp_scoremat(arma::span(0,n-1),arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1))  = X3;
		temp_scoremat(arma::span(0,n-1),arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)).each_col() %= (delta1 % delta2 - commonVec % Lambda03 % arma::exp(h + eta3));
	}

    return -temp_scoremat;
}


// [[Rcpp::export]] 
arma::mat ngradWB_ID_frail_mat_M(const arma::vec& para,
								 const arma::vec& y1,const arma::vec& y2, const arma::vec& delta1, const arma::vec& delta2, 
								 const arma::mat& X1, const arma::mat& X2, const arma::mat& X3){
	//define constants
	int p1 = X1.n_cols; int p2 = X2.n_cols;	int p3 = X3.n_cols;	int n = X1.n_rows;

	double k1 = para(0);
	double a1 = para(1);
	double k2 = para(2);
	double a2 = para(3);
	double k3 = para(4);
	double a3 = para(5);
	double h = para(6);

	//define linear predictors
	arma::vec eta1, eta2, eta3;
	if(p1 == 0){
		eta1 = arma::zeros(n);
	} else{
		eta1 = X1 * para(arma::span(7, 7 + p1 - 1));
	}
	if(p2 == 0){
		eta2 = arma::zeros(n);
	} else{
		eta2 = X2 * para(arma::span(7 + p1, 7 + p1 + p2 - 1));
	}
	if(p3 == 0){
		eta3 = arma::zeros(n);
	} else{
		eta3 = X3 * para(arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1));
	}

	arma::mat temp_scoremat = arma::zeros<arma::mat>(n,p1+p2+p3+7);

	arma::vec AVec = getAVecWB_SM(y1, y2, k1, a1, k2, a2, k3, a3, eta1, eta2, eta3);
	arma::vec commonVec = getCommonVec(delta1, delta2, AVec, h);
	arma::vec Lambda01 = getLambda0gWB(y1, a1, k1);
	arma::vec Lambda02 = getLambda0gWB(y1, a2, k2);
	arma::vec Lambda03y2 = getLambda0gWB(y2, a3, k3);
	arma::vec Lambda03y1 = getLambda0gWB(y1, a3, k3);


	//k1 (what ina calls u5)
	temp_scoremat(arma::span(0,n-1),0) = delta1 - commonVec % Lambda01 % arma::exp(h + eta1);

	//a1 (what ina calls u6)
	temp_scoremat(arma::span(0,n-1),1) = delta1 % (1 + exp(a1) * arma::log(y1)) -
									exp(a1) * commonVec % Lambda01 % arma::exp(h+eta1) % arma::log(y1);
	//k2 (what ina calls u7)
	temp_scoremat(arma::span(0,n-1),2) = (1-delta1) % delta2 - 
									commonVec % Lambda02 % arma::exp(h+eta2);
	//a2 (what ina calls u8)
	temp_scoremat(arma::span(0,n-1),3) = (1-delta1) % delta2 % (1 + exp(a2) * arma::log(y1)) - 
									commonVec % Lambda02 % arma::exp(h+a2+eta2) % arma::log(y1);
	//k3 (what ina calls u9)
	temp_scoremat(arma::span(0,n-1),4) = delta1 % delta2 - 
									commonVec % (Lambda03y2 - Lambda03y1) % arma::exp(h+eta3);

	//a3 (what ina calls u10)
	temp_scoremat(arma::span(0,n-1),5) = delta1 % delta2 % (1 + exp(a3) * arma::log(y2)) - 
									commonVec % arma::exp(h+a3+eta3) % (Lambda03y2 % arma::log(y2) - Lambda03y1 % arma::log(y1)) ;
	//h (what ina calls u1)
	temp_scoremat(arma::span(0,n-1),6) = exp(h)*(delta1 % delta2/(1+exp(h)) + arma::log1p(exp(h) * AVec)/exp(2*h) - commonVec % AVec);

	//beta1 (what ina calls u2)
	if(p1 > 0){
		temp_scoremat(arma::span(0,n-1),arma::span(7, 7 + p1 - 1)) = X1;
		temp_scoremat(arma::span(0,n-1),arma::span(7, 7 + p1 - 1)).each_col() %= (delta1 - commonVec % Lambda01 % arma::exp(h + eta1));		
	}

	//beta2 (what ina calls u3)
	if(p2 > 0){
		temp_scoremat(arma::span(0,n-1),arma::span(7 + p1, 7 + p1 + p2 - 1)) = X2;
		temp_scoremat(arma::span(0,n-1),arma::span(7 + p1, 7 + p1 + p2 - 1)).each_col() %= ( (1-delta1) % delta2 - commonVec % Lambda02 % arma::exp(h + eta2)) ;
	}

	//beta3 (what ina calls u4)
	if(p3 > 0){
		temp_scoremat(arma::span(0,n-1),arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1))  = X3;
		temp_scoremat(arma::span(0,n-1),arma::span(7 + p1 + p2 , 7 + p1 + p2 + p3 - 1)).each_col() %= (delta1 % delta2 - commonVec % (Lambda03y2 - Lambda03y1) % arma::exp(h + eta3));
	}
	
    return -temp_scoremat;
}







