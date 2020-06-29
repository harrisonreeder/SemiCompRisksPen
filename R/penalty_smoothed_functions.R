

#function for the nesterov smoothed approximation of the fused lasso
#Defined in Chen (2012) equation 3.5, with the solution alphastar defined
#in proposition 2
smoothed_fused_lasso_func <- function(para,pen_mat_w_lambda, mu_smooth_fused){
  if(is.null(pen_mat_w_lambda) | is.null(mu_smooth_fused)){
    return(0)
  }
  if(mu_smooth_fused==0){
    return(norm(pen_mat_w_lambda %*% para, type = "1"))
  }
  temp_vec <- pen_mat_w_lambda %*% para
  alphastar <-  sign(temp_vec/mu_smooth_fused)*pmin(1, abs(temp_vec/mu_smooth_fused))

  return( as.vector( t(alphastar) %*% temp_vec - mu_smooth_fused*sum(alphastar^2)/2 ) )
}


#function for the gradient of the nesterov smoothed approximation of the fused lasso
smoothed_fused_lasso_prime_func <- function(para,pen_mat_w_lambda, mu_smooth_fused){
  if(is.null(pen_mat_w_lambda) | is.null(mu_smooth_fused)){
    return(0)
  }
  if(mu_smooth_fused==0){
    return(0)
  }
  temp_vec <- pen_mat_w_lambda %*% para /mu_smooth_fused
  alphastar <-  sign(temp_vec)*pmin(1, abs(temp_vec))

  return( as.vector( t(pen_mat_w_lambda) %*% alphastar ) )
}

