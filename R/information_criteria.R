#' Get Information Criterion
#'
#' Function to compute AIC, BIC and GCV.
#'
#' @param nll final negative loglikelihood value (on the sum scale).
#' @param df Degrees of freedom (typically, number of estimated parameters, though may differ under fusion)
#' @param n number of observations (denominator of nll)
#' @param ic string specifying the criterion selected.
#'
#' @return
#' @export
get_ic <- function(nll, df, n=NA, ic){
  switch(tolower(ic),
         aic=2*nll + 2 * df,
         bic=2*nll + log(n) * df,
         gcv=nll/(n^2*(1-df/n)^2))
}
