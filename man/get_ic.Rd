% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ic.R
\name{get_ic}
\alias{get_ic}
\title{Get Information Criterion}
\usage{
get_ic(nll, df, n = NA, ic)
}
\arguments{
\item{nll}{final negative loglikelihood value (on the sum scale).}

\item{df}{Degrees of freedom (typically, number of estimated parameters, though may differ under fusion)}

\item{n}{number of observations (denominator of nll)}

\item{ic}{string specifying the criterion selected.}
}
\value{
a number representing the requested information criterion.
}
\description{
Function to compute AIC, BIC and GCV.
}
