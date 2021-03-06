% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcsigma.R
\name{calcSigma}
\alias{calcSigma}
\title{Calculate variance-covariance matrix}
\usage{
calcSigma(X1, X2, cov.method, range = 10, smoothness = 2, l = 1,
  tau = 1)
}
\arguments{
\item{X1}{first vector of x values}

\item{X2}{second vector of x values}

\item{cov.method}{either \code{sqexp} (squared exponential) or \code{matern}}

\item{range}{parameter for Matern function. Default is 10.}

\item{smoothness}{parameter for Matern function. Default is 2.}

\item{l}{the length scale parameter. Determines the 'wigglyness' of the fluctuations for squared exponential.}

\item{tau}{the amplitude parameter. Determines the magnitude of fluctuations over time.}
}
\value{
A data frame of observation points (x) and y-values
}
\description{
Calculate variance-covariance matrix between two vectors, based on squared exponential kernal.
}
\examples{
calcSigma(1:5, 1:5, cov.method = "sqexp", l = 1, tau = 1)
}
\seealso{
\code{\link{GetGP}}
}
