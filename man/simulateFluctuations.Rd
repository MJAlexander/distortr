% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulatefluctuations.R
\name{simulateFluctuations}
\alias{simulateFluctuations}
\title{Simulate time series of fluctuations}
\usage{
simulateFluctuations(nyears, prop.sample = 1, method, params, obs.err = T,
  sigma.y = NULL, seed = 123)
}
\arguments{
\item{nyears}{number of years of observation period}

\item{prop.sample}{proportion of total years to 'observe'. Default is 1.}

\item{method}{string to indicate type of method from which to simulate. Options are: "ar", "arma", "splines", "gp".}

\item{params}{list of parameters needed to simulate from method specified. Each method takes different parameter arguments. See documentation for each process for more information.
\itemize{
\item{AR: specify nyears, rho, sigma}
\item{ARMA: specify nyears, phi, theta, sigma.ar}
\item{Splines: specify x.i, degree}
\item{GP: specify covariance function (sqexp or matern). For sqexp specify nyears, tau, l. For matern specify nyears, tau, smoothness and range.}
}}

\item{obs.err}{whether or not to include observational error. Default is TRUE.}

\item{sigma.y}{value of sd of observational error, if included}

\item{seed}{value of random seed}
}
\value{
A data frame of observation times and y-values
}
\description{
Simulate a time series of fluctuations based on specified underlying process. Options are AR(1), ARMA(1,1), P-splines or GPR.
}
\examples{
nyears <- 100
prop.sample <- 0.7
obs.err <- TRUE
sigma.y <- 0.5
seed <- 123
method <- "ar"
params <- list(rho = 0.5, sigma = 1, ystart = NULL, eps0.t = NULL)
res <- simulateFluctuations(nyears, prop.sample, method, params, obs.err, sigma.y)
ggplot(data = res, aes(x = t, y = y)) + geom_line() + geom_point() + theme_bw()
}
\seealso{
\code{\link{GetAR}}, \code{\link{GetARMA}}, \code{\link{GetPSplines}}, \code{\link{GetGP}}
}
