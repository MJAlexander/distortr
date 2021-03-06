% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotdata.R
\name{plotData}
\alias{plotData}
\title{Plot time series of data}
\usage{
plotData(data.df, maintitle = NULL, plot.se = T)
}
\arguments{
\item{data.df}{A data frame of t values and data values, with standard errors}

\item{maintitle}{Title of plot}

\item{plot.se}{Whether or not to plot standard errors}
}
\value{
A plot of time series data with standard errors, if desired.
}
\description{
Plot time series of data
}
\examples{
nyears <- 100
prop.sample <- 0.7
obs.err <- TRUE
sigma.y <- 0.5
seed <- 123
method <- 'splines'
params <- list(sigma.alpha = 1, order = 1)
res <- simulateFluctuations(nyears, prop.sample, method, params, obs.err, sigma.y)
res$se <- 0.1
plotData(res)
}
