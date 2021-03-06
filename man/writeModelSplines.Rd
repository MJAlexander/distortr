% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/writemodelsplines.R
\name{writeModelSplines}
\alias{writeModelSplines}
\title{Write a JAGS model to fit penalized splines regression}
\usage{
writeModelSplines(order = NULL, cs.smoothing = T, nserror.estimated = T,
  file.name = "model.txt")
}
\arguments{
\item{order}{order of penalization (either 1 or 2).}

\item{cs.smoothing}{whether smoothing paramter is country specific. If `FALSE`, smoothing parameter is global.}

\item{nserror.estimated}{whether to estimate non-sampling error. IF `FALSE`, fixed sampling error is inputted.}

\item{file.name}{name of file to be saved. Must be a `.txt` file}
}
\value{
A text file that contains a JAGS model
}
\description{
Write a JAGS model to fit first- and second-order penalized splines regression
}
\examples{
order <- 1
cs.smoothing <- TRUE
nserror.estimated <- TRUE
writeModelSplines(order = order, nserror.estimated = nserror.estimated, cs.smoothing = cs.smoothing)
}
