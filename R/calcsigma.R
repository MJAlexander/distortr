#' Calculate variance-covariance matrix
#'
#' Calculate variance-covariance matrix between two vectors, based on squared exponential kernal.
#'
#' @param X1 first vector of x values
#' @param X2 second vector of x values
#' @param l the length scale parameter. Determines the 'wigglyness' of the fluctuations.
#' @param tau the amplitude parameter. Determines the magnitude of fluctuations over time.
#' @export
#' @return A data frame of observation points (x) and y-values
#' @seealso \code{\link{GetGP}}
#' @examples
#' calcSigma(1:5, 1:5, l = 1, tau = 1)

calcSigma <- function(X1,X2,l=1, tau = 1) {
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- tau^2 * exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
    }
  }
  return(Sigma)
}
