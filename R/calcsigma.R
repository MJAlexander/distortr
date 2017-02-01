#' Calculate variance-covariance matrix
#'
#' Calculate variance-covariance matrix between two vectors, based on squared exponential kernal.
#'
#' @param X1 first vector of x values
#' @param X2 second vector of x values
#' @param method either \code{sqexp} (squared exponential) or \code{matern}
#' @param l the length scale parameter. Determines the 'wigglyness' of the fluctuations for squared exponential.
#' @param tau the amplitude parameter. Determines the magnitude of fluctuations over time.
#' @param range parameter for Matern function. Default is 10.
#' @param smoothness parameter for Matern function. Default is 2.
#' @export
#' @return A data frame of observation points (x) and y-values
#' @seealso \code{\link{GetGP}}
#' @examples
#' calcSigma(1:5, 1:5, method = "sqexp", l = 1, tau = 1)

calcSigma <- function(X1,X2,
                      method,
                      range = 10, smoothness = 2,
                      l=1, tau = 1) {
  if(method=="sqexp"){
    Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
    for (i in 1:nrow(Sigma)) {
      for (j in 1:ncol(Sigma)) {
        Sigma[i,j] <- tau^2 * exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
      }
    }
  }
  if(method=="matern"){
    d=rdist(X1,X2)
    # calculate the covariance values
    Sigma=tau^2*Matern(d,smoothness=smoothness,range=range)
  }
  return(Sigma)
}
