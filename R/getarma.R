#' Simulate time series of fluctuations based on ARMA process
#'
#' Simulate a time series of fluctuations based on AR(1)
#'
#' @param nyears number of years of observation period
#' @param rho value of autocorrelation. Must be value between -1 and 1.
#' @param sigma standard deviation of distortions
#' @param eps0.t values of distortions. Default is NULL. If NULL, they are drawn random from N(0,1)
#' @param ystart value of first observation. Default is NULL. If NULL, it is drawn randomly from other parameters
#' @param seed value of random seed
#' @export
#' @return A vector of values.
#' @seealso \code{\link{GetARMA}}
#' @examples
#' y.ar <- GetAR(50, 0.8, 1)
#' ggplot(data = NULL, aes(x = 1:50, y = y.ar)) + geom_line() + geom_point() + theme_bw() + xlab("time")

GetAR <- function(nyears,
                  rho,
                  sigma,
                  eps0.t = NULL,
                  ystart = NULL,
                  seed = 123){
  if (is.null(eps0.t)){
    set.seed(seed)
    eps0.t <- rnorm(nyears, 0, 1)
  }
  y.t <- rep(NA, nyears)
  if (is.null(ystart)){
    y.t[1] <- sigma/sqrt(1-rho^2)*eps0.t[1]
  } else {
    y.t[1] <- ystart
  }
  for (t in 2:nyears){
    y.t[t] <- rho*y.t[t-1] + sigma*eps0.t[t]
  }
  return(y.t)
}

#' Simulate time series of fluctuations based on ARMA process
#'
#' Simulate a time series of fluctuations based on ARMA(1,1)
#' @param nyears number of years of observation period
#' @param phi value of autocorrelation.
#' @param theta strength of moving average.
#' @param sigma.ar standard deviation of distortions
#' @param ystart value of first observation. Default is NULL. If NULL, it is drawn randomly from other parameters
#' @param seed value of random seed
#' @export
#' @return A vector of values.
#' @seealso \code{\link{GetAR}}
#' @examples
#' y.arma <- GetARMA(50, 0.8, 0.2, 1)
#' ggplot(data = NULL, aes(x = 1:50, y = y.arma)) + geom_line() + geom_point() + theme_bw() + xlab("time")
GetARMA <- function(nyears,
                    phi,
                    theta,
                    sigma.ar,
                    ystart = NULL,
                    seed = 123){

  gamma0 <- ((1-2*phi*theta + theta^2)/(1-phi^2))*sigma.ar^2
  sqrtgamma0 <- sqrt(gamma0)
  arma.t <- rep(NA, nyears)
  if (is.null(ystart)){
    arma1 <- rnorm(1, 0, sqrtgamma0)
    arma.t[1] <- arma1 # mu.t
  } else {
    arma.t[1] <- ystart
  }
  # e.t is the innovation
  e.t <- rep(NA, nyears)
  e.t[1] <- rnorm(1, sigma.ar^2/gamma0*arma.t[1], sqrt(sigma.ar^2*(1-sigma.ar^2/gamma0)))
  for (t in 2:nyears){
    e.t[t] <- rnorm(1, 0, sigma.ar)
    arma.t[t] <- phi*arma.t[t-1] - theta*e.t[t-1] + e.t[t]
  }
  return(arma.t)
}

