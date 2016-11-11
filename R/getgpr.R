#' Simulate time series of fluctuations based on a Gaussian Process
#'
#' Simulate time series of fluctuations based on a Gaussian Process, using a squared exponential covariance function.
#'
#' @param nyears number of years of observation period
#' @param tau the amplitude parameter. Determines the magnitude of fluctuations over time.
#' @param l the length scale parameter. Determines the 'wigglyness' of the fluctuations.
#' @param seed value of random seed
#' @export
#' @return A data frame of observation points (x) and y-values
#' @seealso \code{\link{calcSigma}}
#' @examples
#' res <- GetGP(50, tau = 1, l = 1)
#' ggplot(data = res, aes (x = x, y = y)) + geom_line() + geom_point() + theme_bw()

GetGP <- function(nyears,
                  tau = 1,
                  l = 1,
                  seed = 123){
  x.star <- 1:nyears
  k.xsxs <- calcSigma(x.star,x.star, l = l, tau = tau)
  set.seed(seed)
  y <- mvrnorm(1, rep(0, length(x.star)), k.xsxs)
  f <- data.frame(x = x.star, y = y)
  return(f)
}

