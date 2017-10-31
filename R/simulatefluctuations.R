#' Simulate time series of fluctuations
#'
#' Simulate a time series of fluctuations based on specified underlying process. Options are AR(1), ARMA(1,1), P-splines or GPR.
#'
#' @param nyears number of years of observation period
#' @param prop.sample proportion of total years to 'observe'. Default is 1.
#' @param method string to indicate type of method from which to simulate. Options are: "ar", "arma", "splines", "gp".
#' @param params list of parameters needed to simulate from method specified. Each method takes different parameter arguments. See documentation for each process for more information.
#' \itemize{
#' \item{AR: specify nyears, rho, sigma}
#' \item{ARMA: specify nyears, phi, theta, sigma.ar}
#' \item{Splines: specify x.i, degree}
#' \item{GP: specify covariance function (sqexp or matern). For sqexp specify nyears, tau, l. For matern specify nyears, tau, smoothness and range.}
#' }
#' @param obs.err whether or not to include observational error. Default is TRUE.
#' @param sigma.y value of sd of observational error, if included
#' @param seed value of random seed
#' @export
#' @return A data frame of observation times and y-values
#' @seealso \code{\link{GetAR}}, \code{\link{GetARMA}}, \code{\link{GetPSplines}}, \code{\link{GetGP}}
#' @examples
#' nyears <- 100
#' prop.sample <- 0.7
#' obs.err <- TRUE
#' sigma.y <- 0.5
#' seed <- 123
#' method <- "ar"
#' params <- list(rho = 0.5, sigma = 1, ystart = NULL, eps0.t = NULL)
#' res <- simulateFluctuations(nyears, prop.sample, method, params, obs.err, sigma.y)
#' ggplot(data = res, aes(x = t, y = y)) + geom_line() + geom_point() + theme_bw()

simulateFluctuations <- function(
  nyears, ##<< number of years of observation period
  prop.sample = 1, ##<< proportion of total years to 'observe'
  method, ##<< string to indicate type of method from which to simulate. options: "ar", "arma", "splines", "gp"
  params, ##<< list of parameters needed to simulate from method specified
  obs.err = T, ##<< whether or not to include observational error
  sigma.y = NULL, ##<< value of sd of observational error, if included
  seed = 123 ##<< set seed for random draws
){
  # make parameter list into variables
  for (i in 1:length(params)){
    assign(names(params)[i], params[[i]])
  }

  ## AR(1)
  if(method=="ar"){
    mu.t <- GetAR(nyears, rho, sigma, eps0.t, ystart, seed)
    set.seed(seed+1)
    if(obs.err){
      y.t <- mu.t + rnorm(nyears, 0, sigma.y)
    }
    else{
      y.t <- mu.t
    }
    t.i <- sort(sample(seq(1,nyears), size = round(nyears*prop.sample)))
    y.i <- y.t[t.i]
    d <- data.frame(t = t.i, y = y.i)
  }

  ## ARMA(1,1)
  if(method=="arma"){
    mu.t <- GetARMA(nyears, phi, theta, sigma.ar, ystart, seed)
    set.seed(seed+2)
    if(obs.err){
      y.t <- mu.t + rnorm(nyears, 0, sigma.y)
    }
    else{
      y.t <- mu.t
    }
    t.i <- sort(sample(seq(1,nyears), size = round(nyears*prop.sample)))
    y.i <- y.t[t.i]
    d <- data.frame(t = t.i, y = y.i)
  }

  ## Psplines
  if(method=="splines"){
    x.t <- 1:nyears
    res <- GetSplines(x.t)
    K <- length(res$knots.k)
    mu.t <- GetPSplines(res$B.ik, sigma.alpha, order, seed)
    set.seed(seed+3)
    if(obs.err){
      y.t <- mu.t + rnorm(nyears, 0, sigma.y)
    }
    else{
      y.t <- mu.t
    }
    t.i <- sort(sample(seq(1,nyears), size = round(nyears*prop.sample)))
    y.i <- y.t[t.i]
    d <- data.frame(t = t.i, y = y.i)
  }

  ## GP
  if(method=="gp"){
    if(cov.method=="sqexp") res <- GetGP(nyears, cov.method = "sqexp", l = l, tau = tau, seed = seed)
    if(cov.method=="matern") res <- GetGP(nyears, cov.method = "matern", smoothness = smoothness, range = range, tau=tau, seed = seed)
    x.t <- 1:nyears
    mu.t <- res$y
    set.seed(seed+4)
    if(obs.err){
      y.t <- mu.t + rnorm(nyears, 0, sigma.y)
    }
    else{
      y.t <- mu.t
    }
    t.i <- sort(sample(seq(1,nyears), size = round(nyears*prop.sample)))
    y.i <- y.t[t.i]
    d <- data.frame(t = t.i, y = y.i)
  }
  return(d)
}

