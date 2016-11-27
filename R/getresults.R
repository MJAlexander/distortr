#' Get results from MCMC estimation
#'
#' @param mod A JAGS model object
#' @param method The method of smoothing to implement (choices: ar, arma, splines, gp)
#' @param order The order of splines penalization (either 1 or 2)
#' @param alpha Significance level. Default is 5\%.
#' @export
#' @return A data frame of x values, estimates and uncertainty intervals.
#' @seealso \code{\link{runMCMC, plotResults}}
#' @examples
#' mod <- runMCMC(df = res, nyears = 100, method = "splines", order = 1,
#' nchains = 4, nburnin = 100, niter = 100+3000, nthin = 3)
#' df.mu <- getResults(mod, method = "splines")

getResults <- function(mod,
                       method,
                       order = NULL,
                       alpha = 0.05){
  # get out the mu estimates
  if(method %in% c("ar", "arma", "splines")){
    mu.st <- mod$BUGSoutput$sims.list[["mu.t"]]
  }
  if(method=="gp"){
    mu.st <- mod$BUGSoutput$sims.list[["G"]]
  }

  # get df of results
  mu.qt <- apply(mu.st, 2, quantile, c(alpha/2, 0.5, (1-alpha/2)))
  df.mu <- data.frame(time = 1:nyears, t(mu.qt))
  colnames(df.mu) <- c("time", "lower", "median", "upper")
  return(df.mu)
}
