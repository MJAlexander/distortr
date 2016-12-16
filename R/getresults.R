#' Get results from MCMC estimation
#'
#' @param mod A JAGS model object
#' @param method The method of smoothing to implement (choices: ar, arma, splines, gp)
#' @param alpha.level Significance level. Default is 5\%.
#' @export
#' @return A data frame of x values, estimates and uncertainty intervals.
#' @seealso \code{\link{runMCMC}, \link{plotResults}}

getResults <- function(mod,
                       method,
                       alpha.level = 0.05){
  # get out the mu estimates
  if(method %in% c("ar", "arma", "splines")){
    mu.st <- mod$BUGSoutput$sims.list[["mu.t"]]
  }
  if(method=="gp"){
    mu.st <- mod$BUGSoutput$sims.list[["G"]]
  }

  # get df of results
  mu.qt <- apply(mu.st, 2, quantile, c(alpha.level/2, 0.5, (1-alpha.level/2)))
  df.mu <- data.frame(time = 1:nyears, t(mu.qt))
  colnames(df.mu) <- c("time", "lower", "median", "upper")
  return(df.mu)
}
