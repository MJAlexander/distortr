#' Get Gaussian Process covariance data to run JAGS country model
#'
#' Function to get Gaussian Process variance-covariance data to run JAGS country model
#'
#' @param nyears.c vector with number of years of observations for each country
#' @param niso number of countries
#' @param cov.method either \code{sqexp} (squared exponential) or \code{matern}
#' @param range parameter for Matern function. Default is 10.
#' @param smoothness parameter for Matern function. Default is 2.
#' @export
#' @return An array which contains either the distance matrix for each country (squared exponential) or the covariance matrix for each country (matern).

getGPData <- function(nyear.c,
                      niso,
                      cov.method = "matern",
                      range = 10, smoothness = 2){

  if(!(cov.method %in% c("matern", "sqexp"))){
    stop("cov.method must either be matern or sqexp")
  }
  max.nyears <- max(nyears.c)
  Dist.c <- array(NA, c(max.nyears, max.nyears, niso))

  for(i in 1:niso){
    Dist.c[1:nyears.c[i],1:nyears.c[i],i] <- rdist(1:nyears.c[i])
  }

  Sigma.corr.c <-array(NA, c(max.nyears, max.nyears, niso))

  for(i in 1:niso){
    Sigma.corr.c[1:nyears.c[i],1:nyears.c[i],i] <- calcSigma(1:nyears.c[i], 1:nyears.c[i],
                                                           cov.method = "matern", range = range, smoothness = smoothness)
  }
  if(cov.method=="matern"){
    return(Sigma.corr.c)
  }
  if(cov.method=="sqexp"){
    return(Dist.c)
  }

}
