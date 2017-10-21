#' Get results from MCMC estimation
#'
#' @param mod A JAGS model object
#' @param method The method of smoothing to implement (choices: ar, arma, splines, gp)
#' @param time.trend For global models, whether a time.trend is included.
#' @param iso.number The index number for country of interest, if model is global.
#' @param nyears The number of observation years for country of interest
#' @param startyear The year of first observation for country of interest, if model is global.
#' @param alpha.level Significance level. Default is 5\%.
#' @export
#' @return A data frame of x values, estimates and uncertainty intervals.
#' @seealso \code{\link{runMCMC}, \link{plotResults}}

getResults <- function(mod,
                       method,
                       time.trend = NULL,
                       iso.number = NULL,
                       nyears = NULL,
                       startyear = NULL,
                       alpha.level = 0.05){
  # get out the mu estimates
  if(method %in% c("ar", "arma", "splines")){
    if(grepl("mu.t", paste(dimnames(mod$BUGSoutput$sims.array)[[3]], collapse = " "))){
      # only one time series
      type <- "sim"
      mu.st <- mod$BUGSoutput$sims.list[["mu.t"]]
    }
    if(grepl("mu.ct", paste(dimnames(mod$BUGSoutput$sims.array)[[3]], collapse = " "))){
      # many countries, so must specify iso number
      type <- "global"
      if(is.null(iso.number)){
        stop("Must specify country index.")
      }
      if(is.null(nyears)){
        stop("Must specify number of years of observations.")
      }
      if(is.null(startyear)){
        stop("Must specify start year.")
      }
      mu.st <- mod$BUGSoutput$sims.array[,,paste0("mu.ct[",iso.number,",",1:nyears,"]")]
      mu.st <- sapply(1:dim(mu.st)[3], function(i) c(mu.st[,,i]))
      if(method %in% c("ar", "arma")){
        if(time.trend){
          gamma.s <- mod$BUGSoutput$sims.list[["gamma"]][,iso.number]
          gamma.st <- sapply(1:nyears, function (i) gamma.s*i)
          mu.st <- mu.st+ gamma.st
        }
      }
    }
  }
  if(method=="gp"){
    if(dimnames(mod$BUGSoutput$sims.array)[[3]][1]=="G[1]"){
      type <- "sim"
      mu.st <- mod$BUGSoutput$sims.list[["G"]]
    }
    if(dimnames(mod$BUGSoutput$sims.array)[[3]][1]=="G[1,1]"){
      type <- "global"
      if(is.null(iso.number)){
        stop("Must specify country index.")
      }
      if(is.null(nyears)){
        stop("Must specify number of years of observations.")
      }
      if(is.null(startyear)){
        stop("Must specify start year.")
      }
      mu.st <- mod$BUGSoutput$sims.array[,,paste0("G[",1:nyears,",",iso.number,"]")]
      mu.st <- sapply(1:dim(mu.st)[3], function(i) c(mu.st[,,i]))
    }
  }

  # get df of results
  mu.qt <- apply(mu.st, 2, quantile, c(alpha.level/2, 0.5, (1-alpha.level/2)))
  if(type =="sim"){
    df.mu <- data.frame(time = 1:nyears, t(mu.qt))
  }
  if(type=="global"){
    df.mu <- data.frame(time = (1:nyears)+(startyear-1), t(mu.qt))
  }
  colnames(df.mu) <- c("time", "lower", "median", "upper")
  return(df.mu)
}
