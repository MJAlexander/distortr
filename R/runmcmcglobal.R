#' Run MCMC estimation for hierarichal models with multiple countries
#'
#' Run MCMC estimation of time series data for multiple countries using JAGS.
#'
#' @param method The method of smoothing to implement (choices: ar, arma, splines, gp)
#' @param input.data List of required input data. See \code{processData}, \code{getSplinesData} and \code{getGPData} to get required data in compatible form.
#' @param order The order of splines penalization (either 1 or 2)
#' @param matern.cov Whether or not to use Matern covariance function. Default is \code{TRUE}.
#' @param cs.arma whether ARMA parameter(s) are country specific. If `FALSE`, parameter is global.
#' @param cs.smoothing whether smoothing paramter is country specific. If `FALSE`, smoothing parameter is global.
#' @param time.trend if `TRUE` a linear time trend is estimated.
#' @param nserror.estimated whether to estimate non-sampling error. IF `FALSE`, fixed sampling error is inputted.
#' @param nchains Number of MCMC chains
#' @param nburnin Number of iterations to throw away as burn in.
#' @param niter Number of total iterations.
#' @param nthin Degree of thinning of MCMC chains
#' @param model.file.path Text file which contains the model to be fitted. If \code{NULL}, the text file is drawn from the \code{models} folder.
#' @param model.save.file.path Path to save model, if written.
#' @export
#' @return A JAGS model object
#' @seealso \code{\link{processData}, \link{getSplinesData}, \link{getGPData}}
#' @examples
#' nyears <- 100
#' prop.sample <- 0.7
#' obs.err <- T
#' sigma.y <- 0.5
#' seed <- 123
#' method <- 'splines'
#' params <- list(sigma.alpha = 1, order = 1)
#' df <- simulateFluctuations(nyears, prop.sample, method, params, obs.err, sigma.y)
#' df$se <- 1
#' mod <- runMCMC(df = df, nyears = 100, method = "splines", order = 1,
#' nchains = 4, nburnin = 100, niter = 100+3000, nthin = 3)

runMCMCGlobal <- function(method,
                          input.data,
                          order = NULL,
                          matern.cov=TRUE,
                          cs.arma = NULL,
                          cs.smoothing = TRUE,
                          time.trend = FALSE,
                          nserror.estimated = TRUE,
                          nchains = 3,
                          nburnin = 1000,
                          niter = 2000,
                          nthin = 1,
                          model.file.path = NULL,
                          model.save.file.path = "R/model.txt"){

  # extract all data from list
  for(i in 1:length(input.data)){
    assign(names(input.data)[i], input.data[[i]])
  }
  # model path to run
  if(is.null(model.file.path)){
    model.path.to.run <- model.save.file.path
  }
  else{
    model.path.to.run <- model.file.path
  }

  # base data. will be added to.
  jags.data <- list(y.ci = y.ci, gett.ci = (gett.ci - startyear.c+1), niso = niso, n.c =n.c,
                    nyears.c= nyears.c, se.ci = se.ci, sigma.y = sigma.y,
                    source.ci = source.ci,
                    region.c = region.c, nregions = nregions)

  if(nserror.estimated){
    jags.data$nsources <- nsources
  }
  if(method=="ar"){
    if(is.null(cs.arma)){
      stop("Need to specify cs.arma.")
    }
    if(is.null(model.file.path)){
      # write model based on desired parameters
      writeModelAR(cs.arma = cs.arma,
                   cs.smoothing = cs.smoothing,
                   time.trend = time.trend,
                   nserror.estimated = nserror.estimated,
                   file.name = model.save.file.path)
    }

    parnames <- c("mu.ct", "loglike.ci", "yrep.ci", "beta", "sigma", "rho", "sigma.y",
                  "mu.beta", "sigma.beta", "mu.beta.global", "sigma.beta.global")

    if(time.trend){
      parnames <- c(parnames, "gamma", "mu.gamma", "sigma.gamma", "mu.gamma.global", "sigma.gamma.global")
    }
    if(cs.arma){
      parnames <- c(parnames, "mu.rho", "sigma.rho", "mu.rho.global", "sigma.rho.global")
    }
    if(cs.smoothing){
      parnames <- c(parnames, "mu.logsigma", "sigma.logsigma", "mu.logsigma.global", "sigma.logsigma.global")
    }
  }

  if(method=="arma"){
    if(is.null(cs.arma)){
      stop("Need to specify cs.arma.")
    }
    if(is.null(model.file.path)){
      # write model based on desired parameters
      writeModelARMA(cs.arma = cs.arma,
                     cs.smoothing = cs.smoothing,
                     time.trend = time.trend,
                     nserror.estimated = nserror.estimated,
                     file.name = model.save.file.path)
    }

    parnames <- c("mu.ct", "loglike.ci", "yrep.ci", "beta", "eta", "rho", "theta", "sigma.y",
                  "mu.beta", "sigma.beta", "mu.beta.global", "sigma.beta.global")

    if(time.trend){
      parnames <- c(parnames, "gamma", "mu.gamma", "sigma.gamma", "mu.gamma.global", "sigma.gamma.global")
    }
    if(cs.arma){
      parnames <- c(parnames, "mu.rho", "sigma.rho", "mu.rho.global", "sigma.rho.global",
                    "mu.theta", "sigma.theta", "mu.theta.global", "sigma.theta.global")
    }
    if(cs.smoothing){
      parnames <- c(parnames, "mu.eta", "sigma.eta", "mu.eta.global", "sigma.eta.global")
    }
  }

  if(method=="splines"){
    if(is.null(order)){
      stop("Order of penalization must be specified.")
    }
    if(is.null(model.file.path)){
      # write model based on desired parameters
      writeModelSplines(order = order,
                        cs.smoothing = cs.smoothing,
                        nserror.estimated = nserror.estimated,
                        file.name = model.save.file.path)
    }
    if(order==1){
      jags.data$Z.tkc <- Z.tkc
      jags.data$H.c <- H.c
      jags.data$H <- H
    }
    if(order==2){
      jags.data$Z.tkc <- Z.tkc
      jags.data$BG.tdc <- BG.tdc
      jags.data$H.c <- H.c
      jags.data$H <- H
      jags.data$D <- 2
    }
    parnames <- c("mu.ct",  "loglike.ci", "yrep.ci", "beta.d", "sigma.delta", "sigma.y",
                  "mu.beta", "sigma.beta",
                  "mu.beta.global", "sigma.beta.global")
    if(cs.smoothing){
      parnames <- c(parnames, "chi.delta", "psi.delta", "mu.chi.delta", "sigma.chi.delta")
    }
  }

  if(method=="gp"){
    if(is.null(model.file.path)){
      # write model based on desired parameters
      if(matern.cov==TRUE) cov.fun <- "matern"
      if(matern.cov==FALSE) cov.fun <- "sqexp"
      writeModelGP(cov.fun = cov.fun,
                        time.trend = time.trend,
                        cs.smoothing = cs.smoothing,
                        nserror.estimated = nserror.estimated,
                        file.name = model.save.file.path)
    }

    if(matern.cov==TRUE){
      jags.data$Sigma.corr <- Sigma.corr.c
    }
    if(matern.cov==FALSE){
      jags.data$Dist <- Dist.c
      jags.data$kappa <- 2
    }

    parnames <- c("G", "loglike.ci", "yrep.ci", "beta", "sigma.g", "sigma.y",
                  "mu.beta", "sigma.beta", "mu.beta.global", "sigma.beta.global")

    if(time.trend){
      parnames <- c(parnames, "gamma", "mu.gamma", "sigma.gamma", "mu.gamma.global", "sigma.gamma.global")
    }

    if(cs.smoothing){
      parnames <- c(parnames, "chi", "sigma.psi", "chi.global", "sigma.chi.global")
    }

    if(matern.cov==FALSE){
      parnames <- c(parnames, "p", "mu.p", "sigma.p", "mu.p.global", "sigma.p.global")
    }

  }

  # remove values of sigma.y if it is to be estimated
  if(nserror.estimated){
    jags.data[["sigma.y"]] <- NULL
  }
  ## run the model
  cat("Running model.\n")
  mod <- jags(data = jags.data,
              parameters.to.save= c(parnames),
              n.chains = nchains,
              n.burnin = nburnin,
              n.iter = niter,
              n.thin =nthin,
              model.file = model.path.to.run)
  if(max(mod$BUGSoutput$summary[, c("Rhat")])>1.1){
    cat("Something hasn't converged.\n")
    cat(paste("Max Rhat is", max(mod$BUGSoutput$summary[, c("Rhat")]), "\n"))

  }
  return(mod)
}


