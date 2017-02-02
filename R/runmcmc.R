#' Run MCMC estimation
#'
#' Run MCMC estimation of time series using JAGS.
#'
#' @param df A dataframe of x and y observations, and standard errors around ys
#' @param nyears number of years of observations
#' @param method The method of smoothing to implement (choices: ar, arma, splines, gp)
#' @param order The order of splines penalization (either 1 or 2)
#' @param matern.cov Whether or not to use Matern covariance function. Default is \code{TRUE}.
#' @param obs.err is TRUE if standard errors are observed
#' @param measurement.err is TRUE if there is assumed to be measurement error
#' @param nchains Number of MCMC chains
#' @param nburnin Number of iterations to throw away as burn in.
#' @param niter Number of total iterations.
#' @param nthin Degree of thinning of MCMC chains
#' @param model.file.path Text file which contains the model to be fitted. If \code{NULL}, the text file is drawn from the \code{models} folder.
#' @export
#' @return A JAGS model object
#' @seealso \code{\link{getResults}, \link{plotResults}}
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

runMCMC <- function(df,
                    nyears,
                    method,
                    order = NULL,
                    matern.cov=TRUE,
                    obs.err = FALSE,
                    measurement.err = FALSE,
                    nchains = 4,
                    nburnin = 1000,
                    niter = 1000+30000,
                    nthin =30,
                    model.file.path = NULL){

  if(obs.err){
    nu.i = 1/(df$se)^2
  }
  else{
    nu.i <- rep(0,nrow(df))
  }
  if(method=="ar"){
    if(is.null(model.file.path)){
      ifelse(measurement.err, model.file.path <- "R/models/model_ar.txt",
             model.file.path <- "R/models/model_ar_nme.txt")
    }
    jags.data <- list(y.i = df$y, gett.i = df$t, nyears=nyears, n = length(df$t))
    parnames <- c("sigma", "rho", "mu.t")
    if(measurement.err==TRUE){
      jags.data$nu.i <- nu.i
      parnames <- c(parnames, "sigma.y")
    }
  }

  if(method=="arma"){
    if(is.null(model.file.path)){
      ifelse(measurement.err, model.file.path <- "R/models/model_arma.txt",
             model.file.path <- "R/models/model_arma_nme.txt")
    }
    jags.data <- list(y.i = df$y, gett.i = df$t,nyears=nyears, n = length(df$t))
    parnames <- c("sigma.ar", "phi", "theta", "mu.t")
    if(measurement.err==TRUE){
      jags.data$nu.i <- nu.i
      parnames <- c(parnames, "sigma.y")
    }
  }

  if(method=="splines"){
    if(is.null(order)){
      stop("Order of penalization must be specified.")
    }
    if(is.null(model.file.path)){
      ifelse(measurement.err, model.file.path <- paste0("R/models/model_splines_", order,".txt"),
             model.file.path <- model.file.path <- paste0("R/models/model_splines_", order,"_nme.txt"))
    }
    x.t <- 1:nyears
    sp <- GetSplines(x.t)
    K <- length(sp$knots.k)
    B.tk <- sp$B.ik
    jags.data <- list(y.i = df$y, gett.i = df$t, nyears=nyears, n = length(df$t), K = K, B.tk = B.tk)
    parnames <- c("alpha.k", "sigma.alpha", "mu.t")
    if(measurement.err==TRUE){
      jags.data$nu.i <- nu.i
      parnames <- c(parnames, "sigma.y")
    }
  }

  if(method=="gp"){
    if(matern.cov==TRUE){
      if(is.null(model.file.path)){
        ifelse(measurement.err, model.file.path <- "R/models/model_gp_matern.txt",
               model.file.path <- "R/models/model_gp_matern_nme.txt")
      }
      ## can just calculate Sigma up to the amplitude here because in CODEm the parameters are set
      Sigma.corr <- calcSigma(1:nyears, 1:nyears, cov.method = "matern")
      parnames <- c("beta0","sigma.g","mu.y", "G")
      jags.data <- list(y.i = df$y, gett.i = df$t, nyears=nyears, n = length(df$t), Sigma.corr = Sigma.corr)
      if(measurement.err==TRUE){
        jags.data$nu.i <- nu.i
        parnames <- c(parnames, "sigma.y")
      }
    }
    if(matern.cov==FALSE){
      if(is.null(model.file.path)){
        ifelse(measurement.err, model.file.path <- "R/models/model_gp.txt",
               model.file.path <- "R/models/model_gp_nme.txt")
      }
      ## need to calculate distance matrix
      Dist <- rdist(1:nyears)
      ## currently using the reparameterized version
      parnames <- c("beta0","sigma.g","p","mu.y", "G")
      jags.data <- list(y.i = df$y, gett.i = df$t, nyears=nyears, n = length(df$t),kappa=2, Dist = Dist)
      if(measurement.err==TRUE){
        jags.data$nu.i <- nu.i
        parnames <- c(parnames, "sigma.y")
      }
    }
  }

  ## run the model
  cat("Running model.\n")
  mod <- jags.parallel(data = jags.data,
                       parameters.to.save= c(parnames),
                       n.chains = nchains,
                       n.burnin = nburnin,
                       n.iter = niter,
                       n.thin =nthin,
                       model.file = model.file.path)
  if(max(mod$BUGSoutput$summary[, c("Rhat")])>1.1) cat("Something hasn't converged.\n")
  return(mod)
}


