#' Run MCMC estimation
#'
#' Run MCMC estimation of time series using JAGS.
#'
#' @param df A dataframe of x and y observations
#' @param nyears number of years of observations
#' @param method The method of smoothing to implement (choices: ar, arma, splines, gp)
#' @param order The order of splines penalization (either 1 or 2)
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
#' res <- simulateFluctuations(nyears, prop.sample, method, params, obs.err, sigma.y)
#' mod <- runMCMC(df = res, nyears = 100, method = "splines", order = 1,
#' nchains = 4, nburnin = 100, niter = 100+3000, nthin = 3)

runMCMC <- function(df,
                    nyears,
                    method,
                    order = NULL,
                    nchains = 4,
                    nburnin = 1000,
                    niter = 1000+30000,
                    nthin =30,
                    model.file.path = NULL){

  if(method=="ar"){
    if(is.null(model.file.path)){
      model.file.path <- "R/models/model_ar.txt"
    }
    jags.data <- list(y.i = res$y, gett.i = res$t, nyears=nyears, n = length(res$t))
    parnames <- c("sigma", "rho", "sigma.y", "mu.t")
  }

  if(method=="arma"){
    if(is.null(model.file.path)){
      model.file.path <- "R/models/model_arma.txt"
    }
    jags.data <- list(y.i = res$y, gett.i = res$t, nyears=nyears, n = length(res$t))
    parnames <- c("sigma.ar", "phi", "theta", "sigma.y", "mu.t")
  }

  if(method=="splines"){
    if(is.null(order)){
      stop("Order of penalization must be specified.")
    }
    if(order==1){
      if(is.null(model.file.path)){
        model.file.path <- "R/models/model_splines_1.txt"
      }
    }
    if(order==2){
      if(is.null(model.file.path)){
        model.file.path <- "R/models/model_splines_2.txt"
      }
    }
    x.t <- 1:nyears
    sp <- GetSplines(x.t)
    K <- length(sp$knots.k)
    B.tk <- sp$B.ik
    jags.data <- list(y.i = res$y, gett.i = res$t, nyears=nyears, n = length(res$t), K = K, B.tk = B.tk)
    parnames <- c("alpha.k", "sigma.alpha", "sigma.y", "mu.t")
  }

  if(method=="gp"){
    if(is.null(model.file.path)){
      model.file.path <- "R/models/model_gp_4.txt"
    }
    ## need to calculate distance matrix
    Dist <- rdist(1:nyears)
    ## currently using the reparameterized version
    parnames <- c("beta0","sigma.y","sigma.g","p","mu.y", "G")
    jags.data <- list(y.i = res$y, gett.i = res$t, nyears=nyears, n = length(res$t),kappa=2, Dist = Dist)
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


