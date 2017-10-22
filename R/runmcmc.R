#' Run MCMC estimation
#'
#' Run MCMC estimation for time series using JAGS. Can either be run on a single time series of a set of time series.
#'
#' @param input.data Input data to JAGS.
#' If single country, this is a dataframe of x and y observations, and standard errors around ys.
#' If a global run, this is a list of required input data. See \code{processData}, \code{getSplinesData} and \code{getGPData} to get required data in compatible form.
#' @param method The method of smoothing to implement (choices: ar, arma, splines, gp)
#' @param nyears For single country runs: number of years of observations
#' @param obs.err For single country runs: is TRUE if standard errors are observed
#' @param measurement.err For single country runs: is TRUE if there is assumed to be measurement error
#' @param cs.arma For global runs: whether ARMA parameter(s) are country specific. If `FALSE`, parameter is global.
#' @param cs.smoothing For global runs: whether smoothing paramter is country specific. If `FALSE`, smoothing parameter is global.
#' @param time.trend For global runs: if `TRUE` a linear time trend is estimated.
#' @param nserror.estimated For global runs: whether to estimate non-sampling error. IF `FALSE`, fixed sampling error is inputted.
#' @param order The order of splines penalization (either 1 or 2)
#' @param I Knot spacing for splines
#' @param matern.cov Whether or not to use Matern covariance function if \code{method=="gp"}. Default is \code{TRUE}.
#' @param nchains Number of MCMC chains
#' @param nburnin Number of iterations to throw away as burn in.
#' @param niter Number of total iterations.
#' @param nthin Degree of thinning of MCMC chains
#' @param model.file.path Text file which contains the model to be fitted. If \code{NULL}, the text file is drawn from the \code{models} folder.
#' @param model.save.file.path For global runs: path to save model, if written.
#' @export
#' @return A JAGS model object
#' @seealso \code{\link{getResults}, \link{plotResults}}
#' @examples
#' nyears <- 100
#' prop.sample <- 0.7
#' obs.err <- TRUE
#' sigma.y <- 0.5
#' seed <- 123
#' method <- 'splines'
#' params <- list(sigma.alpha = 1, order = 1)
#' df <- simulateFluctuations(nyears, prop.sample, method, params, obs.err, sigma.y)
#' df$se <- 1
#' mod <- runMCMC(input.data = df, nyears = 100, method = "splines", order = 1,nchains = 4, nburnin = 100, niter = 100+3000, nthin = 3)


runMCMC <- function(input.data,
                    method,
                    nyears = NULL,
                    obs.err = TRUE,
                    measurement.err = TRUE,
                    cs.arma = NULL,
                    cs.smoothing = TRUE,
                    time.trend = FALSE,
                    nserror.estimated = TRUE,
                    order = NULL,
                    I = 2.5,
                    matern.cov=TRUE,
                    nchains = 3,
                    nburnin = 1000,
                    niter = 2000,
                    nthin = 1,
                    model.file.path = NULL,
                    model.save.file.path = "R/model.txt"){

  if(is.null(dim(input.data))){
    mod <- runMCMCGlobal(method = method,
                         input.data = input.data,
                         order = order,
                         matern.cov=matern.cov,
                         cs.arma = cs.arma,
                         cs.smoothing = cs.smoothing,
                         time.trend = time.trend,
                         nserror.estimated = nserror.estimated,
                         nchains = nchains,
                         nburnin = nburnin,
                         niter = niter,
                         nthin = nthin,
                         model.file.path = model.file.path,
                         model.save.file.path = model.save.file.path)
  }
  if(length(dim(input.data))==2){
    mod <- runMCMCCountry(df = input.data,
                          nyears = nyears,
                          method = method,
                          order = order,
                          I = I,
                          matern.cov=matern.cov,
                          obs.err = obs.err,
                          measurement.err = measurement.err,
                          nchains = nchains,
                          nburnin = nburnin,
                          niter = niter,
                          nthin = nthin,
                          model.file.path = model.file.path)
  }
  return(mod)
}
