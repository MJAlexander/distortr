#' Run validation of MCMC model
#'
#' Calculate some model validation measures for MCMC models. These are calculated by leaving out some of the available sample. The function returns values for root-meean-squared-error, coverage and interval score.
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
#' @param leave.out.method If equal to \code{random}, oberservations are left out at random. If equal to \code{recent} the most recent x\% of observations are left out.
#' @param leave.out.percent The percent of observations to leave out. Default is 20\%
#' @param nreps The number of times to repeat the model validation if the leave out method is random
#' @param alpha.level Significance level of uncertainty intervals. Default is 5\%.
#' @param alpha.n Significant level of interval score. Default is 10\%.
#' @param seed Value of random seed.
#' @export
#' @return A list containing values of RMSE, a vector of coverage values for each repetition, and a vector of interval scores for each repetition.
#' @seealso \code{\link{runMCMC}, \link{simulateFluctuations}}
#' @examples
#' nyears <- 100
#' prop.sample <- 0.7
#' obs.err <- T
#' sigma.y <- 0.5
#' seed <- 123
#' method <- 'splines'
#' params <- list(sigma.alpha = 1, order = 1)
#' df <- simulateFluctuations(nyears, prop.sample, method, params, obs.err, sigma.y)
#' validation.results <- runModelValidation(df = res, nyears = 100, method = "splines", order = 1,
#' nchains = 4, nburnin = 100, niter = 100+3000, nthin = 3,
#' leave.out.method = "recent", leave.out.percent = 20)

runModelValidation <- function(df,
                               nyears,
                               method,
                               order = NULL,
                               nchains = 4,
                               nburnin = 1000,
                               niter = 1000+30000,
                               nthin =30,
                               model.file.path = NULL,
                               leave.out.method,
                               leave.out.percent = 20,
                               nreps = 10,
                               alpha.level = 0.05,
                               alpha.n = 0.1,
                               seed = 123){
  if(leave.out.method=="random"){
    list.sampled <- list()
    list.lo <- list()
    set.seed(seed)
    # need to do this many times and take an average of the validation measures
    for(i in 1:nreps){
      rows.to.sample <- sort(sample(1:nrow(df), round(nrow(df)*(100-leave.out.percent)/100)))
      rows.to.lo <- (1:nrow(df) %in% rows.to.sample)==FALSE
      list.sampled[[i]] <- df[rows.to.sample,]
      list.lo[[i]] <- df[rows.to.lo,]
    }
  }

  if(leave.out.method=="recent"){
    nreps <- 1
    list.sampled <- list()
    list.lo <- list()
    nobs.to.keep <- round(nrow(df)*(100-leave.out.percent)/100)
    list.sampled[[1]] <- df[1:nobs.to.keep,]
    list.lo[[1]] <- df[((nobs.to.keep+1):nrow(df)),]
  }

  mod.full <- runMCMC(df = df, nyears = nyears,
                      method = method, order = order,
                      nchains = nchains, nburnin = nburnin,
                      niter = niter, nthin = nthin)
  df.mu.full <- getResults(mod.full, method = method, alpha.level = alpha.level)

  rmse.s <- c()
  coverage.s <- c()
  intscores.s <- c()
  for(i in 1:nreps){
    mod.lo <- runMCMC(df = list.sampled[[i]], nyears = nyears,
                      method = method, order = order,
                      nchains = nchains, nburnin = nburnin,
                      niter = niter, nthin = nthin)
    df.mu.lo <- getResults(mod.lo, method = method, alpha.level = alpha.level)

    # root mean squared error
    rmse <- sum(sqrt((df.mu.full$median - df.mu.lo$median)^2))
    rmse.s <- c(rmse.s, rmse)

    # coverage and sharpness
    y.lo <- list.lo[[i]]$y
    t.lo <- list.lo[[i]]$t
    in.ui <- c()
    int.scores <- c()
    for(j in list.lo[[i]]$t){
      # coverage
      in.ui.j <- y.lo[t.lo==j] > df.mu.lo$lower[df.mu.lo$time==j] & y.lo[t.lo==j] < df.mu.lo$upper[df.mu.lo$time==j]
      in.ui <- c(in.ui, in.ui.j)
      # sharpness
      n.j <- (df.mu.lo$upper[df.mu.lo$time==j] -
                df.mu.lo$lower[df.mu.lo$time==j]) + 2/alpha.n*(df.mu.lo$lower[df.mu.lo$time==j] -
                                                                 y.lo[t.lo==j])*(y.lo[t.lo==j] < df.mu.lo$lower[df.mu.lo$time==j]) + 2/alpha.n*(y.lo[t.lo==j] -
                                                                                                                                                  df.mu.lo$upper[df.mu.lo$time==j])*(y.lo[t.lo==j] > df.mu.lo$upper[df.mu.lo$time==j])
      int.scores <- c(int.scores, n.j)
    }
    coverage.s <- c(coverage.s, mean(in.ui))
    intscores.s <- c(intscores.s, mean(int.scores))
  }
  return(list(rmse= rmse.s, coverage = coverage.s, int.score = intscores.s))
}


