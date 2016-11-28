#' Plot results from MCMC estimation
#'
#' @param data.df A data frame of x values and data values
#' @param res.df A data frame of x values, estimates and uncertainty intervals
#' @param method The method of smoothing to implement (choices: ar, arma, splines, gp)
#' @param order The order of splines penalization (either 1 or 2)
#' @param save.plot Whether or not to save plot. Default is \code{FALSE}
#' @param save.file.path Directory to save file
#' @export
#' @return A plot of time series data and estimates.
#' @seealso \code{\link{runMCMC, getResults}}
#' @examples
#' #' nyears <- 100
#' prop.sample <- 0.7
#' obs.err <- T
#' sigma.y <- 0.5
#' seed <- 123
#' method <- 'splines'
#' params <- list(sigma.alpha = 1, order = 1)
#' res <- simulateFluctuations(nyears, prop.sample, method, params, obs.err, sigma.y)
#' mod <- runMCMC(df = res, nyears = 100, method = "splines", order = 1, nchains = 4, nburnin = 100, niter = 100+3000, nthin = 3)
#' df.mu <- getResults(mod, method = "splines")
#' plotResults(res, df.mu, method = "splines", order = 1)

plotResults <- function(data.df, res.df,
                        method,
                        order = NULL,
                        save.plot = F,
                        save.file.path = "output"){
  if(method=="splines" & is.null(order)) stop("Order of penalization must be specified.")
  # plot results
  p <- ggplot(data = data.df, aes(x = t, y = y)) + geom_point() + theme_bw()+
    geom_line(data = df.mu, aes(x = time, y = median), color = "red")+
    geom_ribbon(data=df.mu,aes(x=time,y=NULL,ymin=lower, ymax=upper), alpha = 0.2, fill = "red")
  print(p)
  if(save.plot==T){
    cat("Saving plot.\n ")
    dir.create(file.path(save.file.path, "/plots/"), showWarnings = FALSE)
    ggsave(filename = paste0(save.file.path, "/plots/", method, order,".pdf"), plot = p)
  }
}
