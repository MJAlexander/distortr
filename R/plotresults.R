#' Plot results from MCMC estimation
#'
#' @param data.df A data frame of t values and data values, with standard errors
#' @param res.df A data frame of t values, estimates and uncertainty intervals
#' @param res2.df A data frame of t values, estimates and uncertainty intervals. Default is NULL
#' @param maintitle Title of plot
#' @param plot.se Whether or not to plot standard errors
#' @param save.plot Whether or not to save plot. Default is \code{FALSE}
#' @param save.file.name Directory/file name where to save file
#' @export
#' @return A plot of time series data and estimates.
#' @seealso \code{\link{runMCMC}, \link{getResults}}
#' @examples
#' nyears <- 100
#' prop.sample <- 0.7
#' obs.err <- TRUE
#' sigma.y <- 0.5
#' seed <- 123
#' method <- 'splines'
#' params <- list(sigma.alpha = 1, order = 1)
#' res <- simulateFluctuations(nyears, prop.sample, method, params, obs.err, sigma.y)
#' res$se <- 0.1
#' mod <- runMCMC(input.data = res, nyears = 100, method = "splines", order = 2, nchains = 2, nburnin = 100, niter = 100+3000, nthin = 3)
#' df.mu <- getResults(mod, method = "splines", nyears = nyears)
#' plotResults(res, df.mu, method = "splines", order = 1, save.file.name = "test.pdf")

plotResults <- function(data.df, res.df, res2.df = NULL,
                        method,
                        order = NULL,
                        maintitle = NULL,
                        plot.se = T,
                        save.plot = T,
                        save.file.name = "output/plots/myplot.pdf", ...){

  data.df <- as.data.frame(data.df)
  # plot results
  p <- plotData(data.df = data.df, maintitle = maintitle, plot.se = plot.se)
  p <- p +
    geom_line(data = res.df, aes(x = time, y = median), color = "red")+
    geom_ribbon(data = res.df, aes(x=time, y=NULL, ymin=lower, ymax=upper),
                alpha = 0.2, fill = "red", inherit.aes=FALSE)
  if(!is.null(res2.df)){
    p <- p + geom_line(data = res2.df, aes(x = time, y = median), color = "blue")+
      geom_ribbon(data=res2.df,aes(x=time,y=NULL,ymin=lower, ymax=upper),
                  alpha = 0.2, fill = "blue", inherit.aes=FALSE)
  }
  if(save.plot==T){
    cat("Saving plot.\n ")
    dir.create(file.path("output/plots/"), showWarnings = FALSE)
    ggsave(filename = save.file.name, plot = p)
  }
  return(p)
}
