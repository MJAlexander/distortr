#' Plot results from MCMC estimation
#'
#' @param data.df A data frame of x values and data values, with standard errors
#' @param res.df A data frame of x values, estimates and uncertainty intervals
#' @param method The method of smoothing to implement (choices: ar, arma, splines, gp)
#' @param order The order of splines penalization (either 1 or 2)
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
#' mod <- runMCMC(df = res, nyears = 100, method = "splines", order = 1,
#' nchains = 4, nburnin = 100, niter = 100+3000, nthin = 3)
#' df.mu <- getResults(mod, method = "splines")
#' plotResults(res, df.mu, method = "splines", order = 1)

plotResults <- function(data.df, res.df,
                        method,
                        order = NULL,
                        maintitle = NULL,
                        plot.se = T,
                        save.plot = T,
                        save.file.name = paste0("output/plots/", method, order,".pdf")){
  if(method=="splines" & is.null(order)) stop("Order of penalization must be specified.")
  # plot results
  p <- ggplot(data = data.df, aes(x = t, y = y))
  if(plot.se==T){
    p <- p+ geom_errorbar(data=data.df,aes(x=t,y=NULL,ymin=y-2*se, ymax=y+2*se), width=0.2, color = "grey")
  }
  p <- p + geom_point() + theme_bw()+
    geom_line(data = res.df, aes(x = time, y = median), color = "red")+
    geom_ribbon(data=res.df,aes(x=time,y=NULL,ymin=lower, ymax=upper), alpha = 0.2, fill = "red")+
    ggtitle(maintitle)
  print(p)
  if(save.plot==T){
    cat("Saving plot.\n ")
    dir.create(file.path("output/plots/"), showWarnings = FALSE)
    ggsave(filename = save.file.name, plot = p)
  }
}
