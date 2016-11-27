#' Plot results from MCMC estimation
#'
#' @param res.df A data frame of x values, estimates and uncertainty intervals
#' @param method The method of smoothing to implement (choices: ar, arma, splines, gp)
#' @param order The order of splines penalization (either 1 or 2)
#' @param save.plot Whether or not to save plot. Default is \code{FALSE}
#' @param save.file.path Directory to save file
#' @export
#' @return A plot of time series data and estimates.
#' @seealso \code{\link{runMCMC, getResults}}
#' @examples
#' df.mu <- getResults(mod, method = "splines")
#' plotResults(df.mu, method = "splines", order = 1)

plotResults <- function(res.df,
                        method,
                        order = NULL,
                        save.plot = F,
                        save.file.path = "output"){
  if(method=="splines" & is.null(order)) stop("Order of penalization must be specified.")
  # plot results
  p <- ggplot(data = res, aes(x = t, y = y)) + geom_point() + theme_bw()+
    geom_line(data = df.mu, aes(x = time, y = median), color = "red")+
    geom_ribbon(data=df.mu,aes(x=time,y=NULL,ymin=lower, ymax=upper), alpha = 0.2, fill = "red")
  print(p)
  if(save.plot==T){
    cat("Saving plot.\n ")
    dir.create(file.path(save.file.path, "/plots/"), showWarnings = FALSE)
    ggsave(filename = paste0(save.file.path, "/plots/", method, order,".pdf"), plot = p)
  }
}
