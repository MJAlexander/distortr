#' Plot time series of data
#'
#' @param data.df A data frame of t values and data values, with standard errors
#' @param maintitle Title of plot
#' @param plot.se Whether or not to plot standard errors
#' @export
#' @return A plot of time series data with standard errors, if desired.
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
#' plotData(res)

plotData <- function(data.df,
                        maintitle = NULL,
                        plot.se = T,
                        ...){
  if(!("t" %in% colnames(data.df))){
    # assume column order is t, y, se, source
    colnames(data.df) <- c("t", "y", "se", "source")
  }

  # plot results
  p <- ggplot(data = data.df, aes(x = t, y = y, color = source))
  if(plot.se==T){
    p <- p+ geom_errorbar(data=data.df,aes(x=t,y=NULL,ymin=y-2*se, ymax=y+2*se), width=0.2, color = "grey")
  }
  p <- p + geom_point()+
    ggtitle(maintitle)+
    theme_bw()
  if(length(levels(data.df$source))>1){
    myColors <- brewer.pal(length(levels(data.df$source)),"Set2")
    names(myColors) <- levels(data.df$source)
    colScale <- scale_colour_manual(name = "source", values = myColors)
    p <- p+  colScale + guides(color=guide_legend("Data source"))
  }

  return(p)
}
