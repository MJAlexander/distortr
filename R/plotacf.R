#' Plot sample autocorrelation function of time series
#'
#' Plot sample ACF. The function is different to standard R function as it allows for missing values
#' @param df A data frame of x values and data values, with standard errors
#' @param nyears number of years of observation period
#' @export
#' @return A plot of the sample autocorrelation function, with 95\% confidence intervals.
#' @references Cryer, J. D. 1986. Time series analysis. Boston: Duxbury Press.
#' @examples
#'nyears <- 50
#'prop.sample <- .8
#'obs.err <- F
#'sigma.y <- 1
#'seed <- 1234
#'method <- "ar"
#'params <- list(rho = 0.95, sigma = 0.5, ystart = NULL, eps0.t = NULL)
#'df <- simulateFluctuations(nyears, prop.sample, method, params, obs.err, sigma.y)
#'plotACF(df, nyears)

plotACF <- function(df, nyears){
  t.na <- sapply(1:nyears, function(i) (i %in% df$t)*i)
  t.na <- ifelse(t.na==0, NA, t.na)

  y.na <- rep(NA, nyears)
  for(i in 1:length(t.na)){
    if(is.na(t.na[i])){
      y.na[i] <- NA
    }
    else{
      y.na[i] <- df$y[df$t==t.na[i]]
    }
  }

  y.bar <- mean(y.na, na.rm=T)
  y.demeaned <- (y.na - y.bar)
  b.jk <- matrix(NA, nrow = length(y.na), ncol = length(y.na))
  for (k in 0:(length(y.na)-1)){
    for (j in 1:(length(y.na)-k)){
      b.jk[j,k+1] <- y.demeaned[j]*y.demeaned[j+k]
    }
  }

  denom <- sum(b.jk[,1], na.rm=T)

  r.k <- sapply(1:length(y.na), function(i) sum(b.jk[,i], na.rm=T)/denom)
  r.k[r.k==0] <- NA


  m0 <- sum(!is.na(y.na))
  se.rk <- sapply(1:length(y.na), function(i) sqrt(sum(!is.na(b.jk[,i]))/((m0+2)*m0)))
  se.rk[se.rk==0] <- NA

  p <- ggplot(data = NULL, aes(x = 1:nyears, y = r.k))  +
    geom_segment(aes(xend=1:nyears), yend=0) +
    geom_path(aes(x = 1:nyears, y = 2*se.rk), col = "red", linetype = "dashed") +
    geom_path(aes(x = 1:nyears, y = -2*se.rk), col = "red", linetype = "dashed")+
    xlab("t") + ylab("Sample ACF") + theme_bw()
  return(p)
}
