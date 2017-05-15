#' Calculate WAIC from log-likelihood
#'
#' Calculate Watanabe-Akaike or widely available information criterion from log posterior predictive distibution.
#'
#' @param log_lik log posterior predictive distibution. Must be of dimensions \code{i x s} where i is observations
#' and s is number of posterior draws.
#' @export
#' @return A list that contains values for WAIC, p_waic
#' @note If observations are for many countries, do a \code{cbind} to get a matrix of dimensions \code{i x s}
#' @author Leontine Alkema

waic <- function(log_lik){ # log-like should be the log(p(y_i|theta^(s)) with i=1,..,n rows and s=1,2,..,S columns
  log_lik
  dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw/
    matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
  loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
                    colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n*colVars(pointwise))
  return(list(waic=total["waic"], elpd_waic=total["elpd_waic"],
              p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
              pointwise=pointwise, total=total, se=se))
}

#' Calculate variance of columns
#'
#' Calculate variance of columns of a matrix.
#'
#' @param a a numeric matrix
#' @export
#' @return A numeric value
#' @author Leontine Alkema


colVars <- function(a){
  n <-   dim(a)[[1]];
  c <- dim(a)[[2]];
  return(.colMeans(((a - matrix(.colMeans(a, n, c), nrow = n,
                                ncol = c, byrow = TRUE)) ^ 2), n, c) * n / (n - 1))
}

