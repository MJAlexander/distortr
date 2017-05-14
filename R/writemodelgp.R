#' Write a JAGS model to fit a Gaussian Process model
#'
#' Write a JAGS model to fit Gaussian Process model, with either exponential or Matern covariance function, with or without time trend
#'
#' @param cov.fun either squared exponential ("sqexp") or matern ("matern")
#' @param cs.smoothing whether smoothing paramter is country specific. If `FALSE`, smoothing parameter is global.
#' @param time.trend if `TRUE` a linear time trend is estimated.
#' @param nserror.estimated whether to estimate non-sampling error. IF `FALSE`, fixed sampling error is inputted.
#' @param file.name name of file to be saved. Must be a `.txt` file
#' @export
#' @return A text file that contains a JAGS model
#' @examples
#' cov.fun <- "sqexp"
#' cs.smoothing <- T
#' time.trend <- T
#' nserror.estimated <- T
#' writeModelGP(cov.fun = cov.fun, cs.smoothing = cs.smoothing, time.trend = time.trend, nserror.estimated = nserror.estimated)

writeModelGP <- function( # Write JAGS model out as a .txt file
  cov.fun = "sqexp",
  cs.smoothing = T,
  time.trend =T,
  nserror.estimated = T,
  file.name = "model.txt"
)
{
  # start txt-file:
  cat("model{  ", file = file.path(file.name), fill = T, append = FALSE)
  cat("
      for(c in 1:niso){
      # data
      for (i in 1:n.c[c]){
      y.i[c,i] ~ dmnorm(G[gett.ci[c,i],c],nu.ci[c,i])
      nu.ci[c,i] <- pow((se.ci[c,i]^2+sigma.y[source.ci[c,i]]^2), -1)
      }
      G[1:nyears.c[c],c] ~ dmnorm(mu.g[1:nyears.c[c],c],Sigma.inv[1:nyears.c[c],1:nyears.c[c],c]) ## gaussian process
      Sigma.inv[1:nyears.c[c],1:nyears.c[c],c] <- inverse(Sigma[1:nyears.c[c],1:nyears.c[c],c]) ##Var-Covar Matrix
      beta[c] ~ dnorm(mu.beta[region.c[c]], tau.beta[region.c[c]])
      ", file = file.path(file.name), fill = T, append = T)
  if(cs.smoothing){
    cat("
        for(t in 1:nyears.c[c]){
          Sigma[t,t,c] <- pow(tau.g[c],-1) + 0.00001  ##avoids issue of non positive definite matrix
         ", file = file.path(file.name), fill = T, append = T)
    if(cov.fun == "matern"){
      cat("
        for(j in (t+1):nyears.c[c]) {
            Sigma[t,j,c]<- pow(tau.g[c],-1)*Sigma.corr[t,j,c]
            Sigma[j,t,c] <- Sigma[t,j,c]
          } #End j loop
          } #End t loop
            tau.g[c]<-pow(sigma.g[c],-2)
            sigma.g[c] <- exp(logsigma.g[c])
            logsigma.g[c] ~ dnorm(chi[region.c[c]], psi[region.c[c]])
          } #end c
          for(r in 1:nregions){
            chi[r] ~ dnorm(chi.global, tau.chi.global)
            psi[r] <- pow(sigma.psi[r],-2)
            sigma.psi[r] ~ dunif(0,40)
          }
          chi.global ~ dnorm(0, 0.01)
          tau.chi.global <- pow(sigma.chi.global, -2)
          sigma.chi.global ~ dunif(0, 40)
          ", file = file.path(file.name), fill = T, append = T)
    }
    if(cov.fun == "sqexp"){
      cat("
        for(j in (t+1):nyears.c[c]) {
          Sigma[t,j,c]<- pow(tau.g[c],-1)*(pow(p[c],pow(Dist[t,j,c],kappa)))
          Sigma[j,t,c] <- Sigma[t,j,c]
        } #End j loop
        } #End t loop
          p[c]~dnorm(mu.p[region.c[c]],tau.p[region.c[c]])T(0,1)
          tau.g[c]<-pow(sigma.g[c],-2)
          sigma.g[c] <- exp(logsigma.g[c])
          logsigma.g[c] ~ dnorm(chi[region.c[c]], psi[region.c[c]])
        } #end c
          for(r in 1:nregions){
          mu.p[r] ~ dnorm(mu.p.global, tau.p.global)
          tau.p[r] <- pow(sigma.p[r], -2)
          sigma.p[r] ~ dunif(0, 40)
          chi[r] ~ dnorm(chi.global, tau.chi.global)
          psi[r] <- pow(sigma.psi[r],-2)
          sigma.psi[r] ~ dunif(0,40)
          }
          mu.p.global ~ dnorm(0, 0.01)
          tau.p.global <- pow(sigma.p.global, -2)
          sigma.p.global ~ dunif(0, 40)

          chi.global ~ dnorm(0, 0.01)
          tau.chi.global <- pow(sigma.chi.global, -2)
          sigma.chi.global ~ dunif(0, 40)
          ", file = file.path(file.name), fill = T, append = T)
    }
  }
  if(!cs.smoothing){
    cat("
        for(t in 1:nyears.c[c]){
        Sigma[t,t,c] <- pow(tau.g,-1) + 0.00001  ##avoids issue of non positive definite matrix
      ", file = file.path(file.name), fill = T, append = T)
    if(cov.fun=="sqexp"){
      cat("
        for(j in (t+1):nyears.c[c]) {
          Sigma[t,j,c]<- pow(tau.g,-1)*(pow(p,pow(Dist[t,j,c],kappa)))
          Sigma[j,t,c] <- Sigma[t,j,c]
    } #End j loop
          ", file = file.path(file.name), fill = T, append = T)
    }
    if(cov.fun=="matern"){
      cat("
        for(j in (t+1):nyears.c[c]) {
          Sigma[t,j,c]<- pow(tau.g,-1)*Sigma.corr[t,j,c]
          Sigma[j,t,c] <- Sigma[t,j,c]
    } #End j loop
          ", file = file.path(file.name), fill = T, append = T)
    }
    cat("
    } #End t loop
        } #end c
        tau.g <- pow(sigma.g, -2)
        sigma.g ~ dunif(0,40)
        p ~ dunif(0, 1)
        ", file = file.path(file.name), fill = T, append = T)
  }
  cat("
      for(r in 1:nregions){
      mu.beta[r] ~ dnorm(mu.beta.global, tau.beta.global)
      tau.beta[r] <- pow(sigma.beta[r], -2)
      sigma.beta[r] ~ dunif(0, 40)
      }
      mu.beta.global ~ dnorm(0, 0.01)
      tau.beta.global <- pow(sigma.beta.global, -2)
      sigma.beta.global ~ dunif(0, 40)
      ", file = file.path(file.name), fill = T, append = T)
  if(time.trend){
    cat("
        for(c in 1:niso){
        gamma[c] ~ dnorm(mu.gamma[region.c[c]],tau.gamma[region.c[c]])
        for(t in 1:nyears.c[c]){
          mu.g[t,c] <- beta[c] + gamma[c]*t
        }
        }
        for(r in 1:nregions){
        # hierarchy for gammas
        mu.gamma[r] ~ dnorm(mu.gamma.global, tau.gamma.global)
        tau.gamma[r] <- pow(sigma.gamma[r], -2)
        sigma.gamma[r] ~ dunif(0, 40)
        } # end r
        mu.gamma.global ~ dnorm(0, 0.01)
        tau.gamma.global <- pow(sigma.gamma.global, -2)
        sigma.gamma.global ~ dunif(0, 40)
        ", file = file.path(file.name), fill = T, append = T)
  }
  if(!time.trend){
    cat("
        for(c in 1:niso){
        for(t in 1:nyears.c[c]){
        mu.g[t,c] <- beta[c]
        }
  }
        ", file = file.path(file.name), fill = T, append = T)
  }
  if(nserror.estimated){
    cat("
        for(s in 1:nsources){
        sigma.y[s] ~ dunif(0,40)
        }
        ", file = file.path(file.name), fill = T, append = T)
  }
  # close model file
  cat("} # end model ", file = file.path(file.name), fill = T, append = T)
  return(invisible())
  }

