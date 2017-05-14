#' Write a JAGS model to fit AR(1) model
#'
#' Write a JAGS model to fit AR(1) model, with or without time trend
#'
#' @param cs.arma whether AR parameter is country specific. If `FALSE`, parameter is global.
#' @param cs.smoothing whether smoothing paramter is country specific. If `FALSE`, smoothing parameter is global.
#' @param time.trend if `TRUE` a linear time trend is estimated.
#' @param nserror.estimated whether to estimate non-sampling error. IF `FALSE`, fixed sampling error is inputted.
#' @param file.name name of file to be saved. Must be a `.txt` file
#' @export
#' @return A text file that contains a JAGS model
#' @examples
#' cs.arma <- T
#' cs.smoothing <- T
#' time.trend <- T
#' nserror.estimated <- T
#' writeModelAR(cs.arma = cs.arma, cs.smoothing = cs.smoothing, time.trend = time.trend, nserror.estimated = nserror.estimated)

writeModelAR <- function( # Write JAGS model out as a .txt file
  cs.arma = T,
  cs.smoothing = T,
  time.trend =T,
  nserror.estimated = T,
  file.name = "model.txt"
)
{
  # start txt-file:
  cat("model{  ", file = file.path(file.name), fill = T, append = FALSE)
  if(time.trend){
    cat("
      for(c in 1:niso){
        # data
        for (i in 1:n.c[c]){
        y.ci[c,i] ~ dnorm(yhat.ci[c,i], nu.ci[c,i])
        yhat.ci[c,i] <- mu.ct[c,gett.ci[c,i]] + gett.ci[c,i]*gamma[c]
        nu.ci[c,i] <- pow((se.ci[c,i]^2+sigma.y[source.ci[c,i]]^2), -1)
        }
        ", file = file.path(file.name), fill = T, append = T)
  }
  if(!time.trend){
    cat("
        for(c in 1:niso){
        # data
        for (i in 1:n.c[c]){
        y.ci[c,i] ~ dnorm(yhat.ci[c,i], nu.ci[c,i])
        yhat.ci[c,i] <- mu.ct[c,gett.ci[c,i]]
        nu.ci[c,i] <- pow((se.ci[c,i]^2+sigma.y[source.ci[c,i]]^2), -1)
        }
        ", file = file.path(file.name), fill = T, append = T)
  }
  if(cs.arma&cs.smoothing){
    cat("
      mu.ct[c,1] ~ dnorm(beta[c], tau.stat[c])
      tau.stat[c] <- (1-pow(rho[c],2))/pow(sigma[c],2)
      for (t in 2:nyears.c[c]){
        mu.ct[c,t] ~ dnorm(muhat.ct[c,t], tau[c])
        muhat.ct[c,t] <- rho[c]*mu.ct[c,t-1]
      } # end t
      rho[c] ~ dnorm(mu.rho[region.c[c]], tau.rho[region.c[c]])T(-1,1)
      tau[c] <- pow(sigma[c],-2)
      sigma[c] <- exp(logsigma[c])
      logsigma[c] ~ dnorm(chi.sigma[region.c[c]], tau.sigma[region.c[c]])
      beta[c] ~ dnorm(mu.beta[region.c[c]],tau.beta[region.c[c]])
    } # end c
    for(r in 1:nregions){
      # hierarchy for betas
      mu.beta[r] ~ dnorm(mu.beta.global, tau.beta.global)
      tau.beta[r] <- pow(sigma.beta[r], -2)
      sigma.beta[r] ~ dunif(0, 40)

      # hierarchy for rhos
      mu.rho[r] ~ dnorm(mu.rho.global, tau.rho.global)
      tau.rho[r] <- pow(sigma.rho[r], -2)
      sigma.rho[r] ~ dunif(0, 40)

      # hierarchy for smoother
      chi.sigma[r] ~ dnorm(mu.chi.global, tau.chi.global)
      tau.sigma[r] <- pow(sigma.sigma[r], -2)
      sigma.sigma[r] ~ dunif(0, 40)

    } # end r
    mu.beta.global ~ dnorm(0, 0.01)
    tau.beta.global <- pow(sigma.beta.global, -2)
    sigma.beta.global ~ dunif(0, 40)


    mu.rho.global ~ dnorm(0, 0.01)
    tau.rho.global <- pow(sigma.rho.global, -2)
    sigma.rho.global ~ dunif(0, 40)


    mu.chi.global ~ dnorm(0, 0.01)
    tau.chi.global <- pow(sigma.chi.global, -2)
    sigma.chi.global ~ dunif(0, 40)

        ", file = file.path(file.name), fill = T, append = T)
  }
  if(cs.arma&!cs.smoothing){
    cat("
      mu.ct[c,1] ~ dnorm(beta[c], tau.stat[c])
      tau.stat[c] <- (1-pow(rho[c],2))/pow(sigma,2)
      for (t in 2:nyears.c[c]){
        mu.ct[c,t] ~ dnorm(muhat.ct[c,t], tau)
        muhat.ct[c,t] <- rho[c]*mu.ct[c,t-1]
      } # end t
      rho[c] ~ dnorm(mu.rho[region.c[c]], tau.rho[region.c[c]])T(-1,1)
      beta[c] ~ dnorm(mu.beta[region.c[c]],tau.beta[region.c[c]])
    } # end c
    for(r in 1:nregions){
      # hierarchy for betas
        mu.beta[r] ~ dnorm(mu.beta.global, tau.beta.global)
        tau.beta[r] <- pow(sigma.beta[r], -2)
        sigma.beta[r] ~ dunif(0, 40)

        # hierarchy for rhos
        mu.rho[r] ~ dnorm(mu.rho.global, tau.rho.global)
        tau.rho[r] <- pow(sigma.rho[r], -2)
        sigma.rho[r] ~ dunif(0, 40)

    } # end r
    mu.beta.global ~ dnorm(0, 0.01)
    tau.beta.global <- pow(sigma.beta.global, -2)
    sigma.beta.global ~ dunif(0, 40)

    mu.rho.global ~ dnorm(0, 0.01)
    tau.rho.global <- pow(sigma.rho.global, -2)
    sigma.rho.global ~ dunif(0, 40)

    tau <- pow(sigma, -2)
    sigma ~ dunif(0, 40)
        ", file = file.path(file.name), fill = T, append = T)
  }
  if(!cs.arma&cs.smoothing){
    cat("
        mu.ct[c,1] ~ dnorm(beta[c], tau.stat[c])
        tau.stat[c] <- (1-pow(rho,2))/pow(sigma[c],2)
        for (t in 2:nyears.c[c]){
        mu.ct[c,t] ~ dnorm(muhat.ct[c,t], tau[c])
        muhat.ct[c,t] <- rho*mu.ct[c,t-1]
        } # end t
        tau[c] <- pow(sigma[c],-2)
        sigma[c] <- exp(logsigma[c])
        logsigma[c] ~ dnorm(chi.sigma[region.c[c]], tau.sigma[region.c[c]])
        beta[c] ~ dnorm(mu.beta[region.c[c]],tau.beta[region.c[c]])
  } # end c
        for(r in 1:nregions){
        # hierarchy for betas
        mu.beta[r] ~ dnorm(mu.beta.global, tau.beta.global)
        tau.beta[r] <- pow(sigma.beta[r], -2)
        sigma.beta[r] ~ dunif(0, 40)

        # hierarchy for smoother
        chi.sigma[r] ~ dnorm(mu.chi.global, tau.chi.global)
        tau.sigma[r] <- pow(sigma.sigma[r], -2)
        sigma.sigma[r] ~ dunif(0, 40)

        } # end r
        mu.beta.global ~ dnorm(0, 0.01)
        tau.beta.global <- pow(sigma.beta.global, -2)
        sigma.beta.global ~ dunif(0, 40)

        rho ~ dunif(-1,1)

        mu.chi.global ~ dnorm(0, 0.01)
        tau.chi.global <- pow(sigma.chi.global, -2)
        sigma.chi.global ~ dunif(0, 40)
        ", file = file.path(file.name), fill = T, append = T)
  }
  if(!cs.arma&!cs.smoothing){
    cat("
      mu.ct[c,1] ~ dnorm(beta[c], tau.stat)
        for (t in 2:nyears.c[c]){
        mu.ct[c,t] ~ dnorm(muhat.ct[c,t], tau)
        muhat.ct[c,t] <- rho*mu.ct[c,t-1]
        } # end t
        beta[c] ~ dnorm(mu.beta[region.c[c]],tau.beta[region.c[c]])
      } # end c
      for(r in 1:nregions){
          # hierarchy for betas
            mu.beta[r] ~ dnorm(mu.beta.global[], tau.beta.global)
            tau.beta[r] <- pow(sigma.beta[r], -2)
            sigma.beta[r] ~ dunif(0, 40)

      } # end r
      mu.beta.global ~ dnorm(0, 0.01)
      tau.beta.global <- pow(sigma.beta.global, -2)
      sigma.beta.global ~ dunif(0, 40)

      tau.stat <- (1-pow(rho,2))/pow(sigma,2)
      rho ~ dunif(-1,1)
      tau <- pow(sigma, -2)
      sigma ~ dunif(0, 40)
        ", file = file.path(file.name), fill = T, append = T)
  }
  if(time.trend){
    cat("
        for(c in 1:niso){
          gamma[c] ~ dnorm(mu.gamma[region.c[c]],tau.gamma[region.c[c]])
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

