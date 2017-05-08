#' Write a JAGS model to fit ARMA(1,1) model
#'
#' Write a JAGS model to fit ARMA(1,1) model, with or without time trend
#'
#' @param cs.arma whether autocorrelation and innovations parameter is country specific. If `FALSE`, parameters are global.
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
#' writeModelARMA(cs.arma = cs.arma, cs.smoothing = cs.smoothing, time.trend = time.trend, nserror.estimated = nserror.estimated)

writeModelARMA <- function( # Write JAGS model out as a .txt file
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
        for(c in 1:n.iso){
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
        for(c in 1:n.iso){
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
        mu.ct[c,1] ~ dnorm(beta[c], 1/eta[c])
        e.ct[c,1] ~ dnorm(sigma.ar[c]^2/eta[c]*mu.ct[c,1], 1/(sigma.ar[c]^2*(1-sigma.ar[c]^2/eta[c])))
        for (t in 2:nyears.c[c]){
        e.ct[c,t] ~ dnorm(0, tau.ar[c])
        mu.ct[c,t] <- rho[c]*mu.ct[c,t-1] - theta[c]*e.ct[c,t-1] + e.ct[c,t]
        } # end t
        rho[c] ~ dnorm(mu.rho[region.c[c]], tau.rho[region.c[c]])T(0,1)
        theta[c] ~ dnorm(mu.theta[region.c[c]], tau.theta[region.c[c]])T(-1,0)
        eta[c] <- pow(sqrteta[c],2)
        sqrteta[c] <- exp(logsqrteta[c])
        logsqrteta[c] ~ dnorm(chi.eta[region.c[c]], tau.eta[region.c[c]])
        beta[c] ~ dnorm(mu.beta[region.c[c]],tau.beta[region.c[c]])
        sigma.ar[c] <- sqrt(eta[c] /((1-2*rho[c]*theta[c] + theta[c]^2)/(1-rho[c]^2)))
        tau.ar[c] <- pow(sigma.ar[c], -2)
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

        # hierarchy for thetas
        mu.theta[r] ~ dnorm(mu.theta.global, tau.theta.global)
        tau.theta[r] <- pow(sigma.theta[r], -2)
        sigma.theta[r] ~ dunif(0, 40)

        # hierarchy for smoother
        chi.eta[r] ~ dnorm(mu.chi.global, tau.chi.global)
        tau.eta[r] <- pow(sigma.eta[r], -2)
        sigma.eta[r] ~ dunif(0, 40)

        } # end r
        mu.beta.global ~ dnorm(0, 0.01)
        tau.beta.global <- pow(sigma.beta.global, -2)
        sigma.beta.global ~ dunif(0, 40)

        mu.rho.global ~ dunif(0,1)
        tau.rho.global <- pow(sigma.rho.global, -2)
        sigma.rho.global ~ dunif(0, 40)

        mu.theta.global ~ dunif(-1,0)
        tau.theta.global <- pow(sigma.theta.global, -2)
        sigma.theta.global ~ dunif(0, 40)

        mu.chi.global ~ dnorm(0, 0.01)
        tau.chi.global <- pow(sigma.chi.global, -2)
        sigma.chi.global ~ dunif(0, 40)


        ", file = file.path(file.name), fill = T, append = T)
        }
  if(cs.arma&!cs.smoothing){
    cat("
        mu.ct[c,1] ~ dnorm(beta[c], 1/eta)
        e.ct[c,1] ~ dnorm(sigma.ar[c]^2/eta*mu.ct[c,1], 1/(sigma.ar[c]^2*(1-sigma.ar[c]^2/eta)))
        for (t in 2:nyears.c[c]){
        e.ct[c,t] ~ dnorm(0, tau.ar[c])
        mu.ct[c,t] <- rho[c]*mu.ct[c,t-1] - theta[c]*e.ct[c,t-1] + e.ct[c,t]
        } # end t
        rho[c] ~ dnorm(mu.rho[region.c[c]], tau.rho[region.c[c]])T(0,1)
        theta[c] ~ dnorm(mu.theta[region.c[c]], tau.theta[region.c[c]])T(-1,0)
        beta[c] ~ dnorm(mu.beta[region.c[c]],tau.beta[region.c[c]])
        sigma.ar[c] <- sqrt(eta/((1-2*rho[c]*theta[c] + theta[c]^2)/(1-rho[c]^2)))
        tau.ar[c] <- pow(sigma.ar[c], -2)
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

        # hierarchy for thetas
        mu.theta[r] ~ dnorm(mu.theta.global, tau.theta.global)
        tau.theta[r] <- pow(sigma.theta[r], -2)
        sigma.theta[r] ~ dunif(0, 40)

        } # end r
        mu.beta.global ~ dnorm(0, 0.01)
        tau.beta.global <- pow(sigma.beta.global, -2)
        sigma.beta.global ~ dunif(0, 40)

        mu.rho.global ~ dunif(0,1)
        tau.rho.global <- pow(sigma.rho.global, -2)
        sigma.rho.global ~ dunif(0, 40)

        mu.theta.global ~ dunif(-1,0)
        tau.theta.global <- pow(sigma.theta.global, -2)
        sigma.theta.global ~ dunif(0, 40)

        eta <- pow(sqrteta, 2)
        sqrteta ~ dunif(0, 40)
        ", file = file.path(file.name), fill = T, append = T)
  }
  if(!cs.arma&cs.smoothing){
    cat("
        mu.ct[c,1] ~ dnorm(beta[c], 1/eta[c])
        e.ct[c,1] ~ dnorm(sigma.ar[c]^2/eta[c]*mu.ct[c,1], 1/(sigma.ar[c]^2*(1-sigma.ar[c]^2/eta[c])))
        for (t in 2:nyears.c[c]){
        e.ct[c,t] ~ dnorm(0, tau.ar[c])
        mu.ct[c,t] <- rho[c]*mu.ct[c,t-1] - theta[c]*e.ct[c,t-1] + e.ct[c,t]
        } # end t
        rho[c] ~ dnorm(mu.rho.global, tau.rho.global)T(0,1)
        theta[c] ~ dnorm(mu.theta.global, tau.theta.global)T(-1,0)
        eta[c] <- pow(sqrteta[c],2)
        sqrteta[c] <- exp(logsqrteta[c])
        logsqrteta[c] ~ dnorm(chi.eta[region.c[c]], tau.eta[region.c[c]])
        beta[c] ~ dnorm(mu.beta.global,tau.beta.global)
        sigma.ar[c] <- sqrt(eta[c] /((1-2*rho[c]*theta[c] + theta[c]^2)/(1-rho[c]^2)))
        tau.ar[c] <- pow(sigma.ar[c], -2)
  } # end c
        for(r in 1:nregions){
        # hierarchy for smoother
        chi.eta[r] ~ dnorm(mu.chi.global, tau.chi.global)
        tau.eta[r] <- pow(sigma.eta[r], -2)
        sigma.eta[r] ~ dunif(0, 40)

        } # end r
        mu.beta.global ~ dnorm(0, 0.01)
        tau.beta.global <- pow(sigma.beta.global, -2)
        sigma.beta.global ~ dunif(0, 40)

        mu.rho.global ~ dunif(0,1)
        tau.rho.global <- pow(sigma.rho.global, -2)
        sigma.rho.global ~ dunif(0, 40)

        mu.theta.global ~ dunif(-1,0)
        tau.theta.global <- pow(sigma.theta.global, -2)
        sigma.theta.global ~ dunif(0, 40)

        mu.chi.global ~ dnorm(0, 0.01)
        tau.chi.global <- pow(sigma.chi.global, -2)
        sigma.chi.global ~ dunif(0, 40)

        ", file = file.path(file.name), fill = T, append = T)
}
  if(!cs.arma&!cs.smoothing){
    cat("
        mu.ct[c,1] ~ dnorm(beta[c], 1/eta)
        e.ct[c,1] ~ dnorm(sigma.ar[c]^2/eta*mu.ct[c,1], 1/(sigma.ar[c]^2*(1-sigma.ar[c]^2/eta)))
        for (t in 2:nyears.c[c]){
        e.ct[c,t] ~ dnorm(0, tau.ar[c])
        mu.ct[c,t] <- rho[c]*mu.ct[c,t-1] - theta[c]*e.ct[c,t-1] + e.ct[c,t]
        } # end t
        rho[c] ~ dnorm(mu.rho.global, tau.rho.global)T(0,1)
        theta[c] ~ dnorm(mu.theta.global, tau.theta.global)T(-1,0)
        beta[c] ~ dnorm(mu.beta.global,tau.beta.global)
        sigma.ar[c] <- sqrt(eta /((1-2*rho[c]*theta[c] + theta[c]^2)/(1-rho[c]^2)))
        tau.ar[c] <- pow(sigma.ar[c], -2)
  } # end c
        mu.beta.global ~ dnorm(0, 0.01)
        tau.beta.global <- pow(sigma.beta.global, -2)
        sigma.beta.global ~ dunif(0, 40)

        mu.rho.global ~ dunif(0,1)
        tau.rho.global <- pow(sigma.rho.global, -2)
        sigma.rho.global ~ dunif(0, 40)

        mu.theta.global ~ dunif(-1,0)
        tau.theta.global <- pow(sigma.theta.global, -2)
        sigma.theta.global ~ dunif(0, 40)

        eta <- pow(sqrteta, 2)
        sqrteta ~ dunif(0, 40)
        ", file = file.path(file.name), fill = T, append = T)
    }
  if(time.trend){
    cat("
        for(c in 1:n.iso){
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

