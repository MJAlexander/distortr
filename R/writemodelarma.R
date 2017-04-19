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
#' cs.rho <- T
#' cs.smoothing <- T
#' time.trend <- T
#' nserror.estimated <- T
#' writeModelAR(cs.rho = cs.rho, cs.smoothing = cs.smoothing, time.trend = time.trend, nserror.estimated = nserror.estimated)

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
        nu.ci[c,i] <- pow((se.ci[c,i]^2+sigma.y^2), -1)
        }
        ", file = file.path(file.name), fill = T, append = T)
        }
  if(cs.rho&cs.smoothing){
    cat("
        mu.ct[c,1] ~ dnorm(beta[c], 1/eta[c])
        e.ct[c,1] ~ dnorm(sigma.ar[c]^2/eta[c]*mu.ct[c,1], 1/(sigma.ar[c]^2*(1-sigma.ar[c]^2/eta[c])))
        for (t in 2:nyears.c[c]){
        e.ct[c,t] ~ dnorm(0, tau.ar[c])
        mu.ct[c,t] <- rho[c]*mu.ct[c,t-1] - theta[c]*e.ct[c,t-1] + e.ct[c,t]
        } # end t
        rho[c] ~ dnorm(mu.rho[region.c[c]], sigma.rho[region.c[c]])T(0,1)
        theta[c] ~ dnorm(mu.theta[region.c[c]], sigma.theta[region.c[c]])T(-1,0)
        eta[c] <- pow(sqrteta[c],2)
        sqrteta[c] <- exp(logsqrteta[c])
        logsqrteta[c] ~ dnorm(chi.eta[region.c[c]], psi.eta[region.c[c]])
        beta[c] ~ dnorm(mu.beta[region.c[c]],tau.beta[region.c[c]])
        sigma.ar[c] <- sqrt(eta[c] /((1-2*rho[c]*theta[c] + theta[c]^2)/(1-rho[c]^2)))
        tau.ar[c] <- pow(sigma.ar[c], -2)
  } # end c
        for(r in 1:nregions){
        # hierarchy for betas
        mu.beta[r] ~ dnorm(mu.beta.global, tau.beta.global)
        tau.beta[r] <- pow(sigma.beta[r], -2)
        sigma.beta[r] <- exp(logsigma.beta[r])
        logsigma.beta[r] ~ dnorm(chi.beta, psi.beta)

        # hierarchy for rhos
        mu.rho[r] ~ dnorm(mu.rho.global, tau.rho.global)T(0, 0.99)
        tau.rho[r] <- pow(sigma.rho[r], -2)
        sigma.rho[r] <- exp(logsigma.rho[r])
        logsigma.rho[r] ~ dnorm(chi.rho, psi.rho)

        # hierarchy for thetas
        mu.theta[r] ~ dnorm(mu.theta.global, tau.theta.global)T(-0.99, 0)
        tau.theta[r] <- pow(sigma.theta[r], -2)
        sigma.theta[r] <- exp(logsigma.theta[r])
        logsigma.theta[r] ~ dnorm(chi.theta, psi.theta)

        # hierarchy for smoother
        chi.eta[r] ~ dnorm(mu.chi.global, tau.chi.global)
        psi.eta[r] <- exp(logpsi.eta[r])
        logpsi.eta[r] ~ dnorm(chi.psi, psi.psi)

        } # end r
        mu.beta.global ~ dnorm(0, 0.01)
        tau.beta.global <- pow(sigma.beta.global, -2)
        sigma.beta.global ~ dunif(0, 40)
        chi.beta ~ dnorm(0, 0.01)
        sigma.psi.beta ~ dunif(0, 40)
        psi.beta <- pow(sigma.psi.beta, -2)

        mu.rho.global ~ dunif(0,1)
        tau.rho.global <- pow(sigma.rho.global, -2)
        sigma.rho.global ~ dunif(0, 40)
        chi.rho ~ dnorm(0, 0.01)
        sigma.psi.rho ~ dunif(0, 40)
        psi.rho <- pow(sigma.psi.rho, -2)

        mu.theta.global ~ dunif(-1,0)
        tau.theta.global <- pow(sigma.theta.global, -2)
        sigma.theta.global ~ dunif(0, 40)
        chi.theta ~ dnorm(0, 0.01)
        sigma.psi.theta ~ dunif(0, 40)
        psi.theta <- pow(sigma.psi.theta, -2)

        mu.chi.global ~ dnorm(0, 0.01)
        tau.chi.global <- pow(sigma.chi.global, -2)
        sigma.chi.global ~ dunif(0, 40)
        chi.psi ~ dnorm(0, 0.01)
        sigma.psi.psi ~ dunif(0, 40)
        psi.psi <- pow(sigma.psi.psi, -2)


        ", file = file.path(file.name), fill = T, append = T)
        }
  if(cs.rho&!cs.smoothing){
    cat("
        mu.ct[c,1] ~ dnorm(beta[c], tau.stat[c])
        tau.stat[c] <- (1-pow(rho[c],2))/pow(sigma,2)
        for (t in 2:nyears.c[c]){
        mu.ct[c,t] ~ dnorm(muhat.ct[c,t], tau)
        muhat.ct[c,t] <- rho[c]*mu.ct[c,t-1]
        } # end t
        rho[c] ~ dnorm(mu.rho[region.c[c]], sigma.rho[region.c[c]])T(-1,1)
        beta[c] ~ dnorm(mu.beta[region.c[c]],tau.beta[region.c[c]])
  } # end c
        for(r in 1:nregions){
        # hierarchy for betas
        mu.beta[r] ~ dnorm(mu.beta.global, tau.beta.global)
        tau.beta[r] <- pow(sigma.beta[r], -2)
        sigma.beta[r] <- exp(logsigma.beta[r])
        logsigma.beta[r] ~ dnorm(chi.beta, psi.beta)

        # hierarchy for rhos
        mu.rho[r] ~ dnorm(mu.rho.global, tau.rho.global)
        tau.rho[r] <- pow(sigma.rho[r], -2)
        sigma.rho[r] <- exp(logsigma.rho[r])
        logsigma.rho[r] ~ dnorm(chi.rho, psi.rho)

        } # end r
        mu.beta.global ~ dnorm(0, 0.01)
        tau.beta.global <- pow(sigma.beta.global, -2)
        sigma.beta.global ~ dunif(0, 40)
        chi.beta ~ dnorm(0, 0.01)
        sigma.psi.beta ~ dunif(0, 40)
        psi.beta <- pow(sigma.psi.beta, -2)

        mu.rho.global ~ dnorm(0, 0.01)
        tau.rho.global <- pow(sigma.rho.global, -2)
        sigma.rho.global ~ dunif(0, 40)
        chi.rho ~ dnorm(0, 0.01)
        sigma.psi.rho ~ dunif(0, 40)
        psi.rho <- pow(sigma.psi.rho, -2)

        tau <- pow(sigma, -2)
        sigma ~ dunif(0, 40)
        ", file = file.path(file.name), fill = T, append = T)
  }
  if(!cs.rho&cs.smoothing){
    cat("
        mu.ct[c,1] ~ dnorm(beta[c], tau.stat[c])
        tau.stat[c] <- (1-pow(rho,2))/pow(sigma[c],2)
        for (t in 2:nyears.c[c]){
        mu.ct[c,t] ~ dnorm(muhat.ct[c,t], tau[c])
        muhat.ct[c,t] <- rho*mu.ct[c,t-1]
        } # end t
        tau[c] <- pow(sigma[c],-2)
        sigma[c] <- exp(logsigma[c])
        logsigma[c] ~ dnorm(chi, psi)
        beta[c] ~ dnorm(mu.beta,tau.beta)
  } # end c
        for(r in 1:nregions){
        # hierarchy for betas
        mu.beta[r] ~ dnorm(mu.beta.global, tau.beta.global)
        tau.beta[r] <- pow(sigma.beta[r], -2)
        sigma.beta[r] <- exp(logsigma.beta[r])
        logsigma.beta[r] ~ dnorm(chi.beta, psi.beta)

        # hierarchy for smoother
        chi.sigma[r] ~ dnorm(mu.chi.global, tau.chi.global)
        psi.sigma[r] <- exp(logpsi.sigma[r])
        logpsi.sigma[r] ~ dnorm(chi.psi, psi.psi)

        } # end r
        mu.beta.global ~ dnorm(0, 0.01)
        tau.beta.global <- pow(sigma.beta.global, -2)
        sigma.beta.global ~ dunif(0, 40)
        chi.beta ~ dnorm(0, 0.01)
        sigma.psi.beta ~ dunif(0, 40)
        psi.beta <- pow(sigma.psi.beta, -2)

        mu.chi.global ~ dnorm(0, 0.01)
        tau.chi.global <- pow(sigma.chi.global, -2)
        sigma.chi.global ~ dunif(0, 40)
        chi.psi ~ dnorm(0, 0.01)
        sigma.psi.psi ~ dunif(0, 40)
        psi.psi <- pow(sigma.psi.psi, -2)

        rho ~ dunif(-1,1)
        chi ~ dnorm(0, 0.01)
        ", file = file.path(file.name), fill = T, append = T)
}
  if(!cs.rho&!cs.smoothing){
    cat("
        mu.ct[c,1] ~ dnorm(beta[c], tau.stat)
        for (t in 2:nyears.c[c]){
        mu.ct[c,t] ~ dnorm(muhat.ct[c,t], tau)
        muhat.ct[c,t] <- rho*mu.ct[c,t-1]
        } # end t
        beta[c] ~ dnorm(mu.beta,tau.beta)
  } # end c
        mu.beta ~ dnorm(0, 0.01)
        tau.beta <- pow(sigma.beta, -2)
        sigma.beta ~ dunif(0, 40)
        tau.stat <- (1-pow(rho,2))/pow(sigma,2)
        rho ~ dunif(-1,1)
        tau <- pow(sigma, -2)
        sigma ~ dunif(0, 40)
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
        sigma.gamma[r] <- exp(logsigma.gamma[r])
        logsigma.gamma[r] ~ dnorm(chi.gamma, psi.gamma)
        } # end r
        mu.gamma.global ~ dnorm(0, 0.01)
        tau.gamma.global <- pow(sigma.gamma.global, -2)
        sigma.gamma.global ~ dunif(0, 40)
        chi.gamma ~ dnorm(0, 0.01)
        sigma.psi.gamma ~ dunif(0, 40)
        psi.gamma <- pow(sigma.psi.gamma, -2)
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

