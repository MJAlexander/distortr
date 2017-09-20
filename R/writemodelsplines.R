#' Write a JAGS model to fit penalized splines regression
#'
#' Write a JAGS model to fit first- and second-order penalized splines regression
#'
#' @param order order of penalization (either 1 or 2).
#' @param cs.smoothing whether smoothing paramter is country specific. If `FALSE`, smoothing parameter is global.
#' @param nserror.estimated whether to estimate non-sampling error. IF `FALSE`, fixed sampling error is inputted.
#' @param file.name name of file to be saved. Must be a `.txt` file
#' @export
#' @return A text file that contains a JAGS model
#' @examples
#' order <- 1
#' cs.smoothing <- T
#' nserror.estimated <- T
#' writeModelSplines(order = order, nserror.estimated = nserror.estimated, cs.smoothing = cs.smoothing)


writeModelSplines <- function(
  order = NULL,
  cs.smoothing = T,
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
      y.ci[c,i] ~ dnorm(mu.ct[c,gett.ci[c,i]], nu.ci[c,i])
      yrep.ci[c,i] ~ dnorm(mu.ct[c,gett.ci[c,i]], nu.ci[c,i]) #for validation
      loglike.ci[c,i] <- logdensity.norm(y.ci[c,i], mu.ct[c,gett.ci[c,i]], nu.ci[c,i]) #for WAIC
      nu.ci[c,i] <- pow((se.ci[c,i]^2+sigma.y[source.ci[c,i]]^2), -1)
      }
        ", file = file.path(file.name), fill = T, append = T)
  if(order==1){
    cat("
        # mean value
        for(t in 1:nyears.c[c]){
        mu.ct[c,t] <- beta.d[c] + inprod(Z.tkc[t,,c], delta.hc[,c])
        } # end t
        ", file = file.path(file.name), fill = T, append = T)
    if(cs.smoothing){
      cat("
          for (h in 1:H.c[c]){
            delta.hc[h,c] ~ dnorm(0, tau.delta[c])
          }
          for(h in (H.c[c]+1):H){
            delta.hc[h,c] <- 0
          }
          tau.delta[c] <- pow(sigma.delta[c], -2)
          logsigma.delta[c] ~ dnorm(chi.delta[region.c[c]], psi.delta[region.c[c]])
          sigma.delta[c]<- exp(logsigma.delta[c])
          beta.d[c] ~ dnorm(mu.beta[region.c[c]], tau.beta[region.c[c]])
        }
        for(r in 1:nregions){
          chi.delta[r] ~ dnorm(mu.chi.delta, tau.chi.delta)
          psi.delta[r] <- pow(sigma.psi.delta[r], -2)
          sigma.psi.delta[r] ~ dunif(0,40)
        }
        mu.chi.delta ~ dnorm(0, 0.01)
        tau.chi.delta <- pow(sigma.chi.delta, -2)
        sigma.chi.delta ~ dunif(0, 40)
          ", file = file.path(file.name), fill = T, append = T)
    }
    if(!cs.smoothing){
      cat("
          for (h in 1:H.c[c]){
          delta.hc[h,c] ~ dnorm(0, tau.delta)
          }
          for(h in (H.c[c]+1):H){
          delta.hc[h,c] <- 0
          }
            beta.d[c] ~ dnorm(mu.beta[region.c[c]], tau.beta[region.c[c]])
          }
          tau.delta <- pow(sigma.delta, -2)
          sigma.delta ~ dunif(0, 40)
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
  }
  if(order==2){
    cat("
        # mean value
        for(t in 1:nyears.c[c]){
          mu.ct[c,t] <- inprod(BG.tdc[t,,c], beta.d[,c]) + inprod(Z.tkc[t,,c], delta.hc[,c])
        } # end t
        ", file = file.path(file.name), fill = T, append = T)
    if(cs.smoothing){
      cat("
          for (h in 1:H.c[c]){
            delta.hc[h,c] ~ dnorm(0, tau.delta[c])
          }
          for(h in (H.c[c]+1):H){
            delta.hc[h,c] <- 0
          }
          tau.delta[c] <- pow(sigma.delta[c], -2)
          logsigma.delta[c] ~ dnorm(chi.delta[region.c[c]], psi.delta[region.c[c]])
          sigma.delta[c]<- exp(logsigma.delta[c])
          for(d in 1:D){
            beta.d[d,c] ~ dnorm(mu.beta[d, region.c[c]], tau.beta[d, region.c[c]])
          }
        }
        for(r in 1:nregions){
          chi.delta[r] ~ dnorm(mu.chi.delta, tau.chi.delta)
          psi.delta[r] <- pow(sigma.psi.delta[r], -2)
          sigma.psi.delta[r] ~ dunif(0,40)
          }
          mu.chi.delta ~ dnorm(0, 0.01)
          tau.chi.delta <- pow(sigma.chi.delta, -2)
          sigma.chi.delta ~ dunif(0, 40)

          ", file = file.path(file.name), fill = T, append = T)
}
    if(!cs.smoothing){
      cat("
          for (h in 1:H.c[c]){
          delta.hc[h,c] ~ dnorm(0, tau.delta)
          }
          for(h in (H.c[c]+1):H){
          delta.hc[h,c] <- 0
          }
          for(d in 1:D){
            beta.d[d,c] ~ dnorm(mu.beta[d, region.c[c]], tau.beta[d, region.c[c]])
          }
    }
          tau.delta <- pow(sigma.delta, -2)
          sigma.delta ~ dunif(0, 40)
          ", file = file.path(file.name), fill = T, append = T)
    }
    cat("
        for(d in 1:D){
          for(r in 1:nregions){
            mu.beta[d, r] ~ dnorm(mu.beta.global[d], tau.beta.global[d])
            tau.beta[d,r] <- pow(sigma.beta[d,r], -2)
            sigma.beta[d,r] ~ dunif(0, 40)
              }
          mu.beta.global[d] ~ dnorm(0, 0.01)
          tau.beta.global[d] <- pow(sigma.beta.global[d], -2)
          sigma.beta.global[d] ~ dunif(0, 40)
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
  cat(paste0("Model file written to ", getwd(),"/",file.name, "\n"))
  return(invisible())
}

