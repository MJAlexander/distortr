
runMCMCGlobal <- function(method,
                          input.data,
                          order = NULL,
                          matern.cov=TRUE,
                          cs.arma = NULL,
                          cs.smoothing = TRUE,
                          time.trend = FALSE,
                          nserror.estimated = TRUE,
                          nchains = 4,
                          nburnin = 1000,
                          niter = 1000+30000,
                          nthin = 30,
                          model.file.path = NULL,
                          model.save.file.path = "R/model.txt"){

  for(i in 1:length(input.data)){
    assign(names(input.data)[i], input.data[[i]])
  }

  if(method=="ar"){
    if(is.null(cs.arma)){
      stop("Need to specify cs.arma.")
    }
    if(is.null(model.file.path)){
      # write model based on desired parameters
      writeModelAR(cs.arma = cs.arma,
                   cs.smoothing = cs.smoothing,
                   time.trend = time.trend,
                   nserror.estimated = nserror.estimated,
                   file.name = model.save.file.path)
    }

    jags.data <- list(y.ci = y.ci, gett.ci = (gett.ci - startyear.c+1), niso = niso, n.c =n.c,
                      nyears.c= nyears.c, se.ci = se.ci, sigma.y = sigma.y,
                      source.ci = source.ci, nsources = nsources,
                      region.c = region.c, nregions = nregions)

    if(nserror.estimated){
      jags.data[["sigma.y"]] <- NULL
    }

    parnames <- c("mu.ct", "beta", "sigma", "rho", "sigma.y",
                  "mu.beta", "sigma.beta", "mu.beta.global", "sigma.beta.global")

    if(time.trend){
      parnames <- c(parnames, "gamma", "mu.gamma", "sigma.gamma", "mu.gamma.global", "sigma.gamma.global")
    }
    if(cs.arma){
      parnames <- c(parnames, "mu.rho", "sigma.rho", "mu.rho.global", "sigma.rho.global")
    }
    if(cs.smoothing){
      parnames <- c(parnames, "mu.logsigma", "sigma.logsigma", "mu.logsigma.global", "sigma.logsigma.global")
    }
  }

  # if(method=="arma"){
  # }
  #
  # if(method=="splines"){
  #   if(is.null(order)){
  #     stop("Order of penalization must be specified.")
  #   }
  # }
  #
  # if(method=="gp"){
  #   if(matern.cov==TRUE){
  #   }
  #   if(matern.cov==FALSE){
  #     }
  #   }
  # }

  if(is.null(model.file.path)){
    model.path.to.run <- model.save.file.path
  }
  else{
    model.path.to.run <- model.file.path
  }
  ## run the model
  cat("Running model.\n")
  mod <- jags.parallel(data = jags.data,
                       parameters.to.save= c(parnames),
                       n.chains = nchains,
                       n.burnin = nburnin,
                       n.iter = niter,
                       n.thin =nthin,
                       model.file = model.path.to.run)
  if(max(mod$BUGSoutput$summary[, c("Rhat")])>1.1) cat("Something hasn't converged.\n")
  return(mod)
}


