############# distortr package example code #############
################ real data set (ANC4) ###################


# 0. Load packages --------------------------------------------------------

# if not already installed, to install distortr:
# devtools::install_github("MJAlexander/distortr")
library(distortr)
# these packages will need to be installed using install.packages if not already installed.
library(rjags)
library(R2jags)
library(fields)
library(splines)
library(boot)
library(RColorBrewer)
library(tidyverse)

# 1. Read and process data ------------------------------------------------

# Example data is related to the measure ANC4,
# which is the percentage of women aged 15–49 years attended at least four times during pregnancy by any provider.

d <- anc4

# remove countries with no data
d <- d %>% group_by(iso) %>% mutate(n_obs = n(), n_na = sum(is.na(anc_prop)))
d$no_data <- d$n_obs==d$n_na
d_orig <- d
d <- d[d$no_data==FALSE,]

# impute missing standard errors.
d <- imputeSE(d = d, iso.column = "iso", data.column = "anc_prop", se.column = "se", transform.se.column = "logit_se")
# this adds columns 'se_final' and 'se_imputed'
# d %>% dplyr::select(iso, anc_prop, logit_se, se_final, tr_se_final, se_imputed)

# list of country codes
isos <- unique(d[["iso"]])

# process data to be in the form neccesary for JAGS models.
# The main purpose of this is to change data, standard errors, years etc
# into matrices of dimension c x i where c is country and i is obervation number

## NOTE: by default, models have a region hierarchy.
# If this is not required, then make a new column d$global_region <- 1
# and use this column as the region column below.

input.data <- processData(d = d, iso.column = "iso",
                         data.column = "logit_prop",
                         se.column = "tr_se_final",
                         obsyear.column = "obs_year",
                         region.column = "who_region",
                         source.column = "source_type",
                         end.year = 2015)

# extract data from list

for(i in 1:length(input.data)){
  assign(names(input.data)[i], input.data[[i]])
}




# 2. Plot data ------------------------------------------------------------

# example, plot data for Brazil
iso <- "BRA"

p <- plotData(data.df = d[d$iso==iso,c("obs_year", "logit_prop", "se_final", "source")] , plot.se = T)
p + ylab("ANC4 proportion") + xlab("Year") + ggtitle(paste("Data for", iso))




# 3. Run model ------------------------------------------------------------

## Four method options: ar, arma, splines, gp
## For all methods, have to decide whether
## 1) to have a country-specific smoothing parameter (as opposed to global)
## 2) whether to estimate or input non-sampling error

cs.smoothing <- TRUE
nserror.estimated <- TRUE

# if non-sampling error is to be inputted, sigma.y must be a vector
# of values of non-sampling error for each source type.
# for example, there are three source types with ANC4 (survey, admin, other), so could put
# input.data$sigma.y <- c(0.2, 0.1, 0.5)
# if non-sampling error is to be estimated, sigma.y can just be NA
input.data$sigma.y <- NA

## JAGS parameters (can be changed if not converging - set to be relatively fast)
nchains = 3
nburnin = 1000
niter = 2000
nthin = 1

####### below is an example of each model

# 3a. AR(1) ---------------------------------------------------------------

## For AR(1) method, need to decide
## 1) whether to include time trend (linear)
## 2) whether to have country-specific AR parameters

time.trend <- FALSE
cs.arma <- TRUE

# run the model!
mod_ar <- runMCMC(method = "ar",
                  input.data = input.data,
                  cs.arma = cs.arma, cs.smoothing = cs.smoothing,
                  time.trend = time.trend, nserror.estimated = nserror.estimated,
                  nchains = nchains,
                  nburnin = nburnin,
                  niter = niter,
                  nthin = nthin,
                  model.save.file.path = "model_ar.txt") # the model is saved to a txt file


# 3b. ARMA(1,1) -----------------------------------------------------------

## For ARMA(1,1) method, need to decide
## 1) whether to include time trend (linear)
## 2) whether to have country-specific ARMA parameters

time.trend <- FALSE
cs.arma <- TRUE

# run the model!
mod_arma <- runMCMC(method = "arma",
                  input.data = input.data,
                  cs.arma = cs.arma, cs.smoothing = cs.smoothing,
                  time.trend = time.trend, nserror.estimated = nserror.estimated,
                  nchains = nchains,
                  nburnin = nburnin,
                  niter = niter,
                  nthin = nthin,
                  model.save.file.path = "model_arma.txt")


# 3c. Splines order 1 -----------------------------------------------------

# P-splines with first order penalization has intercept but no time trend
# need additional input data

splines.data <- getSplinesData(nyears.c, niso, order = 1)
input.data.splines <- c(input.data, splines.data)

# note need to specify order
mod_splines1 <- runMCMC(method = "splines",
                    order = 1,
                    input.data = input.data.splines,
                    cs.smoothing = cs.smoothing,
                    nserror.estimated = nserror.estimated,
                    nchains = nchains,
                    nburnin = nburnin,
                    niter = niter,
                    nthin = nthin,
                    model.save.file.path = "model_splines1.txt")

# 3d. Splines order 2 -----------------------------------------------------

# P-splines with first order penalization has intercept and time trend
# need additional input data

splines.data <- getSplinesData(nyears.c, niso, order = 2)
input.data.splines <- c(input.data, splines.data)

# note need to specify order
mod_splines2 <- runMCMC(method = "splines",
                        order = 2,
                        input.data = input.data.splines,
                        cs.smoothing = cs.smoothing,
                        nserror.estimated = nserror.estimated,
                        nchains = nchains,
                        nburnin = nburnin,
                        niter = niter,
                        nthin = nthin,
                        model.save.file.path = "model_splines2.txt")

# 3e. Gaussian Process - squared exponential -----------------------------------------------------

## For GP method, need to decide
## 1) whether to include time trend (linear)
## 2) covariance function (squared exponential or matern)

time.trend <- FALSE

# need extra data related to correlation matrix for GP model
Dist.c <- getGPData(nyear.c, niso, cov.method = "sqexp")
input.data$Dist.c <- Dist.c

mod_gpex <- runMCMC(method = "gp",
                  matern.cov = FALSE,
                  input.data = input.data,
                  cs.smoothing = cs.smoothing,
                  time.trend = time.trend,
                  nserror.estimated = nserror.estimated,
                  nchains = nchains,
                  nburnin = nburnin,
                  niter = niter,
                  nthin = nthin,
                  model.save.file.path = "model_gpex.txt")

# 3e. Gaussian Process - Matern-----------------------------------------------------

## For GP method, need to decide
## 1) whether to include time trend (linear)
## 2) covariance function (squared exponential or matern)

time.trend <- FALSE

## for Matern, need to decide
## 1) range (affects how fast decay of correlation structure is - smaller numbers have faster decay)
## 2) smoothness (differentiability of covariance function)
range <- 10
smoothness <- 2

# need extra data related to correlation matrix for GP model
Sigma.corr.c <- getGPData(nyear.c, niso, cov.method = "matern",
                          range = range, smoothness = smoothness)
input.data$Sigma.corr.c <- Sigma.corr.c

mod_gpmatern <- runMCMC(method = "gp",
                    matern.cov = TRUE,
                    input.data = input.data,
                    cs.smoothing = cs.smoothing,
                    time.trend = time.trend,
                    nserror.estimated = nserror.estimated,
                    nchains = nchains,
                    nburnin = nburnin,
                    niter = niter,
                    nthin = nthin,
                    model.save.file.path = "model_gpmatern.txt")






# 4. Calculate and plot results -------------------------------------------

## for example, for Brazil and AR(1)
iso <- "BRA"
iso.number <- which(isos==iso)

# calculate results
df.res <- getResults(mod_ar, method = "ar", iso.number = iso.number,
                    nyears = nyears.c[iso.number],
                    startyear = startyear.c[iso.number])
# need to transform back from logit scale
df.res$lower <- inv.logit(df.res$lower)
df.res$median <- inv.logit(df.res$median)
df.res$upper <- inv.logit(df.res$upper)

# plot results
p <- plotResults(data.df = d[d$iso==iso,c("obs_year", "anc_prop", "se_final", "source")],
                 res.df = df.res, plot.se=T,
                 save.plot = F, maintitle = isos[iso.number])
pf <- p+ ylim(c(0,1))+ ylab("ANC4 Proportion") + xlab("Year")
print(pf)

# could do this for all countries and save via a loop
# eg. for first order splines
pdf("spline1_results.pdf")
for(i in 1:niso){
  # calculate results
  df.res <- getResults(mod_splines1, method = "splines", iso.number = i,
                       nyears = nyears.c[i],
                       startyear = startyear.c[i])
  df.res$lower <- inv.logit(df.res$lower)
  df.res$median <- inv.logit(df.res$median)
  df.res$upper <- inv.logit(df.res$upper)

  # plot results
  p <- plotResults(data.df = d[d$iso==isos[i],c("obs_year", "anc_prop", "se_final", "source")],
                   res.df = df.res, plot.se=T,
                   save.plot = F, maintitle = isos[i])
  pf <- p+ ylim(c(0,1))+ ylab("ANC4 Proportion") + xlab("Year")
  print(pf)

}
dev.off()

# 5. Calculate WAIC -------------------------------------------------------

# WAIC to compare fit of different models
ll.ar.ts <- t(do.call("cbind", lapply(1:dim(mod_ar$BUGSoutput$sims.list[["loglike.ci"]])[2],
                                   function(i) mod_ar$BUGSoutput$sims.list[["loglike.ci"]][,i,])))
ll.arma.ts <- t(do.call("cbind", lapply(1:dim(mod_arma$BUGSoutput$sims.list[["loglike.ci"]])[2],
                                      function(i) mod_arma$BUGSoutput$sims.list[["loglike.ci"]][,i,])))
ll.spline1.ts <- t(do.call("cbind", lapply(1:dim(mod_splines1$BUGSoutput$sims.list[["loglike.ci"]])[2],
                                        function(i) mod_splines1$BUGSoutput$sims.list[["loglike.ci"]][,i,])))

waic_ar <- waic(log_lik = ll.ar.ts)$total['waic']; waic_ar
waic_arma <- waic(log_lik = ll.arma.ts)$total['waic']; waic_arma
waic_splines1 <- waic(log_lik = ll.spline1.ts)$total['waic']; waic_splines1
