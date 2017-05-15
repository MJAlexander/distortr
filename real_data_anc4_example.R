############# distortr package example code #############
################ real data set (ANC4) ###################


# 0. Load packages --------------------------------------------------------


# load other packages

library(rjags)
library(R2jags)
library(fields)
library(splines)
library(boot)
library(RColorBrewer)
library(tidyverse)




# 1. Read and process data ------------------------------------------------

# read and clean raw data. NOTE: One observation (Brazil, 1999) is removed.
d <- cleanANCData(file.path = "data/who_rhr_anc4_detailed_2017.csv", save.file = TRUE)

# remove countries with no data
d <- d %>% group_by(iso) %>% mutate(n_obs = n(), n_na = sum(is.na(anc_prop)))
d$no_data <- d$n_obs==d$n_na
d_orig <- d
d <- d[d$no_data==FALSE,]

# impute missing standard errors.
d <- imputeSE(d = d, iso.column = "iso", data.column = "anc_prop", se.column = "logit_se")
# this adds columns 'se_final' and 'se_imputed'
# d %>% dplyr::select(iso, anc_prop, logit_se, se_final, se_imputed)

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
                         se.column = "se_final",
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

