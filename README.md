# distortr
This package contains functions to simulate time series of distortions, and fit and validate models. Both the simulation and model fits can be based on one or a combination of ARMA, P-splines and Gaussian processes. 

## How to install
In R:
```R
# install
library(tidyverse)
devtools::install_github("MJAlexander/distortr")
library(distortr)
```


## Simulation

Simulated time series of data can be created following one of these processes:

* AR(1)
* ARMA(1,1)
* Splines regression (unpenalized, first-order penalized, second-order penalized)
* Gaussian process regression

The various parameters associated with each function can be specified, and then the time series are simulated using the `simulateFluctuations` function. This returns a dataframe with `x` and `y` values. 

## Model fitting

Bayesian hierarchical models can be fitted to time series of fluctuations using the JAGS software. Any model from the above list of processes can be chosen. Models are fitted using the `runMCMC` function. Results are obtained and plotted using the `getResults` and `plotResults` functions. 

## Validation

Model performance can be validated using the `runModelValidation` function. This function leaves out a certain number of observations (default is 20%), either at random or the most recent observations. Once the observations are removed, the desired model is rerun and several validation measures are calculated:

* root-mean-squared error
* coverage of uncertainty intervals
* sharpness of uncertainty intervals (interval score)

If the most recent observations are left-out, the results are also returned, which can be plotted.

## Examples

An example workflow:

```R
# check out help files
?simulateFluctuations
?GetAR

# Simulate AR(1) process over 100 years, with 80% of sample observed
# Set the parameters 
nyears <- 100
prop.sample <- 0.8
obs.err <- TRUE
sigma.y <- 0.5
seed <- 1234
method <- "ar"
params <- list(rho = 0.9, sigma = 0.25, ystart = NULL, eps0.t = NULL)
df <- simulateFluctuations(nyears, prop.sample, method, params, obs.err, sigma.y)

# Set measurement error to be 1
df$se <- 1

# Plot data
ggplot(data = df, aes(x = t, y = y)) + geom_point() + theme_bw()+ 
  geom_errorbar(data=df,aes(x=t,y=NULL,ymin=y-2*se, ymax=y+2*se), width=0.2) 

# Fit first-order penalized splines 
mod <- runMCMC(df = df, nyears = 100, method = "splines", order = 1, 
              nchains = 4, nburnin = 100, niter = 100+3000, nthin = 3)

df.mu <- getResults(mod, method = "splines")
plotResults(df, df.mu, method = "splines", order = 1, 
            maintitle = "AR(1) data with splines fit", save.plot = F)

# Run a validation, leaving out most recent 20% of data
rs <- runModelValidation(df = df, nyears = 100, method = "splines", order = 1,
                         nchains = 4, nburnin = 100, niter = 100+3000, nthin = 3, 
                         leave.out.method = "recent", leave.out.percent = 20)

plotResults(df, rs$df.res, rs$df.res.full, method = "splines", order = 1, 
            maintitle = "Splines validation", save.plot = F)

```
## Next steps

The model fitting and validation functions in this package will be modified in order to be fitted on real data. 
