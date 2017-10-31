# distortr
This package aims to aid in the exploration and fitting of temporal smoothing methods. It was designed to be particularly useful in the context of exploring models to estimate demographic indicators. Given data are often sparse or unreliable, especially in the case of developing countries, models that estimate demographic indicators for multiple areas/countries are often hierarchical, incorporating pooling of information across geographies. `distortr` allows for these types of Bayesian hirearchical models to be fit with a variety of different temporal smoothers. 

The package consists of two main parts:

1. Functions to simulate time series of distortions and fit and validate models on simulated data.  
2. Functions to fit Bayesian hierarchical models to datasets with observations from multiple countries. The user can specify whether or not to include a linear trend, and the type of temporal smoother to fit to the data. 

Temporal smoothing methods available for simulation and model fitting are:
* AR(1)
* ARMA(1,1)
* Splines regression (unpenalized, first-order penalized, second-order penalized)
* Gaussian process regression (squared exponential or Matern covariance functions)

## How to install
In R:
```R
# install
library(tidyverse)
devtools::install_github("MJAlexander/distortr")
library(distortr)
```


## Generating simulated data

Simulated time series of data can be created from any of the processes listed above (AR(1), ARMA(1,1), Splines or Gaussian Process). The various parameters associated with each function can be specified, and then the time series are simulated using the `simulateFluctuations` function. This returns a dataframe with `x` and `y` values. The user can also specify how much of the time series is missing. The sample autocorrelation function of the `y` values can be plotted using `plotACF`.

## Using real data

Datasets with observations from multiple countries/areas can be used, with the following columns requred:

- country/area name or code (e.g. country iso code)
- value of observation
- year of observation
- sampling error of observation
- data source (e.g. survey, administrative)

In addition, a region column may also be included (e.g. World Bank region). By default the builtin models include a region heirarchy (i.e. a country within a region within the world). However, models can also be run without the region level. 

Data can be visualized using the `plotData` function.

## Model fitting

Bayesian hierarchical models can be fitted to both simulated and real data using the JAGS software. Any model from the above list of processes can be chosen. Models are fitted using the `runMCMC` function. Results are obtained and plotted using the `getResults` and `plotResults` functions. 

## Validation

For simulated data, model performance can be validated using the `runModelValidation` function. This function leaves out a certain number of observations (default is 20%), either at random or the most recent observations. Once the observations are removed, the desired model is rerun and several validation measures are calculated:

* root-mean-squared error
* coverage of uncertainty intervals
* sharpness of uncertainty intervals (interval score)

If the most recent observations are left-out, the results are also returned, which can be plotted.

## Model Selection

For real data, WAIC can be calculated using the `waic` function. 

## Examples

For an example using real data, please refer to the file [real_data_anc4_example.R](./real_data_anc4_example.R). The raw data used in this example can be found in the [data folder](./data/).

An example workflow using simulated data:

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

# Plot sample ACF function
plotACF(df, nyears)

# Fit first-order penalized splines 
mod <- runMCMC(input.data = df, nyears = 100, method = "splines", order = 1, obs.err = T,
              nchains = 4, nburnin = 100, niter = 100+3000, nthin = 3)

df.mu <- getResults(mod, method = "splines", nyears=100)
plotResults(df, df.mu, method = "splines", order = 1, 
            maintitle = "AR(1) data with splines fit", save.plot = F)

# Run a validation, leaving out most recent 20% of data
rs <- runModelValidation(input.data = df, nyears = 100, method = "splines", order = 1,
                         nchains = 4, nburnin = 100, niter = 100+3000, nthin = 3, 
                         leave.out.method = "recent", leave.out.percent = 20)

plotResults(df, rs$df.res, rs$df.res.full, method = "splines", order = 1, 
            maintitle = "Splines validation", save.plot = F)

```
