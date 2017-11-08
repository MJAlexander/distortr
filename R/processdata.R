#' Process country data to run JAGS model
#'
#' Function to process a country data to run a hierarchical model in JAGS.
#' The function takes a dataframe with country codes, region codes, observations, observation years, standard errors and data source
#' and returns data in form for input into JAGS.
#'
#' @param d a dataframe that contains country codes/names, region codes/names, data column, observation years and standard error column
#' @param iso.column string specifying name of column which contains country codes/names
#' @param data.column string specifying name of column which contains observations
#' @param se.column string specifying name of column which contains standard errors
#' @param obsyear.column string specifying name of column which contains observation years
#' @param region.column string specifying name of column which contains region codes/names. If \code{NULL}, the number of regions is 1 i.e. region is assumed to be the World.
#' @param source.column string specifying name of column which contains data source. If \code{NULL}, the number of sources is assumed to be 1.
#' @param end.year final year for which estimates are required. Default is 2015.
#' @export
#' @return A list which includes matrices of data, observation years, source type and standard errors of dimensions \code{c x i}
#' where \code{c} is the number of countries and \code{i} is the maximum number of observations. The list also includes counts of countries, regions, source types, years by country
#' and starting years of each country, which are required to run the JAGS model.


processData <- function(d,
                        iso.column,
                        data.column,
                        se.column,
                        obsyear.column,
                        region.column = NULL, # if NULL, add a column of 1s
                        source.column = NULL, # if NULL, add a column of 1s
                        end.year = 2015
                        ){

  # number of observations per country
  d <- d %>% group_by_(iso.column) %>% mutate(n_obs = n()) %>% mutate_(n_na = paste0("sum(is.na(", data.column, "))"))
  # check to see if there are any countries with no data
  if(sum(d$n_obs==d$n_na)) stop("Please remove countries with no data.\n ")

  # countries and number of countries
  isos <- unique(d[,iso.column])[[1]]
  niso <- length(isos)

  # number of observations per country
  n.c <- (d[ !duplicated(d[, iso.column]), "n_obs"])[[1]]

  # region of country
  if(is.null(region.column)){
    region.c <- rep(1, niso)
    region.lookup <- data.frame(region = "World", region_number = 1)
    nregions <- 1
  }
  else{
    region.c <- as.numeric(as.factor((d[!duplicated(d[,iso.column]), region.column])[[1]]))
    regions <- levels(as.factor((d[!duplicated(d[,iso.column]), region.column])[[1]]))
    region.lookup <- data.frame(region = regions, region_number = 1:6) # for reference
    nregions <- length(unique(region.c))
  }

  # number of sources
  if(is.null(source.column)){
    nsources <- 1
  }
  else{
    nsources <- length(unique(d[,source.column])[[1]])
  }

  ## make c x i matrices of data, ses, observation years, source type
  y.ci <- matrix(NA, ncol = max(d$n_obs), nrow = niso)
  se.ci <- matrix(NA, ncol = max(d$n_obs), nrow = niso)
  gett.ci <- matrix(NA, ncol = max(d$n_obs), nrow = niso)
  source.ci <- matrix(NA, ncol = max(d$n_obs), nrow = niso)

  for(i in 1:niso){
    y.ci[i,1:d$n_obs[d[,iso.column] == isos[i]][1]] <- d[d[,iso.column]==isos[i], data.column][[1]]
    se.ci[i,1:d$n_obs[d[,iso.column] == isos[i]][1]] <- d[d[,iso.column]==isos[i], se.column][[1]]
    gett.ci[i,1:d$n_obs[d[,iso.column] == isos[i]][1]] <- d[d[,iso.column]==isos[i], obsyear.column][[1]]
    if(is.null(source.column)){
      source.ci[i,1:d$n_obs[d[,iso.column] == isos[i]][1]] <- 1
    }
    else{
      source.ci[i,1:d$n_obs[d[,iso.column] == isos[i]][1]] <- d[d[,iso.column]==isos[i], source.column][[1]]
    }
  }

  # start year by country, nyears of observations by country
  startyear.c <- apply(gett.ci, 1, min, na.rm=T)
  nyears.c <- sapply(1:length(startyear.c), function (i) length(startyear.c[i]:end.year))

  return(list(y.ci = y.ci,
              se.ci = se.ci,
              gett.ci = gett.ci,
              source.ci = source.ci,
              isos = isos,
              niso = niso,
              n.c = n.c,
              region.c = region.c,
              nregions = nregions,
              region.lookup = region.lookup,
              nsources = nsources,
              startyear.c = startyear.c,
              nyears.c = nyears.c
          ))
}
