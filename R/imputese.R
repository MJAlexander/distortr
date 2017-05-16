#' Impute missing standard errors
#'
#' Impute missing standard errors by country based on other standard errors in that country.
#' Imputed value is assumed to be the maximum standard error observed in that country multiplied by a specificed factor.
#' If all standard errors within a country are missing, a pre-specified value can be imputed.
#'
#' @param d a dataframe that contains country codes/names, data column and standard error column
#' @param iso.column string specifying name of column which contains country codes/names
#' @param data.column string specifying name of column which contains observations
#' @param se.column string specifying name of column which contains standard errors
#' @param mult.factor multiplication factor to multiply the maximum observed standard error
#' @param impute.value value of standard error to impute if no information available in country.
#' @export
#' @return A data frame of final standard error values and a column indicating whether the standard error value was imputed.



imputeSE <- function(d,
                     iso.column,
                     data.column,
                     se.column,
                     transform.se.column = NULL,
                     mult.factor = 1.5,
                     impute.value = 0.025,
                     transform.impute.value = 0.1){
  isos <- unique(d[[iso.column]])
  # fill in miss SEs based on max within each country. If that does not exist, choose an imputation value
  d$se_final <- d[[se.column]]
  d$tr_se_final <- d[[transform.se.column]]
  d$se_imputed <- 0 # binary
  for(i in 1:length(isos)){
    if(sum(d[d$iso==isos[i],data.column], na.rm = T)!=0){ #if there are data
      if(sum(is.na(d[d$iso == isos[i],se.column]))!=0){ # if there are ses missing
        if(sum(is.na(d[d$iso == isos[i],se.column]))<length(d[d$iso == isos[i],se.column][[1]])){ #if there are some observed SEs
          d[d$iso == isos[i]&is.na(d$se),"se_final"] <- max(d[d$iso == isos[i],se.column], na.rm=T)*mult.factor
          d[d$iso == isos[i]&is.na(d$se),"tr_se_final"] <- max(d[d$iso == isos[i],transform.se.column], na.rm=T)*mult.factor
          d[d$iso == isos[i]&is.na(d$se),"se_imputed"] <- 1
        }
        else{ # fill with pre-chosen value
          d[d$iso == isos[i]&is.na(d$se),"se_final"] <- impute.value
          d[d$iso == isos[i]&is.na(d$se),"tr_se_final"] <- transform.impute.value
          d[d$iso == isos[i]&is.na(d$se),"se_imputed"] <- 1
        }
      }
    }
  }
  d$se_final[d$se_final==0] <- 0.001
  return(d)
}
