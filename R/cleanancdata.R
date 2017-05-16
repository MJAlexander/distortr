#' Clean example dataset
#'
#' Clean the example dataset used to illustrate code functionality.
#' The dataset is from WHO and contains observations for the percentage of women aged 15–49 years attended at least four times during pregnancy by any provider.
#' Data are available for 150 countries.
#'
#' @param file.path file path to raw data file
#' @param save.file whether or not to save cleaned data. Default is \code{TRUE}.
#' @export
#' @return cleaned data set

cleanANCData <- function(file.path = "data/who_rhr_anc4_detailed_2017.csv",
                          save.file = TRUE){
  # read raw data
  d <- read_csv(file.path)

  # rename columns
  names(d) <- c("iso", "who_region", "country", "anc_percent", "see_notes",
                "coverage_start_year", "coverage_end_year", "survey_start_year", "survey_end_year",
                "sample_size", "source", "source_verified", "notes", "source_type", "age_if_not_15_49")

  # clean some variables to ensure correct format
  d$anc_percent <- as.numeric(d$anc_percent)
  d$sample_size[d$sample_size==0]<- NA
  d$sample_size <- as.numeric(d$sample_size)
  d$iso <- gsub(" ", "", d$iso)


  ## make some new variables to be used for modeling
  d$anc_prop <- d$anc_percent/100 # consider proportion rather than percentage
  d$se <- sqrt(d$anc_prop*(1-d$anc_prop)/d$sample_size) # sampling error
  # replace 0/1s with values slightly above/below to avoid Inf
  d$anc_prop[d$anc_prop==1] <- 0.9999
  d$anc_prop[d$anc_prop==0] <- 0.0001
  # define observation
  d$obs_year <- floor((d$coverage_start_year+d$coverage_end_year)/2)
  # source of data
  d$source <- factor(d$source_type, levels = 1:3, labels = c("survey", "admin", "other"))

  # we will model on the logit sacle
  d$logit_prop <- logit(d$anc_prop)
  d$logit_se <- 1/(d$anc_prop*(1-d$anc_prop))*d$se # by delta method

  # remove weird brazil point for now.
  d <- d[!(d$iso=="BRA"&d$obs_year==1999),]

  if(save.file){
    anc4 <- d
    dir.create(file.path("data/"), showWarnings = FALSE)
    save(anc4, file = "./data/anc4.RData")
  }
  return(d)
}
