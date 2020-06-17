library(data.table)
library(magrittr)

load("demo.RData")
data <- cohort_unrelated
group <- "Group"
time <- "Time"
subject <- "Subject"
response <- "Fixations"
responseGroup <- "LookType"
concave <-  TRUE 
rho.0 <-  0.9 
cor <-  TRUE 
cores <-  1
verbose <-  FALSE
curveType <- list('Cohort' = 'logistic', 'Unrelated_Cohort' = 'doubleGauss')

# factors >:(
ff <- sapply(data, is.factor)
data[ff] <- lapply(data[ff], as.character)
rm(ff)
data <- as.data.frame(data)

## possible additions
# refits - number of times to jitter initial parameters

## for each fit, give numeric value, i.e., 1 - AR1 R2 > 0.95, 2 - AR1 R2 > 0.85, etc. 
## this will make refitting SO much easier
## it will make dropping subjects easier too!
## Definitely going to be using data table for this


bdots(data, group, subject, time, response, responseGroup = "LookType")
bdots <- function(data, # dataset 
                  group, # groups for subjects
                  subject, # subjects
                  time, # column for time
                  response, # response vector
                  responseGroup = NULL, # response vector groups, i.e., looktype (this may not be there in general, i.e., tumr)
                  curveType = list('Cohort' = 'logistic', 'Unrelated_Cohort' = 'doubleGauss'), # logistic, doubleGauss, etc. maybe have length match responseGroup length?
                  concave = NULL, # doubleGauss only concavity
                  cor = TRUE, # autocorrelation?
                  rho.0 = 0.9, # autocor value
                  cores = 1, # cores to use
                  verbose = FALSE) {
  
  ## Need to validate inputs after things work ##
  ## Also kill all factors ##
  
  ## Cheat around DT reference, conditionally set key for subset
  dat <- data.table()
  dat$subject <- data[[subject]] %>% as.character()
  dat$time <- data[[time]] %>% as.numeric()
  dat$response <- data[[response]] %>% as.numeric()
  
  ## Assign group only if argument passed
  if(!is.null(group)) {
    dat$group <- data[[group]] %>% as.character()
  } else {
    dat$group <- "none"
  }
  
  ## Assign responseGroup only if argument passed
  if(!is.null(responseGroup)) {
    dat$responseGroup <- data[[responseGroup]] %>% as.character()
  } else {
    dat$responseGroup <- "none"
  }
  
  ## Set key
  setkeyv(dat, c("subject", "group", "responseGroup"))
  
  ## Loop through responseGroup in parellel
  respGroups <- unique(dat$responseGroup)
  
  ## response group is top level, as it's likely largest subset of data
  ## this is more efficient for passing in parallel
  ## This could be optimized if, say, cores > length(respGroup)
  ## in which case, pass parallel to subject groups
  ## Also need to go through and pass correct curve type
  if(length(respGroups) == 1) {
    fits <- bdots_fit(dat, curveType, concave, cor, rho.0, cores, verose)
  } else {
    ## Portion to do in parallel
    fits <- lapply(respGroups, function(rg) {
      bdots_fit(dat[responseGroup == rg, ], curveType[[rg]], concave, cor, rho.0, cores, verbose)
    })
  }

  
  
  return(dat)

}




test <- bdots(cohort_unrelated, 
              group = "DB_cond", 
              subject = "Subject", 
              time = "Time", 
              response = "Fixations", 
              responseGroup = "LookType")

## For testing above
dat <- test


































