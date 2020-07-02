library(data.table)
library(magrittr)


## possible additions
# refits - number of times to jitter initial parameters

## for each fit, give numeric value, i.e., 1 - AR1 R2 > 0.95, 2 - AR1 R2 > 0.85, etc.
## this will make refitting SO much easier
## it will make dropping subjects easier too!
## Definitely going to be using data table for this



bdotsFit <- function(data, # dataset
                     subject, # subjects
                     time, # column for time
                     y, # response vector
                     group, # groups for subjects
                     curveType = c("doubleGauss"), # logistic, doubleGauss, etc. maybe have length match responseGroup length?
                     concave = NULL, # doubleGauss only concavity
                     cor = TRUE, # autocorrelation?
                     rho.0 = 0.9, # autocor value
                     cores = 1, # cores to use, 0 == use all available
                     verbose = FALSE) {

  ## data.table parallelizes itself w/ sensible default (50% available)
  setDTthreads(max(getDTthreads(), cores))

  ## Need to validate inputs after things work ##
  ## Also kill all factors ##


  ## Cheat around DT reference, conditionally set key for subset
  ## Can possibly avoid this set, as was done below
  dat <- data.table()
  dat$subject <- data[[subject]] %>% as.character()
  dat$time <- data[[time]] %>% as.numeric()
  dat$y <- data[[y]] %>% as.numeric()

  ## Set group variables in data.table
  for(gg in seq_along(group)) {
    set(dat, j = group[gg], value = data[[group[gg]]])
  }

  f <- function(dat, ct, cc, rho) {
    q1 <- list(res = list(a = dat, b = class(dat), c = cc, curveType = rho))
    return(q1)
  }

  ## Set key
  setkeyv(dat, c("subject", group))

  ## HOLY SHIT THIS WORKS
  ## DOUBLE HOLY SHIT, USING SD PASSES ENTIRE DATA TABLE SUBSET !!
  ## just replace f with bdotsFitter, once finished, goddamn, doing my job for me
  ## .SD going in only contains exactly what is needed to fit function
  rr <- dat[, list(f(.SD, curveType, concave, cor, rho.0, verbose)), by = group,
            .SDcols = c("subject", "time", "y")]



  return(dat)

}

#### LIVE TESTING ####
load("~/bdots/demo.RData")
data <- cohort_unrelated
group = c("Group", "LookType")
time <- "Time"
subject <- "Subject"
y <- "Fixations"
concave <-  TRUE
rho.0 <-  0.9
cor <-  TRUE
cores <-  1
verbose <-  FALSE
curveType <- c("doubleGauss")

## made within function above (at one points)
tdat <- rr$V1[1][[1]][[1]]

test <- bdots(cohort_unrelated,
              subject = "Subject",
              time = "Time",
              group = c("Group", "LookType"),
              y = "Fixations")

## For testing above
dat <- test


































