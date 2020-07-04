library(data.table)
library(magrittr)

## currently, only actually using
# data.table
# parallel
# nlme
# we could almost remove data.table....
# since we don't have any operations  w/n groups
# and the dat[, list(f...), by = , .SDcols = ] isn't parallelized
# However, it will be nice for subsetting things like
# result[fitCode == xyz, this thing]
# or plot(result[lkajsdfl, alkdsjf])
## Keeping dependencies down is dope. Eff you stringr

## possible additions
# refits - number of times to jitter initial parameters

## for each fit, give numeric value, i.e., 1 - AR1 R2 > 0.95, 2 - AR1 R2 > 0.85, etc.
## this will make refitting SO much easier
## it will make dropping subjects easier too!
## Definitely going to be using data table for this

# pair verbose with message

## If rho null but cor = TRUE, couldn't we guess at rho?

bdotsFit <- function(data, # dataset
                     subject, # subjects
                     time, # column for time
                     y, # response vector
                     group, # groups for subjects
                     curveType = c("doubleGauss"), # logistic, doubleGauss, etc. maybe have length match responseGroup length?
                     concave = NULL, # doubleGauss only concavity
                     cor = TRUE, # autocorrelation?
                     rho = 0.9, # autocor value
                     cores = 1, # cores to use, 0 == use all available
                     verbose = FALSE) {

  ## data.table parallelizes itself w/ sensible default (50% available)
  setDTthreads(max(getDTthreads(), cores))

  ## Need to validate inputs after things work ##
  ## Also kill all factors ##


  ## Possibility to clean up nonsense below and remove
  ## magrittr
  # ## Set group variables in data.table
  # vn <- c("y", "time", "subject", group)
  # on <-c(y, time, subject, group)
  #
  # dat <- data.table(matrix(NA, ncol = length(vn)))
  # setNames(dat) <- vn
  # for(nm in seq_along(vn)) {
  #   set(dat, j = vn[nm], value = data[[on[nm]]])
  # }

  ## Cheat around DT reference, conditionally set key for subset
  ## Can possibly avoid this set, as was done below
  # and get rid of magrittr
  # maybe get rid of data table :(
  dat <- data.table()
  dat$subject <- data[[subject]]
  dat$time <- data[[time]] %>% as.numeric()
  dat$y <- data[[y]] %>% as.numeric()

  ## Set group variables in data.tableb+
  for(gg in seq_along(group)) {
    set(dat, j = group[gg], value = data[[group[gg]]])
  }


  ## Set key
  setkeyv(dat, c("subject", group))

  ## HOLY SHIT THIS WORKS
  ## DOUBLE HOLY SHIT, USING SD PASSES ENTIRE DATA TABLE SUBSET !!
  ## just replace f with bdotsFitter, once finished, goddamn, doing my job for me
  ## .SD going in only contains exactly what is needed to fit function
  ## dude,  fuck. Make by = c(group, "subject"), and do every thing at once
  ## and in super duper parallel
  # (above) <- nope, data.table doesn't dispatch like that :(
  # rr <- dat[, list(bdotsFitter(.SD, curveType, concave, cor, rho, verbose)), keyby = c("subject", group),
            # .SDcols = c("time", "y")]


 ## Unfortunatley, this is faster x3 (and far less sexy xInf)
 cl <- makePSOCKcluster(4)
 clusterExport(cl, c("newdat", "curveType", "concave", "cor", "rho", "jitter", "verbose",
                     "curveFitter", "estDgaussCurve", "dgaussPars", "gnls", "corAR1"), envir = .GlobalEnv)
 clusterEvalQ(cl, library(data.table))
 clusterEvalQ(cl, library(nlme))
 newdat <- split(dat, by = c("subject", group))
 system.time(res <- parLapply(cl, newdat, bdotsFitter)) # this also needs arguments
 stopCluster(cl)

 #tt <- strsplit(names(newdat), "\\.")
 # this gives me same as before (good), but still need to
 # do processing on 'result'
 #x <- names(newdat)[1]
 tt <- lapply(names(newdat), function(x) {
   result <- res[[x]]
   x <- strsplit(x, "\\.")
   dat1 <- as.data.table(matrix(x[[1]], ncol = length(x[[1]])))
   names(dat1) <-  c("subject", group)

   # I don't really like that this has to be like this
   # it's a length 1 list where result is a list of length 3, but w/e for now
   # since assigning class, could also make method `[` or `[[`
   dat1$fit <- list(result)

   # These got named wierd things from function, so use number. ew
   dat1$R2 <- result['R2']
   dat1$AR1 <- (result['fitCode'] < 2)
   dat1$fitCode <- result['fitCode']
   dat1
 })

 ## This alone is a few seconds (2.097 vs 0.021)
 system.time(tt1 <- Reduce(rbind, tt))
 system.time(tt2 <- rbindlist(tt))

 ## What else should be returned? anything in sub level stuff?
 # perhaps a curve type, so list(curvetype, tt1), for example? what else?
 # think of summary functions


  return(dat)

}
system.time(doubleGauss.fit())
#### LIVE TESTING ####
## values in R/testing_environment.R
load(file = "~/packages/bdots/data/test_env.RData")

test <- bdotsFit(cohort_unrelated,
              subject = "Subject",
              time = "Time",
              group = c("Group", "LookType"),
              y = "Fixations")

## For testing above
dat <- test


































