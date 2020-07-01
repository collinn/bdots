## Actual bdots fitting function, not exported to user
curveType <- "doubleGauss"
concave <-  TRUE 
rho.0 <-  0.9 
cor <-  TRUE 
cores <-  1
verbose <-  FALSE
dat <- dat[responseGroup == "Cohort", ]
bdots_fit <- function(dat, curveType, concave, cor, rho.0, cores, verbose) {
  print(curveType)
  print(head(dat))
  
  time <- sort(unique(dat$time))
  
  ## Aggregate duplicated outcomes
  if(any(dat[, .N, by = .(subject, group, time, response)]$N > 1)) {
    warning("Some subjects have multiple observations for unique time. These will be averaged")
    dat[, response := mean(response), by = .(subject, time, group)]
    dat <- unique(dat, by = c("subject", "group", "time"))
  }
}

## Can be placed in helper functions for doubleGauss
dgauss <- function(time, mu, ht, sig1, sig2, base1, base2) {
  (time < mu) * (exp(-1 * (time - mu) ^ 2 / (2 * sig1 ^ 2)) * (ht - base1) +
                   base1) + (mu <= time) * (exp(-1 * (time - mu) ^ 2 / (2 * sig2 ^ 2)) * (ht - base2) +
                                              base2)
}
