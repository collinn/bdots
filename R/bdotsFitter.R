## Actual bdots fitting function, not exported to user

## test environment
load(file = "~/packages/bdots/data/test_env.RData")
dat <- tdat

bdotsFitter <- function(dat, curveType, concave, cor, rho.0, cores, verbose) {


  time <- sort(unique(dat$time))

  ## Aggregate duplicated y values
  if(any(dat[, .N, by = .(subject, time)]$N > 1)) {
    warning("Some subjects have multiple observations for unique time. These will be averaged")
    dat[, y := mean(y), by = .(subject, time)]
    dat <- unique(dat, by = c("subject", "time"))
  }

  ## Requires subject as factor, but then gives us a list of subjects here
  ## alternative could be calling something, by subject, like we did above using by = .()
  ## For now, lets just pretend that we have a single subject object
  subjectDat <- dat[subject == "1", ]

  ## For boot strapping, each subject will need vector of coefficients
  ## curveType used for fitting
  ## standard deviation of these parameters
}

## Can be placed in helper functions for doubleGauss
dgauss <- function(time, mu, ht, sig1, sig2, base1, base2) {
  (time < mu) * (exp(-1 * (time - mu) ^ 2 / (2 * sig1 ^ 2)) * (ht - base1) +
                   base1) + (mu <= time) * (exp(-1 * (time - mu) ^ 2 / (2 * sig2 ^ 2)) * (ht - base2) +
                                              base2)
}
