## doubleGauss related functions

## pdf for doubleGauss function
dgauss <- function(time, mu, ht, sig1, sig2, base1, base2) {
  (time < mu) * (exp(-1 * (time - mu) ^ 2 / (2 * sig1 ^ 2)) * (ht - base1) +
                   base1) + (mu <= time) * (exp(-1 * (time - mu) ^ 2 / (2 * sig2 ^ 2)) * (ht - base2) +
                                              base2)
}

## misc functions for finding parameters
dgaussPars <- function(time, y, conc) {
  mu <- ifelse(conc, time[which.max(y)], time[which.min(y)])
  ht <- ifelse(conc, max(y), min(y))
  base1 <- ifelse(conc, min(y[time < mu]), max(y[time < mu]))
  base2 <- ifelse(conc, min(y[time > mu]), max(y[time > mu]))
  
  ## A little more involved
  y1 <- y - base1
  y1 <- rev(y1[time <= mu])
  time1 <- rev(time[time <= mu])
  totalY1 <- sum(y1)
  sigma1 <- mu - time1[which.min(abs((pnorm(1) - pnorm(-1)) * totalY1 - cumsum(y1)))]
  
  y2 <- y - base2
  y2 <- rev(y2[time >= mu])
  time2 <- rev(time[time >= mu])
  totalY2 <- sum(y2)
  sigma2 <- time2[which.min(abs((pnorm(1) - pnorm(-1)) * totalY2 - cumsum(y2)))] - mu
  
  return(c(mu = mu, ht = ht, sig1 = sigma1, sig2 = sigma2, 
           base1 = base1, base2 = base2))
}

## I guess let's make jitter a numeric argument, since jitter > 0 implies jitter == TRUE
estDgaussCurve <- function(time, y, rho, conc, params = NULL, 
                           cor = TRUE, get.cov.only = FALSE, jitter = FALSE) {
  if (is.null(params)) {
    params <- dgaussPars(time, y, conc)
  } else {
    if (length(params) != 6) stop("doubleGauss requires 6 parameters be specified for refitting")
    if (!all(names(params) %in% c("mu", "ht", "sig1", "sig2", "base1", "base2"))) {
      stop("doubleGauss parameters for refitting must be correctly labeled")
    }
  }
  
  ff <- quote(y ~ (time < mu) * (exp(-1 * (time - mu) ^ 2
                  / (2 * sig1 ^ 2)) * (ht - base1) + base1) 
                  + (mu <= time) * (exp(-1 * (time - mu) ^ 2 
                  / (2 * sig2 ^ 2)) * (ht - base2) + base2))

  if(get.cov.only) {
    fit <- gnls(eval(ff), start = params, data = data.frame(time, y), 
                correlation = corAR1(rho), 
                control = gnlsControl(maxIter = 0, nlsMaxIter = 0, msMaxIter = 0, returnObject = TRUE))
    cor <- TRUE
  } else {
    if (cor) {
      fit <- tryCatch(gnls(eval(ff), start = params, correlation = corAR1(rho)), 
                      error = function(e) NULL)
      
      if (is.null(fit)) {
        attempts <- jitters
        while (attempts > 0 & is.null(fit)) {
          attempts <- attempts - 1
          params <- jitter(params)
          fit <- tryCatch(gnls(eval(ff), start = params, correlation = corAR1(rho)), 
                          error = function(e) NULL)
        }
        if (is.null(fit)) cor <- FALSE
      }
    }
    
    if (!cor) {
      fit <- tryCatch(gnls(eval(ff), start = params), error = function(e) NULL)
      
      if (is.null(fit)) {
        attempts <- jitters
        while (attempts > 0 & is.null(fit)) {
          attempts <- attempts - 1
          params <- jitter(params)
          fit <- tryCatch(gnls(eval(ff), start = params), error = function(e) NULL)
          }
        }
      }
    }
    list(fit = fit, cor = cor)
}
