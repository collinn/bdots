# pdf for logistic function
# logistic_f()

logisticPars <- function(dat) {
  time <- dat$time
  y <- dat$y

  mini <- min(y)
  peak <- max(y)
  r <- (peak - mini)
  cross <- time[which.min(abs(.5 * r - y))]

  # slope
  q75 <- .75 * r + mini
  q25 <- .25 * r + mini
  time75 <- time[which.min(abs(q75 - y))]
  time25 <- time[which.min(abs(q25 - y))]
  slope <- (q75 - q25) / (time75 - time25)

  return(c(mini = mini, peak = peak, slope = slope, cross = cross))
}

estLogisticCurve <- function(dat, rho, params = NULL,
                              get.cov.only = FALSE, refits = FALSE) {

  if (is.null(params)) {
    params <- logisticPars(dat)
  } else {
    print(length(params))
    print(params)
    if (length(params) != 4) stop("logistic requires 4 parameters be specified for refitting")
    if (!all(names(params) %in% c("mini", "peak", "slope", "cross"))) {
      stop("logistic parameters for refitting must be correctly labeled")
    }
  }

  ff <- quote(y ~ mini + (peak - mini) / (1 + exp(4 * slope * (cross - (time)) / (peak - mini))))

  fit <- curveFitter(dat, ff, params, rho, refits, get.cov.only, ...)
  return(list(fit = fit, ff = ff))
}
