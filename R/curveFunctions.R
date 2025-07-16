#' Logistic curve function for nlme
#'
#' Logistic function used in fitting nlme curve for observations
#'
#' @param dat subject data to be used
#' @param y outcome variable
#' @param time time variable
#' @param params \code{NULL} unless user wants to specify starting parameters for gnls
#' @param ... just in case
#'
#' @details \code{y ~ mini + (peak - mini) / (1 + exp(4 * slope * (cross - (time)) / (peak - mini)))}
#' @export
logistic <- function(dat, y, time, params = NULL, ...) {

  logisticPars <- function(dat, y, time, ...) {
    time <- dat[[time]]
    y <- dat[[y]]

    # idx <- order(time)
    # time <- time[idx]
    # y <- y[idx]

    ## Remove cases with zero variance
    if (var(y) == 0) {
      return(NULL)
    }

    mini <- min(y)
    peak <- max(y)
    r <- (peak - mini)
    cross <- time[which.min(abs(0.5*r - y))]

    # slope
    q75 <- .75 * r + mini
    q25 <- .25 * r + mini
    time75 <- time[which.min(abs(q75 - y))]
    time25 <- time[which.min(abs(q25 - y))]

    # need to not accidentally get Inf which happens when dividing by zero
    tr <- time75 - time25
    tr <- ifelse(tr == 0, median(time), tr)

    slope <- (q75 - q25) / tr

    return(c(mini = mini, peak = peak, slope = slope, cross = cross))
  }

  if (is.null(params)) {
    params <- logisticPars(dat, y, time)
  } else {
    if (length(params) != 4) stop("logistic requires 4 parameters be specified for refitting")
    if (!all(names(params) %in% c("mini", "peak", "slope", "cross"))) {
      stop("logistic parameters for refitting must be correctly labeled")
    }
  }
  ## Return NA list if var(y) is 0
  if (is.null(params)) {
    return(NULL)
  }
  y <- str2lang(y)
  time <- str2lang(time)
  ff <- bquote(.(y) ~ mini + (peak - mini) / (1 + exp(4 * slope * (cross - (.(time))) / (peak - mini))))
  attr(ff, "parnames") <- names(params)
  return(list(formula = ff, params = params))
}


## -----------

#' Double Gauss curve function for nlme
#'
#' Double Gauss function used in fitting nlme curve for observations
#'
#' @param dat subject data to be used
#' @param y outcome variable, character vector
#' @param time time variable, character vector
#' @param concave Boolean
#' @param params \code{NULL} unless user wants to specify starting parameters for gnls
#' @param ... just in case
#'
#' @details User should only have to worry about setting concavity
#' of this function
#'
#' \code{y ~ (time < mu) * (exp(-1 * (time - mu) ^ 2
#' / (2 * sig1 ^ 2)) * (ht - base1) + base1)
#' + (mu <= time) * (exp(-1 * (time - mu) ^ 2
#'                          / (2 * sig2 ^ 2)) * (ht - base2) + base2)}
#' @export
doubleGauss <- function(dat, y, time, params = NULL, concave = TRUE, ...) {

  if (is.null(params)) {
    params <- dgaussPars(dat, y, time, concave)
  } else {
    if (length(params) != 6) stop("doubleGauss requires 6 parameters be specified for refitting")
    if (!all(names(params) %in% c("mu", "ht", "sig1", "sig2", "base1", "base2"))) {
      stop("doubleGauss parameters for refitting must be correctly labeled")
    }
  }
  ## Return NA list if var(y) is 0
  if (is.null(params)) {
    return(NULL)
  }

  y <- str2lang(y)
  time <- str2lang(time)
  ff <- bquote(.(y) ~ (.(time) < mu) * (exp(-1 * (.(time) - mu) ^ 2
                                            / (2 * sig1 ^ 2)) * (ht - base1) + base1)
               + (mu <= .(time)) * (exp(-1 * (.(time) - mu) ^ 2
                                        / (2 * sig2 ^ 2)) * (ht - base2) + base2))
  attr(ff, "parnames") <- names(params)
  return(list(formula = ff, params = params))
}

## Take this out of doubleGauss function proper
dgaussPars <- function(dat, y, time, conc) {
  time <- dat[[time]]
  y <- dat[[y]]

  ## Remove cases with zero variance
  if (var(y) == 0) {
    return(NULL)
  }

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
  y2 <- y2[time >= mu]
  time2 <- time[time >= mu]
  totalY2 <- sum(y2)
  sigma2 <- time2[which.min(abs((pnorm(1) - pnorm(-1)) * totalY2 - cumsum(y2)))] - mu

  return(c(mu = mu, ht = ht, sig1 = sigma1, sig2 = sigma2,
           base1 = base1, base2 = base2))
}

## -----------

#' Polynomial curve function for nlme
#'
#' Polynomial function used in fitting nlme curve for observations
#'
#' @param dat subject data to be used
#' @param y outcome variable
#' @param time time variable
#' @param degree degree of polynomial
#' @param raw Boolean, use raw polynomials?
#' @param params \code{NULL} unless user wants to specify starting parameters for gnls
#' @param ... just in case
#'
#' @details It's recommended that one uses raw polynomials for this function for
#' numerical stability. As inference is not performed on the parameters themselves,
#' this should have minimial consequences
#'
#' @details \code{y ~ mini + (peak - mini) / (1 + exp(4 * slope * (cross - (time)) / (peak - mini)))}
#' @export
polynomial <- function(dat, y, time, degree, raw = TRUE, params = NULL, ...) {

  polyPars <- function(dat, y, time, degree, raw) {
    ## Remove cases with zero variance
    if (var(dat[[y]]) == 0) {
      return(NULL)
    }
    pp <- lm(dat[[y]] ~ poly(dat[[time]], degree = degree, raw = raw))
    setNames(coef(pp), c(paste0("beta", seq(degree + 1L))))
  }

  if (is.null(params)) {
    params <- polyPars(dat, y, time, degree, raw)
  } else {
    if (length(params) != degree + 1L) stop("poly requires that number of params matches degree + 1 (for intercept)")
    if (!all(names(params) %in% paste0("beta", seq(degree + 1L)))) {
      stop("polynomial parameters must be labeled beta1, ..., beta[degree + 1]")
    }
  }
  ## Return NA list if var(y) is 0
  if (is.null(params)) {
    return(NULL)
  }

  time_names <- paste0("(Time^", seq(degree + 1L) - 1L , ")")

  ff <- paste(names(params), time_names, sep = "*", collapse = "+")
  ff <- str2lang(ff)
  y <- str2lang(y)
  ff <- bquote(.(y) ~ .(ff))
  attr(ff, "parnames") <- names(params)
  return(list(formula = ff, params = params))
}


#' Linear curve function
#'
#' Linear function used in fitting nlme curve for observations
#'
#' @param dat subject data to be used
#' @param y outcome variable
#' @param time time variable
#' @param params \code{NULL} unless user wants to specify starting parameters for gnls
#' @param ... just in case
#'
#' @details Don't use this function please
#'
#' @details \code{y ~ slope*time + intercept}
#' @export
linear <- function(dat, y, time, params = NULL, ...) {
  linearPars <- function(dat, y, time, ...) {
    time <- dat[[time]]
    y <- dat[[y]]

    ## Remove cases with zero variance
    if (var(y) == 0) {
      return(NULL)
    }
    mm <- (max(y) - min(y)) / max(time)
    bb <- mean(y) - mm * mean(time)
    return(c(intercept = bb, slope = mm))
  }

  if (is.null(params)) {
    params <- linearPars(dat, y, time)
  } else {
    if (length(params) != 2) stop("linear requires 2 parameters be specified for refitting")
    if (!all(names(params) %in% c("intercept", "slope"))) {
      stop("linear parameters for refitting must be correctly labeled")
    }
  }
  ## Return NA list if var(y) is 0
  if (is.null(params)) {
    return(NULL)
  }
  y <- str2lang(y)
  time <- str2lang(time)
  ff <- bquote(.(y) ~ slope * .(time) + intercept)
  attr(ff, "parnames") <- names(params)
  return(list(formula = ff, params = params))
}


#' Exponential curve function
#'
#' Exponential function used in fitting nlme curve for observations
#'
#' @param dat subject data to be used
#' @param y outcome variable
#' @param time time variable
#' @param params \code{NULL} unless user wants to specify starting parameters for gnls
#' @param ... just in case
#'
#' @details Remove any values of zero, or jitter, before using with bdotsFit
#'
#' @details \code{y ~ x_0 exp(k beta)}
#' @export
expCurve <- function(dat, y, time, params = NULL, ...) {
  estExpPars <- function(dat, y, time) {
    tt <- lm(log(dat[[y]]) ~ dat[[time]])
    x0 <- exp(coef(tt)[1])
    k <- coef(tt)[2]
    names(x0) <- names(k) <- NULL
    return(c(x0 = x0, k = k))
  }

  if (any(dat[[y]] <= 0)) {
    message("Subjects with values <= 0 will not be fit")
    return(NULL)
  }

  if (is.null(params)) {
    params <- estExpPars(dat, y, time)
  } else {
    # put checks here
  }
  y <- str2lang(y)
  time <- str2lang(time)
  ff <- bquote(.(y) ~ x0 * exp(.(time) * k))
  attr(ff, "parnames") <- names(params)
  return(list(formula = ff, params = params))
}


#' Logistic saccade curve function for nlme
#'
#' And even cooler logistic fitting function, will think of a better name later
#'
#' @param dat subject data to be used
#' @param y outcome variable
#' @param time time variable
#' @param params \code{NULL} unless user wants to specify starting parameters for gnls
#' @param startSamp how many samples from distribution should we use to investigate
#' @param ... just in case
#'
#' @details \code{y ~ mini + (peak - mini) / (1 + exp(4 * slope * (cross - (time)) / (peak - mini)))}
#' @export
logistic_sac <- function(dat, y, time, params = NULL, startSamp = 8,...) {

  logisticPars <- function(dat, y, time, startSamp, ...) {
    time <- dat[[time]]
    y <- dat[[y]]

    ## Remove cases with zero variance
    if (var(y) == 0) {
      return(NULL)
    }

    ## Sensible distribution of starting pars
    spars <- structure(list(fn = c(1L, 1L, 1L, 1L),
                            param = c("mini", "peak", "slope", "cross"),
                            mean = c(0.115, 0.885, 0.0016, 765), sd = c(0.12, 0.12, 0.00075, 85), min = c(0, 0.5, 0.0009, 300),
                            max = c(0.3, 1, 0.01, 1100)), row.names = c(NA, -4L), class = c("data.table", "data.frame"))

    ## function def
    fn <- function(p, t) {
      b0 <- p[1] # base
      b1 <- p[2] # max
      sl <- p[3] # slope
      xo <- p[4] # crossover
      b0 + (b1-b0) / (1 + exp(4*sl*((xo-t)/(b1-b0))))
    }

    tryPars <- vector("list", length = startSamp)

    ## Get starting pars
    for (i in seq_len(startSamp)) {
      maxFix <- 2
      while (maxFix > 1 | maxFix < 0.6) { # added minimum independent of mbob
        tryPars[[i]] <- Inf
        while (any(spars[, tryPars[[i]] <= min | tryPars[[i]] >= max ])) {
          tryPars[[i]] <- spars[, rnorm(length(tryPars[[i]]))*sd + mean]
        }
        maxFix <- max(fn(tryPars[[i]], time))
      }
    }

    ## Find which is better fit
    r2 <- vector("numeric", length = startSamp)
    for (i in seq_len(startSamp)) {
      yhat <- fn(tryPars[[i]], time)
      r2[i] <- mean((y-yhat)^2)
    }

    finalPars <- tryPars[[which.min(r2)]]
    names(finalPars) <- c("mini", "peak", "slope", "cross")

    return(finalPars)
  }

  if (is.null(params)) {
    params <- logisticPars(dat, y, time, startSamp)
  } else {
    if (length(params) != 4) stop("logistic requires 4 parameters be specified for refitting")
    if (!all(names(params) %in% c("mini", "peak", "slope", "cross"))) {
      stop("logistic parameters for refitting must be correctly labeled")
    }
  }
  ## Return NA list if var(y) is 0
  if (is.null(params)) {
    return(NULL)
  }
  y <- str2lang(y)
  time <- str2lang(time)
  ff <- bquote(.(y) ~ mini + (peak - mini) / (1 + exp(4 * slope * (cross - (.(time))) / (peak - mini))))
  attr(ff, "parnames") <- names(params)
  return(list(formula = ff, params = params))
}





#' Double Gauss curve function for nlme
#'
#' Double Gauss function used in fitting nlme curve for observations but now even
#' better since it samples across a sensible distribution for starting parameters
#'
#' @param dat subject data to be used
#' @param y outcome variable, character vector
#' @param time time variable, character vector
#' @param params \code{NULL} unless user wants to specify starting parameters for gnls
#' @param startSamp how many samples from distribution should we use to investigate
#' @param ... just in case
#'
#' @details User should only have to worry about setting concavity
#' of this function
#'
#' \code{y ~ (time < mu) * (exp(-1 * (time - mu) ^ 2
#' / (2 * sig1 ^ 2)) * (ht - base1) + base1)
#' + (mu <= time) * (exp(-1 * (time - mu) ^ 2
#'                          / (2 * sig2 ^ 2)) * (ht - base2) + base2)}
#' @export
doubleGauss_sac <- function(dat, y, time, params = NULL, startSamp = 8, ...) {

  dgaussPars <- function(dat, y, time, startSamp = 8) {
    time <- dat[[time]]
    y <- dat[[y]]

    ## Remove cases with zero variance
    if (var(y) == 0) {
      return(NULL)
    }

    spars <- structure(list(fn = c(2L, 2L, 2L, 2L, 2L, 2L),
                            param = c("mu",  "ht", "sig1", "sig2", "base1", "base2"),
                            mean = c(630, 0.18, 130, 250, 0.05, 0.05),
                            sd = c(77, 0.05, 30, 120, 0.015, 0.015),
                            min = c(300, 0.05, 50, 50, 0, 0),
                            max = c(1300, 0.35, 250, 400, 0.15, 0.15)),
                       row.names = c(NA, -6L), class = c("data.table","data.frame"))

    ## function def
    fn <- function(p, t) {
      mu <- p[1]
      ht <- p[2]
      s1 <- p[3]
      s2 <- p[4]
      b1 <- p[5]
      b2 <- p[6]
      lhs <- (t < mu) * ((ht-b1) * exp((t - mu)^2/(-2*s1^2)) + b1)
      rhs <- (t >= mu) * ((ht-b2) * exp((t - mu)^2/(-2*s2^2)) + b2)
      lhs+rhs
    }

    tryPars <- vector("list", length = startSamp)

    for (i in seq_len(startSamp)) {
      maxFix <- 2
      while (maxFix > 1) {
        tryPars[[i]] <- Inf
        while (any(spars[, tryPars[[i]] <= min | tryPars[[i]] >= max ])) {
          tryPars[[i]] <- spars[, rnorm(length(tryPars[[i]]))*sd + mean]
        }
        maxFix <- max(fn(tryPars[[i]], time))
      }
    }

    r2 <- vector("numeric", length = startSamp)
    for (i in seq_len(startSamp)) {
      yhat <- fn(tryPars[[i]], time)
      r2[i] <- mean((y-yhat)^2)
    }

    finalPars <- tryPars[[which.min(r2)]]
    names(finalPars) <- c("mu", "ht", "sig1", "sig2", "base1", "base2")

    return(finalPars)
  }

  if (is.null(params)) {
    params <- dgaussPars(dat, y, time, startSamp)
  } else {
    if (length(params) != 6) stop("doubleGauss requires 6 parameters be specified for refitting")
    if (!all(names(params) %in% c("mu", "ht", "sig1", "sig2", "base1", "base2"))) {
      stop("doubleGauss parameters for refitting must be correctly labeled")
    }
  }
  ## Return NA list if var(y) is 0
  if (is.null(params)) {
    return(NULL)
  }

  y <- str2lang(y)
  time <- str2lang(time)
  ff <- bquote(.(y) ~ (.(time) < mu) * (exp(-1 * (.(time) - mu) ^ 2
                                            / (2 * sig1 ^ 2)) * (ht - base1) + base1)
               + (mu <= .(time)) * (exp(-1 * (.(time) - mu) ^ 2
                                        / (2 * sig2 ^ 2)) * (ht - base2) + base2))
  attr(ff, "parnames") <- names(params)
  return(list(formula = ff, params = params))
}
