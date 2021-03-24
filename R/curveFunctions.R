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

    mini <- min(y)
    peak <- max(y)
    r <- (peak - mini)
    cross <- time[which.min(abs(0.5*r - y))]

    # slope
    q75 <- .75 * r + mini
    q25 <- .25 * r + mini
    time75 <- time[which.min(abs(q75 - y))]
    time25 <- time[which.min(abs(q25 - y))]
    slope <- (q75 - q25) / (time75 - time25)

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
#' Logistic function used in fitting nlme curve for observations
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

  time_names <- paste0("I(Time^", seq(degree + 1L) - 1L , ")")

  ff <- paste(names(params), time_names, sep = "*", collapse = "+")
  ff <- str2lang(ff)
  y <- str2lang(y)
  ff <- bquote(.(y) ~ .(ff))
  attr(ff, "parnames") <- names(params)
  return(list(formula = ff, params = params))
}
