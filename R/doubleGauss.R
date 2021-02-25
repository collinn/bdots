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

  ## Make functions subroutines so that
  # 1) pass less shit into parallel
  # 2) users will need to make inside anyways
  dgaussPars <- function(dat, y, time, conc = concave) {
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

  if (is.null(params)) {
    params <- dgaussPars(dat, y, time, concave)
  } else {
    if (length(params) != 6) stop("doubleGauss requires 6 parameters be specified for refitting")
    if (!all(names(params) %in% c("mu", "ht", "sig1", "sig2", "base1", "base2"))) {
      stop("doubleGauss parameters for refitting must be correctly labeled")
    }
  }

  ## Maybe have functions that just return this formula for
  # ease of update for end user
  y <- str2lang(y)
  time <- str2lang(time)
  ff <- bquote(.(y) ~ (.(time) < mu) * (exp(-1 * (.(time) - mu) ^ 2
                                     / (2 * sig1 ^ 2)) * (ht - base1) + base1)
              + (mu <= .(time)) * (exp(-1 * (.(time) - mu) ^ 2
                                    / (2 * sig2 ^ 2)) * (ht - base2) + base2))
  attr(ff, "parnames") <- names(params)
  return(list(formula = ff, params = params))
  }
