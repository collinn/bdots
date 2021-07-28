#' DoubleGauss2 curve function for nlme
#'
#' DoubleGauss2 function used in fitting nlme curve for observations
#'
#' @param dat subject data to be used
#' @param y outcome variable, character vector
#' @param time time variable, character vector
#' @param concave Boolean
#' @param params \code{NULL} unless user wants to specify starting parameters for gnls
#' @param ... just in case
#'
#' @details User should only have to worry about setting concavity
#' of this function. Presently only work for time series scaled out to 2000ms
#'
#' \code{y ~ (time < mu) * (exp(-1 * (time - mu) ^ 2
#' / (2 * sig1 ^ 2)) * (ht - base1) + base1)
#' + (mu <= time) * (exp(-1 * (time - mu) ^ 2
#'                          / (2 * sig2 ^ 2)) * (ht - base2) + base2)}
#' @export
doubleGauss2 <- function(dat, y, time, params = NULL, concave = TRUE, ...) {

  if (is.null(params)) {
    params <- dgaussPars2(dat, y, time, concave)
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


dgaussPars2 <- function(dat, y, time, conc = TRUE) {
  time <- dat[[time]]
  y <- dat[[y]]

  ## Remove cases with zero variance
  if (var(y) == 0) {
    return(NULL)
  }

  pars <- setNames(rep(0, 6), c("mu", "ht", "sig1", "sig2", "base1", "base2"))

  ## Based on max height of y
  idx <- ifelse(conc, which.max(y), which.min(y))
  pars['mu'] <- time[idx]
  pars['ht'] <- y[idx]

  ## Set upper and lower bounds
  lb <- c(0, 0, 4*(time[2] - time[1]), 4*(time[2] - time[1]), 0, 0)
  ub <- c(max(time), pars['ht'] * 1.25, max(time)/2, max(time)/2, pars['ht'], pars['ht'])

  ## determine starting base
  n <- length(y)
  if (idx < 0.1*n) {
    pars['base1'] <- mean(y[1:idx])
  } else {
    pars['base1'] <- mean(y[1:idx][1:round(0.1 * n)])
  }

  if (n - idx < 0.1*n) {
    pars['base2'] <- mean(y[(idx + 1):n])
  } else {
    pars['base2'] <- mean(y[(idx + 1):n][1:round(0.1 * n)])
  }

  ## Determine sigmas

  ## Scaling function for parameters (see where these times come from)
  # or when used? dgls in matlab just undoes it
  scalePars <- function(pars, dir = 1) {
    scale_mat <- matrix(c(2000, 1, 500, 500, 1, 1, rep(-0.5, 6)),
                        nrow = 2, byrow = TRUE)
    if (dir == 1) {
      pars <- (pars + scale_mat[2, ]) / scale_mat[1, ]
    } else {
      pars <- pars * scale_mat[1, ] - scale_mat[2, ]
    }
    pars
  }

  ## least squares from doublegauss
  dgls <- function(p, time, y) {
    p <- scalePars(p, dir = 2)
    which_gauss <- time < p['mu']
    y1 <- exp(-1 * (time - p['mu'])^2 / (2 * p['sig1']^2)) * (p['ht'] - p['base1']) + p['base1']
    y2 <- exp(-1 * (time - p['mu'])^2 / (2 * p['sig2']^2)) * (p['ht'] - p['base2']) + p['base2']
    yh <- which_gauss * y1 + (1 - which_gauss) * y2
    mean((y - yh)^2)
  }


  search_size <- 20
  for (i in 3:4) {
    sig_vals <- seq(lb[i], ub[i], length.out = search_size)
    par_attempts <- matrix(rep(pars, search_size),
                           byrow = TRUE, ncol = length(pars))
    par_attempts[, i] <- sig_vals
    ## idx for other sig start vals
    oidx <- 12/i # 3 to 4, 4 to 3
    par_attempts[, oidx] <- (lb[oidx] + ub[oidx]) / 2
    res <- apply(par_attempts, 1, function(p) {
      p <- scalePars(p)
      names(p) <- names(pars)
      dgls(p, time, y)
    })
    pars[i] <- sig_vals[which.min(res)]
  }

  ## Take 10 stabs at computing start values
  res <- matrix(rep(0, 10 * (length(pars) + 1)), nrow = 10)
  optns <- list(maxit = 1000)
  s_lb <- scalePars(lb)
  s_ub <- scalePars(ub)
  for (attempt in 1:10) {
    startp <- pars
    if (attempt == 2) {
      startp[1] <- startp[1] * .85
    } else if (attempt == 3) {
      startp[1] <- startp[1] * 1.15
    } else if (attempt == 4) {
      startp[3:4] <- startp[3:4] * 0.75
    } else if (attempt == 5) {
      startp[3:4] <- startp[3:4] * 1.5
    } else if (attempt != 1) {
      vv <- runif(6, 0.9, 1.1)
      startp <- startp * vv
    }
    s_pars <- scalePars(startp)
    newpars <- optim(s_pars, fn = dgls, time = time, y = y,
                     method = "L-BFGS-B", lower = s_lb, upper = s_ub, control = optns, hessian = TRUE)
    res[attempt, ] <- c(scalePars(newpars$par, 2), newpars$value)
  }
  ## Starting value with par ests
  best_idx <- which.min(res[, length(pars) + 1])
  best_pars <- setNames(res[best_idx, 1:length(pars)], names(pars))
  return(best_pars)

}
