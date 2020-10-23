#' Effective Alpha Functional
#'
#' Functional that returns function for computing effective alpha for given
#' parameters and distribution
#'
#' @param rho Correlation coefficient
#' @param n Number of observations
#' @param df Degrees of freedom if \code{method = "t"}
#' @param method Character string. Determines distribution for adjusted alpha
#' can be either \code{"norm"} for normal distribution or \code{"t"} for t-dist
#'
#' @importFrom mvtnorm pmvt pmvnorm
effectiveAlpha_f <- function(rho, n = 10, df = NULL, method = "norm") {

  sigma <- diag(rep(1, n))
  sigma <- rho ^ abs(row(sigma) - col(sigma))

  ## Determine method, assign functional for pmvtnorm or pmvt
  if(method == "t") {
    if(is.null(df)) {
      warning("df not supplied for t distribution, using df = n - 1")
      df <- n - 1
    }
    f <- function(k) {
      out <- 1 - mvtnorm::pmvt(lower = -k, upper = k,
                      delta = rep(0, n), df = df, corr = sigma)[1]
      out
    }
  } else {
    if (method != "norm") warning("invalid method supplied, using normal approximation")
    f <- function(k) {
      out <- 1 - mvtnorm::pmvnorm(lower = -k, upper = k,
                         mean = rep(0, n), corr = sigma)[1]
    }
  }
  f
}


#' fwerAlpha
#'
#' Family wise alpha calculation
#' @param rho Correlation coefficient
#' @param k Bounds of non-critical region
#' @param n Number of observations
#'
#' @details Returns effective alpha, given number of tests and the correlation
#' coefficient. This isn't explicitly checked, but there is no reason this function
#' should take any non-scalar values. Derivation of this can be found on pg 12
#' of Jake's 'Detecting time-specific differences'. This function performs the
#' expression \deqn{1 - P(I_t)P(I_t \ | \ I_{t-1})^{N-1}}
#'
#'
#' @import mvtnorm
#' @importFrom stats pnorm
fwerAlpha <- function(rho, k, n = 10) {
  p1 <- 2 * pnorm(k) - 1 # int_{-k}^k f(x) dx, f ~ N(0, 1)
  mean <- c(0, 0)
  sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  p2 <- pmvnorm(lower = c(-k, -k), upper = c(k, k),
                mean = mean, sigma = sigma,
                algorithm = GenzBretz(abseps = .Machine$double.eps))[1]
  p3 <- p2 / p1
  out <- 1 - p3 ^ (n - 1) * p1
}



