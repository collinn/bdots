## This file contains all of the components of finding the modified alpha, including:
# 1. ar1Solver - determines observed autocorrelation
# 2. ar2ScoreDerivative and ar1Loglikelihood - helpers to ar1Solver
# 3. effectiveAlpha_f - functional used in computing effectiveAlpha
# 4. fwerAlpha - determines familywise alpha calculation
# 5. findModifiedAlpha - primary function used to compute modified alpha (uses functions above)

#' Compute AR1 correlation coefficient
#'
#' Computes value for AR1 correlation coefficient for use in \code{p_adjust}
#'
#' @param t A numeric vector of t-statistics
#'
#' @return Estimated AR1 correlation coefficient
#'
#' @examples
#' t <- rt(1e3, df = 1)
#' rho <- ar1Solver(t)
#'
#'
#' @seealso \code{\link[bdots]{p_adjust}}
#' @importFrom stats ar
#' @export
ar1Solver <- function(t) {
  ## Not sure what motivates these values
  n <- length(t)
  v1 <- sum(t[2:n]*t[1:(n-1)])
  v2 <- n - 2 * sum(t^2) + t[n]^2

  phiEst <- polyroot(c(v1, v2, v1, -n))
  phiEst <- Re(phiEst[abs(Im(phiEst)) < 0.0001]) # Discard complex roots
  phiEst <- phiEst[phiEst >= 0]
  phiEst <- phiEst[ar1ScoreDerivative(phiEst, t) < 0] # Discard roots where not concave
  phiEst <- phiEst[phiEst < 1]
  if(length(phiEst) > 1) {
    phiEst <- phiEst[which.max(ar1Loglikelihood(phiEst, t))]
  } else if (length(phiEst) == 0) {
    warning("AR1 solver could not find a proper root. Using Yule-Walker solver.")
    phiEst <- ar(t, order.max = 1, aic = FALSE)$ar
  }
  phiEst
}

## Is this the score? The derivative of the score?
## Need to ask Jake about these ar1 functions
ar1ScoreDerivative <- function(phi, y) {
  phi2 <- phi ^ 2
  n <- length(y)

  ## Not sure where these values come from either
  v1 <- (1 + phi2) / (1 - phi2) ^ 2 * n
  v2 <- (1 + 3 * phi2) / (1 - phi2)^3 * (2 * sum(y ^ 2) - y[n] ^ 2)
  v3 <- 2 * phi * (3 + phi2) / (1 - phi2) ^ 3 * sum(y[2:n] * y[1:(n - 1)])
  v1 - v2 + v3
}

## Ditto above. Clearly it's normal, though
ar1Loglikelihood <- function(phi, y) {
  n <- length(y)
  v1 <- - n / 2 * log(2 * pi * (1 - phi ^ 2))
  v2 <- y[1] ^ 2 + sum((y[2:n] - phi * y[1:(n - 1)]) ^ 2)
  v3 <- 2 * (1 - phi ^ 2)

  v1 - v2 / v3
}


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


#' Find modified alpha
#'
#' find modified alpha
#'
#' @param rho correlation coefficient
#' @param n number of observations
#' @param df degrees of freedom if method == "t"
#' @param alpha starting alpha from which to adjust
#' @param errorAcc acceptable error for alphastar
#' @param gradDiff gradient steps in algorithm
#' @param cores number of cores. Default is zero, or half of what's available
#' @param verbose will probably remove this
#' @param method either "t" or "norm"
#'
#' @import parallel
#' @import stats
findModifiedAlpha <- function(rho, n, df, alpha = 0.05, errorAcc = 0.001,
                              gradDiff = ifelse(cores > 3, 0.5, 0.1), cores = 0,
                              verbose = FALSE, method = "t") {

  ## consequence of pmv* functions
  if (n > 1000) stop("Modified alpha requires that n <= 1000")

  if (cores < 1) {
    cores <- detectCores()/2
  }

  ## Functional that takes argument `k` for bounds
  effectiveAlpha <- effectiveAlpha_f(rho, n, df, method = method)

  if (Sys.info()['sysname'] == "Darwin") {
    cl <- makePSOCKcluster(cores, setup_strategy = "sequential")
  } else {
    cl <- makePSOCKcluster(cores)
  }

  ## Pg 12 of detecting time-specific differences, FWER alpha
  minVal <- qt(1 - alpha / 2, df)
  maxVal <- qt(1 - alpha / 2, df) * 2

  ## This makes sure that effectiveAlpha.normApprox(k) is in (min.val, max.val)
  while (fwerAlpha(rho, maxVal, n) >= alpha) maxVal <- maxVal * 2

  ## Critical value for desired alpha
  k <- uniroot(function(k) fwerAlpha(rho, k, n) - alpha, interval = c(minVal, maxVal))$root

  # alphaStar_vec <- mclapply(c(k, k - gradDiff, k + gradDiff), effectiveAlpha)

  alphaStar_vec <- parLapply(cl, c(k, k + gradDiff), effectiveAlpha)
  alphaStar_vec <- unlist(alphaStar_vec, use.names = FALSE)

  ## What is the minimized error between this and alpha?
  errorMin <- min(abs(alphaStar_vec - alpha))

  ## Only avoid iteration if min gives acceptable error. Otherwise,
  ## we still want to center this around alphaStar at k
  ## Actually, we can improve on this by avoiding weird assignment all together

  ## For now, only iterate with Newton's method + log linearization
  ## That 2 should be a 3 if we do k, k - gd, k + gd
  while(errorMin > errorAcc) {
    gradEst <- log(alphaStar_vec[2] / alphaStar_vec[1]) / gradDiff
    k <- k - log(alphaStar_vec[1] / alpha) / gradEst
    # alphaStar_vec <- mclapply(c(k, k - gradDiff, k + gradDiff), effectiveAlpha)
    alphaStar_vec <- parLapply(cl, c(k, k + gradDiff), effectiveAlpha)
    alphaStar_vec <- unlist(alphaStar_vec, use.names = FALSE)
    errorMin <- min(abs(alphaStar_vec - alpha))
  }
  stopCluster(cl)
  out <- pt(k, df, lower.tail = FALSE) * 2
  out
}



