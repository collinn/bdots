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

  cl <- makePSOCKcluster(cores)

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














