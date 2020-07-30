## Replacement for find.mod.alpha
## jury is out still on whether or not this should be exported to user
## since we can just as easily wrap it in p.adjust and return the modified
## alpha value as an attribute, because, really, who needs it directly?

## For Testing ##
# x <- numeric(500)
# x[1] <- rnorm(1)
# for(i in 2:500) x[i] <- rnorm(1, x[i - 1] * .9, 1 - .9 ^ 2)
# rho.est <- ar(x, FALSE, order.max = 1)$ar
# rho <- rho.est; n <- 500; df <- 49; verbose <- FALSE; alpha <- 0.05;
# grad.diff <- .5; error.acc <- 0.001; cores <- 1; method <- "t"

findModifiedAlpha <- function(rho, n, df, alpha = 0.05, errorAcc = 0.001,
                              gradDiff = ifelse(cores > 3, 0.5, 0.1), cores = 1,
                              verbose = FALSE, method = "t") {

  ## consequence of pmv* functions
  if (n > 1000) stop("Modified alpha requires that n <= 1000")

  ## Functional that takes argument `k` for bounds
  effectiveAlpha <- effectiveAlpha_f(rho, n, df, method = method)

  ## Pg 12 of detecting time-specific differences, FWER alpha
  minVal <- qt(1 - alpha / 2, df); maxVal <- qt(1 - alpha / 2, df) * 2

  ## This makes sure that effectiveAlpha.normApprox(k) is in (min.val, max.val)
  while (fwerAlpha(rho, maxVal, n) >= alpha) maxVal <- maxVal * 2

  ## Critical value for desired alpha
  k <- uniroot(function(k) fwerAlpha(rho, k, n) - alpha, interval = c(minVal, maxVal))$root

  ## This can be made in parallel later
  alphaStar_vec <- vapply(c(k, k - gradDiff, k + gradDiff),
                          effectiveAlpha, numeric(1))

  ## What is the minimized error between this and alpha?
  errorMin <- min(abs(alphaStar_vec - alpha))

  ## Only avoid iteration if min gives acceptable error. Otherwise,
  ## we still want to center this around alphaStar at k
  ## Actually, we can improve on this by avoiding weird assignment all together

  ## For now, only iterate with Newton's method + log linearization
  while(errorMin > errorAcc) {
    gradEst <- log(alphaStar_vec[2] / alphaStar_vec[1]) / gradDiff
    k <- k - log(alphaStar_vec[1] / alpha) / gradDiff
    alphaStar_vec <- vapply(c(k, k - gradDiff, k + gradDiff),
                            effectiveAlpha, numeric(1))
    errorMin <- min(abs(alphaStar_vec - alpha))
  }

  out <- pt(k, df, lower.tail = FALSE) * 2

  out
}














