## Replacement for find.mod.alpha
## jury is out still on whether or not this should be exported to user
## since we can just as easily wrap it in p.adjust and return the modified
## alpha value as an attribute, because, really, who needs it directly?

## For Testing ##
# x <- numeric(750)
# x[1] <- rnorm(1)
# for(i in 2:750) x[i] <- rnorm(1, x[i - 1] * .9, 1 - .9 ^ 2)
# rho.est <- ar(x, FALSE, order.max = 1)$ar
# rho <- rho.est; n <- 500; df <- 74; verbose <- FALSE; alpha <- 0.05;
# gradDiff <- .5; errorAcc <- 0.001; cores <- 1; method <- "t"

## Best way to speed this up with parallel - 2 cores, just split them up
# it's not faster to check (k - gD, k, k + gD), just (k, k + gD) suffices.
# However, with 3 cores, Halley's method is on point.
# Though I find the idea of branching by cores somewhat disagreeable
# If not necessarily faster, this code here is SIGNIFICANTLY shorter and easier to understand
# don't lose much running PSOCKCluster with 1 core, so we will simply check 3 or not
## after testing, turns out Halley's method doesn't help all that much more, but it's
# already coded up, and there may be edge cases where it does matter
findModifiedAlpha <- function(rho, n, df, alpha = 0.05, errorAcc = 0.001,
                              gradDiff = ifelse(cores > 3, 0.5, 0.1),
                              cores = 1, method = "t") {

  ## consequence of pmv* functions
  if (n > 1000) stop("Modified alpha requires that n <= 1000")

  ## Functional that takes argument `k` for bounds
  effectiveAlpha <- effectiveAlpha_f(rho, n, df, method = method)

  ## Pg 12 of detecting time-specific differences, FWER alpha
  minVal <- qt(1 - alpha / 2, df); maxVal <- qt(1 - alpha / 2, df) * 2
  while (fwerAlpha(rho, maxVal, n) >= alpha) maxVal <- maxVal * 2
  k <- uniroot(function(k) fwerAlpha(rho, k, n) - alpha, interval = c(minVal, maxVal))$root

  ## Don't lose much with only 1 core, and removes ambiguity
  cl <- makePSOCKcluster(rep("localhost", cores))
  clusterExport(cl, c("effectiveAlpha", "pmvt", "pmvnorm"))

  if (cores > 2) {
    aStar <- parSapply(cl, c(k, k + gradDiff, k - gradDiff), effectiveAlpha)
  } else {
    aStar <- parSapply(cl, c(k, k + gradDiff), effectiveAlpha)
  }
  ## What is the minimized error between this and alpha?
  errorMin <- min(abs(aStar - alpha))

  while(errorMin > errorAcc) {
    if (cores > 2) { # Halley's method
      gradEst <- (log(aStar[2]) - log(aStar[3])) / (2 * gradDiff)
      hessEst <- (log(aStar[2]) - 2 * log(aStar[1]) + log(aStar[3])) / gradDiff ^ 2
      k <- k - (log(aStar[1]) - log(alpha)) / gradEst *
        (1 - (log(aStar[1]) - log(alpha)) * hessEst / (2 * gradEst ^ 2)) ^ (-1)
      aStar <- parSapply(cl, c(k, k + gradDiff, k - gradDiff), effectiveAlpha)
    } else { # Newton's method
      gradEst <- log(aStar[2] / aStar[1]) / gradDiff
      k <- k - log(aStar[1] / alpha) / gradEst
      aStar <- parSapply(cl, c(k, k + gradDiff), effectiveAlpha)
    }
    errorMin <- min(abs(aStar - alpha))
  }

  stopCluster(cl)
  out <- pt(k, df, lower.tail = FALSE) * 2

  out
}
