## This replaces effectiveAlphaMvNorm and
## effectiveAlphaMvT, while also allowing rho, k both > 1

## below is what I'm going to call fwerAlpha to differentiate.
## This is what is described on pg 12 of detecting time-specific differences
## to determine actual FWER alpha. The distinction/possible confusion for these
## functions shouldn't matter too much, as they won't be exported for the user
## but rather will both be contained in some variant of find.mod.alpha (for use in p.adjust)

## To test
# n <- 10
# rho <- c(0.9, 0.99)
# k <- qt(1 - (alpha <- 0.05)/(1:3), df = 9)
# method <- "norm"

# rh0 = correlation coefficient
# k = non-critical region (since we take 1 - pmvt(-k, k))
# n = number observations, df = degrees of freedom for "t"
# methods = c("norm", "t") for dist to use
## Returns effective alpha, given number of tests and correlation coefficient
## The critical value that we need to determine adjusted alpha comes from
## fwerAlpha. subsequent iterations through find.mod.alpha uses the effectiveAlpha
# function to verify that that the effective alpha is indeed 0.05 (or whatever)
# this is iterated through, updating the value k, determined to be the root of
# fwerAlpha - alpha, with find.alpha.mod finally returning the correct alpha
# with a call to pt(k that was found, df, etc)

## Adjustment here - as this is only going to be called in `find.mod.alpha`, we
# can make some assumptoins. 1, rho will never change and 2, we can ensure that length(k) == 1
# so let's make this a functional instead.
## no longer takes k as argument, but return functional that does
## benefit here is we only have to c om pute sigma once (since rho never changes)
effectiveAlpha_f <- function(rho, n = 10, df = NULL, method = "norm") {

  sigma <- diag(rep(1, n))
  sigma <- abs(row(sigma) - col(sigma))
  sigma <- rho ^ sigma

  ## Determine method, assign functional for pmvtnorm or pmvt
  if(method == "t") {
    if(is.null(df)) {
      warning("df not supplied for t distribution, using df = n - 1")
      df <- n - 1
    }
    f <- function(k) {
      out <- 1 - pmvt(lower = -k, upper = k,
                      delta = rep(0, n), df = df, corr = sigma)[1]
      out
    }
  } else if (method != "norm") {
    warning("invalid method supplied, using normal approximation")
    method <- "norm"
    f <- function(k) {
      out <- 1 - pmvnorm(lower = -k, upper = k,
                         mean = rep(0, n), corr = sigma)[1]
    }
  } else {
    f <- function(k) {
      out <- 1 - pmvnorm(lower = -k, upper = k,
                         mean = rep(0, n), corr = sigma)[1]
      out
    }
  }
  f
}


## Pg 12 of detecting time-specific differences, FWER alpha
## This is 1 - P(I_t)P(I_t|I_{t-1})^{N-1}

# n <- 10
# rho <- c(0.9, 0.99)
# k <- qt(1 - (alpha <- 0.05)/(1:3), df = 9)

# rh0 = correlation coefficient
# k = non-critical region (since we take 1 - pmvt(-k, k))
# n = number observations
## Returns effective alpha, given number of tests and correlation coefficient
## Also, there is no reason this shouldn't take scalar values
## (the original went about this recursively)
fwerAlpha <- function(rho, k, n = 10) {
  p1 <- 2 * pnorm(k) - 1 # int_{-k}^k f(x) dx, f ~ N(0, 1)
  mean <- c(0, 0)
  sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  p2 <- pmvnorm(lower = c(-k, -k), upper = c(k, k),
                mean = mean, sigma = sigma,
                algorithm=GenzBretz(abseps=.Machine$double.eps))[1]
  p3 <- p2 / p1
  out <- 1 - p3 ^ (n - 1) * p1
}


## Older version
# effectiveAlpha <- function(rho, k, n = 10, df = NULL, method = "norm") {
#
#   pars <- matrix(c(rep(rho, times = length(k)),
#                    rep(k, each = length(rho))), ncol = 2)
#   sigma <- diag(rep(1, n))
#   sigma <- abs(row(sigma) - col(sigma))
#   out <- vector("numeric", length = nrow(pars))
#
#   if(method == "t") {
#     if(is.null(df)) {
#       warning("df not supplied for t distribution, using df = n - 1")
#       df <- n - 1
#     }
#     for(i in seq_len(nrow(pars))) {
#       pp <- pars[i, ]
#       out[i] <- 1 - pmvt(lower = -pp[2], upper = pp[2],
#                          delta = rep(0, n), df = df, corr = pp[1] ^ sigma)[1]
#     }
#   } else {
#     if(method != "norm") warning("invalid method supplied, using normal approximation")
#
#     for(i in seq_len(nrow(pars))) {
#       pp <- pars[i, ]
#       out[i] <- 1 - pmvnorm(lower = -pp[2], upper = pp[2],
#                             mean = rep(0, n), corr = pp[1] ^ sigma)[1]
#     }
#   }
#   out
# }



