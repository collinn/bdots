#' Curve Fitter
#'
#' Used in bdotsFit
#'
#' @param dat data used in building curve
#' @param ff formula used in buildilng curve
#' @param params starting parameters
#' @param rho correlation coefficient
#' @param numRefits number of refit attempts
#' @param get.cov.only idk that this is actually used
#' @param ... don't know that this is used, can maybe get rid of it
#'
#' @import data.table
#' @import nlme
curveFitter <- function(dat, ff, params, rho, numRefits = 0, get.cov.only = NULL, ...) {
  ## Curve fitter
  # not exported to user
  # called from estDgauss, estLogisitc, etc, w/e functions are available
  # BC FUCK, THE USER CAN PASS THEIR OWN FORMULA FOR GLNS

  ## all of these arguments must be filled in
  ## If this isn't working for them, in bdotsFit, there can be a (...) argument
  # where they can copy this function out, change it, then pass it back in (functional shit is so cool)
  # Add argument for minimum cutoffs?

  ## If exported to user, need error checking

  if (!is.null(get.cov.only) && get.cov.only) {
    fit <- gnls(eval(ff), start = params, data = dat,
                correlation = corAR1(rho),
                control = gnlsControl(maxIter = 0, nlsMaxIter = 0, msMaxIter = 0, returnObject = TRUE))
  } else {
    if (rho) { # if rho != 0
      fit <- tryCatch(gnls(eval(ff), data = dat, start = params, correlation = corAR1(rho)),
                      error = function(e) NULL)

      if (is.null(fit)) {
        attempts <- numRefits
        while (attempts > 0 & is.null(fit)) {
          attempts <- attempts - 1
          params <- jitter(params)
          fit <- tryCatch(gnls(eval(ff), data = dat, start = params, correlation = corAR1(rho)),
                          error = function(e) NULL)
        }
        if (is.null(fit)) rho <- 0
      }
    }

    if (!rho) {
      fit <- tryCatch(gnls(eval(ff), data = dat, start = params), error = function(e) NULL)

      if (is.null(fit)) {
        attempts <- numRefits
        while (attempts > 0 & is.null(fit)) {
          attempts <- attempts - 1
          params <- jitter(params)
          fit <- tryCatch(gnls(eval(ff), data = dat, start = params), error = function(e) NULL)
        }
      }
    }
  }
  ## cor can be determined by existence of fit$modelStruct$corStuct
  # and in this function, it can be implied by rho
  fit
}
