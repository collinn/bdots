## Curve fitter
# not exported to user
# called from estDgauss, estLogisitc, etc, w/e functions are available
# BC FUCK, THE USER CAN PASS THEIR OWN FORMULA FOR GLNS

## all of these arguments must be filled in
## If this isn't working for them, in bdotsFit, there can be a (...) argument
# where they can copy this function out, change it, then pass it back in (functional shit is so cool)
# Add argument for minimum cutoffs?
curveFitter <- function(dat, rho, cor, get.cov.only = NULL, refits = 0, ff, params) {

  ## Here's the thing - everything below is going to be used in dgauss, logistic,
  # poly, literally whatever else. That is the benefit of establish ff above, because
  # the rest of these are never going to change. Aw hell yeah!
  if(!is.null(get.cov.only) && get.cov.only) {
    fit <- gnls(eval(ff), start = params, data = data.frame(time, y),
                correlation = corAR1(rho),
                control = gnlsControl(maxIter = 0, nlsMaxIter = 0, msMaxIter = 0, returnObject = TRUE))
    cor <- TRUE
  } else {
    if (cor) {
      fit <- tryCatch(gnls(eval(ff), data = dat, start = params, correlation = corAR1(rho)),
                      error = function(e) NULL)

      if (is.null(fit)) {
        attempts <- refits
        while (attempts > 0 & is.null(fit)) {
          attempts <- attempts - 1
          params <- jitter(params)
          fit <- tryCatch(gnls(eval(ff), data = dat, start = params, correlation = corAR1(rho)),
                          error = function(e) NULL)
        }
        if (is.null(fit)) cor <- FALSE
      }
    }

    if (!cor) {
      fit <- tryCatch(gnls(eval(ff), data = dat, start = params), error = function(e) NULL)

      if (is.null(fit)) {
        attempts <- refits
        while (attempts > 0 & is.null(fit)) {
          attempts <- attempts - 1
          params <- jitter(params)
          fit <- tryCatch(gnls(eval(ff), data = dat, start = params), error = function(e) NULL)
        }
      }
    }
  }
  ## cor can be determined by existence of fit$modelStruct$corStuct
  fit
  #list(fit = fit, cor = cor)
}

