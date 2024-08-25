## Contains bdotsFitter and curveFitter, used to fit in bdotsFit and bdotsRefit

## -----------------

#' Fits Individual Subject Curve
#'
#' The one subject version of bdotsFit
#'
#' @param dat data for single subject/group combo
#' @param curveType this is actually a function. Should rename
#' @param rho correlation coefficient
#' @param numRefits number of refit attempts
#' @param verbose not used
#' @param getCovOnly only find covariance matrix from starting parameter values
#' @param params starting parameters, if wanting to add manually
#' @param splitVars variables used to identify group. Might combine with datVarNames
#' @param datVarNames character vector indicating reponse and time values from parent call
#' @param ... not used
#'
#' @details Of particular note here is that rho = 0 is proxy for assuming that
#' ar = FALSE
#'
#'
#' @import data.table
bdotsFitter <- function(dat, curveType, rho, numRefits = 0,
                        verbose, getCovOnly = NULL,
                        params = NULL, splitVars = NULL,
                        datVarNames = NULL, ...) {

  ## variables used for subsetting dat
  y <- datVarNames[['y']]
  time <- datVarNames[['time']]

  arggs <- as.list(environment())

  ## what I maybe want to make is named_dots
  v <- list(...)
  m <- dots(...)
  names(v) <- m
  arggs <- c(arggs, v)
  arggs <- compact(arggs)
  ## Don't like name but w/e
  for(nn in intersect(names(formals(curveType)), names(arggs))) {
    formals(curveType)[[nn]] <- arggs[[nn]]
  }
  res <- curveType()

  ## This occurs in situation with var(y) == 0
  if (is.null(res)) {
    return(NULL)
  } else {
    ff <- res[['formula']]
    params <- res[['params']]
    fit <- curveFitter(dat, ff, params, rho, numRefits, getCovOnly, ...)
  }

  ## Return DT for failed fit
  if (is.null(fit)) {
    fn <- do.call(c, dat[1, splitVars, with = FALSE])
    dt <- as.data.table(matrix(c(fn, c("fit", "R2", "AR1", "fitCode")),
                               ncol = length(fn) + 4))
    names(dt) <- c(names(fn), c("fit", "R2", "AR1", "fitCode"))
    dt[['fit']] <- I(vector("list", length = 1L))
    dt$R2 <- NA
    dt$AR1 <- FALSE
    dt$fitCode <- 6L
    attr(dt, "formula") <- ff
    return(dt)
  }

  SSE <- sum(resid(fit)^2)
  SSY <- sum((dat[[y]] - mean(dat[[y]]))^2)
  R2 <- 1 - SSE/SSY


  hasCor <- !is.null(fit$modelStruct$corStruct)

  ## If ar = FALSE, fitCode should only be 0-3, 6
  if (rho > 0) {
    fitCode <- 3L*(!hasCor) + 1L*(R2 < 0.95)*(R2 > 0.8) + 2L*(R2 < 0.8)
  } else {
    fitCode <- 1L*(R2 < 0.95)*(R2 > 0.8) + 2L*(R2 < 0.8)
  }

  ## This is an UNREASONABLY large object. makes up most of size of bdobject
  if (hasCor) {
    attr(fit$modelStruct$corStruct, "factor") <- NULL
  }

  ## Make return DT
  fn <- do.call(c, dat[1, splitVars, with = FALSE])
  dt <- as.data.table(matrix(c(fn, c("fit", "R2", "AR1", "fitCode")),
                             ncol = length(fn) + 4))
  names(dt) <- c(names(fn), c("fit", "R2", "AR1", "fitCode"))
  #dt$fit <- I(list(fit))
  dt$fit <- list(fit)
  dt$R2 <- R2
  dt$AR1 <- (fitCode < 3)
  dt$fitCode <- fitCode

  ## Cleanest way to get formula, for now
  attr(dt, "formula") <- ff
  dt
}


#' Curve Fitter
#'
#' Used in bdotsFit
#'
#' @param dat data used in building curve
#' @param ff formula used in buildilng curve
#' @param params starting parameters
#' @param rho correlation coefficient
#' @param numRefits number of refit attempts
#' @param getCovOnly only find covariance matrix from starting parameter values
#' @param ... don't know that this is used, can maybe get rid of it
#'
#' @import data.table
#' @import nlme
curveFitter <- function(dat, ff, params, rho, numRefits = 0, getCovOnly = NULL, ...) {

  if (!is.null(getCovOnly) && getCovOnly) {
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
          nudgeVal <- runif(length(params), -0.05, 0.05) * (numRefits - attempts)
          n_params <- params * (1 + nudgeVal)
          fit <- tryCatch(gnls(eval(ff), data = dat, start = n_params, correlation = corAR1(rho)),
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
          nudgeVal <- runif(length(params), -0.05, 0.05) * (numRefits - attempts)
          n_params <- params * (1 + nudgeVal)
          fit <- tryCatch(gnls(eval(ff), data = dat, start = n_params), error = function(e) NULL)
        }
      }

      ## As last resort, have potentially bad fit (also meaning no more NULL fitCode)
      if (is.null(fit)) {

        ## This last resort should also be able to specify AR1 or not
        if (rho) {
          fit <- tryCatch(gnls(eval(ff), start = params, data = dat,
                      correlation = corAR1(rho),
                      control = gnlsControl(maxIter = 0, nlsMaxIter = 0,
                                            msMaxIter = 0, returnObject = TRUE)),
                      error = function(e) NULL)
        } else {
          fit <- tryCatch(gnls(eval(ff), start = params, data = dat,
                      control = gnlsControl(maxIter = 0, nlsMaxIter = 0,
                                            msMaxIter = 0, returnObject = TRUE)),
                      error = function(e) NULL)
        }
      }
    }
  }
  ## cor can be determined by existence of fit$modelStruct$corStuct
  # and in this function, it can be implied by rho
  fit
}
