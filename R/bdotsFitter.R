#' Fits Individual Subject Curve
#'
#' The one subject version of bdotsFit
#'
#' @param dat data for single subject/group combo
#' @param curveType this is actually a function. Should rename
#' @param rho correlation coefficient
#' @param numRefits number of refit attempts
#' @param verbose not used
#' @param get.cov.only holdover from old bdots. Follow up with Jake to see if used
#' @param params starting parameters, if wanting to add manually
#' @param splitVars variables used to identify group. Might combine with datVarNames
#' @param datVarNames character vector indicating reponse and time values from parent call
#' @param ... not used
#'
#' @import data.table
bdotsFitter <- function(dat, curveType, rho, numRefits = 0,
                        verbose, get.cov.only = NULL,
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
  ff <- res[['formula']]
  params <- res[['params']]
  fit <- curveFitter(dat, ff, params, rho, numRefits, get.cov.only, ...)


  ## Return DT for failed fit
  if (is.null(fit)) {
    fn <- do.call(c, dat[1, splitVars, with = FALSE])
    dt <- as.data.table(matrix(c(fn, c("fit", "R2", "AR1", "fitCode")),
                               ncol = length(fn) + 4))
    names(dt) <- c(names(fn), c("fit", "R2", "AR1", "fitCode"))
    dt$fit <- I(list(NULL))
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
  fitCode <- 3L*(!hasCor) + 1L*(R2 < 0.95)*(R2 > 0.8) + 2L*(R2 < 0.8)

  ## This is an UNREASONABLY large object. makes up most of size of bdobject
  if (hasCor) {
    attr(fit$modelStruct$corStruct, "factor") <- NULL
  }

  ## Make return DT
  fn <- do.call(c, dat[1, splitVars, with = FALSE])
  dt <- as.data.table(matrix(c(fn, c("fit", "R2", "AR1", "fitCode")),
                             ncol = length(fn) + 4))
  names(dt) <- c(names(fn), c("fit", "R2", "AR1", "fitCode"))
  dt$fit <- I(list(fit))
  dt$R2 <- R2
  dt$AR1 <- (fitCode < 3)
  dt$fitCode <- fitCode

  ## Cleanest way to get formula, for now
  attr(dt, "formula") <- ff
  dt
}


