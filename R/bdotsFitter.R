#' Fits Individual Subject Curve
#'
#' The one subject version of bdotsFit
#'
#' @param dat data for single subject/group combo
#' @param curveType this is actually a function. Should rename
#' @param rho correlation coefficient
#' @param numRefits number of refit attempts
#' @param verbose not used
#' @param thenames not used (but otherwise in conjunction with verbose)
#' @param get.cov.only holdover from old bdots. Follow up with Jake to see if used
#' @param params starting parameters, if wanting to add manually
#' @param splitVars variables used to identify group. Might combine with datVarNames
#' @param datVarNames character vector indicating reponse and time values
#' @param ... not used
#'
#' @import data.table
bdotsFitter <- function(dat, curveType, rho, numRefits = 0,
                        verbose, thenames = NULL, get.cov.only = NULL,
                        params = NULL, splitVars = NULL,
                        datVarNames = NULL, ...) {

  ## Actual bdots fitt
  ## bdotsFitter takes a dat of a single individual and their curve type
  # It will be simple and robust, and any changes at the subject level
  # can be done here, in one place

  ## Actually, this function CAN be exported. Once the person knows which thing
  # failed, they can inspect it OR they can take names from that list and just pass
  # their own subset portion of the dataset to do it themselves.
  # The modularity of this is coming along nicely, I think
  # for(i in seq_along(newdat)) {
  #   bdotsFitter(newdat[[i]])
  # }
  #
  #bdotsFitter(dat = newdat[[1]], "doubleGauss", concave, rho = 0.9, verbose = FALSE)


  ## See if this  works
  # splitVars used later to add name to what this returns (in attribute)
  ## And actually, this thing should just return a finished DT
  ## This has been moved into bdotsFit
  # if(any(dat[, .N, by = .(time)]$N > 1)) {
  #   warning("Some subjects have multiple observations for unique time. These will be averaged")
  #   dat[, y := mean(y), by = .(time)]
  #   dat <- unique(dat, by = c("time"))
  # }

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


  ## Practice verbose (put at top)
  #if (FALSE) message(thenames)

  if (is.null(fit)) {
    fn <- do.call(c, dat[1, splitVars, with = FALSE])
    dt <- as.data.table(matrix(c(fn, c("fit", "R2", "AR1", "fitCode")),
                               ncol = length(fn) + 4))
    names(dt) <- c(names(fn), c("fit", "R2", "AR1", "fitCode"))
    dt$fit <- I(list(NULL))
    dt$R2 <- NA
    dt$AR1 <- FALSE
    dt$fitCode <- 6
    attr(dt, "formula") <- ff
    return(dt)
  }

  SSE <- sum(resid(fit)^2)
  SSY <- sum((dat[[y]] - mean(dat[[y]]))^2)
  R2 <- 1 - SSE/SSY

  hasCor <- !is.null(fit$modelStruct$corStruct)
  fitCode <- 3*(!hasCor) + 1*(R2 < 0.95)*(R2 > 0.8) + 2*(R2 < 0.8)

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


