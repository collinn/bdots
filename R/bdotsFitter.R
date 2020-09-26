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
bdotsFitter <- function(dat, curveType, rho, numRefits = 0,
                        verbose, thenames = NULL, get.cov.only = NULL,
                        params = NULL, splitVars = NULL,
                        datVarNames = NULL, ...) {


  ## This has been moved into bdotsFit
  # if(any(dat[, .N, by = .(time)]$N > 1)) {
  #   warning("Some subjects have multiple observations for unique time. These will be averaged")
  #   dat[, y := mean(y), by = .(time)]
  #   dat <- unique(dat, by = c("time"))
  # }

  ## variables used for subsetting dat
  y <- datVarNames[['y']]
  time <- datVarNames[['time']]

  #curveType <- quote(curveType)
  #print(class(curveType))
  arggs <- as.list(environment())
  v <- list(...)
  m <- dots(...)
  names(v) <- m
  arggs <- c(arggs, v)
  arggs <- compact(arggs)
  ## Don't like name but w/e
  for(nn in names(arggs)) {
    formals(curveType)[[nn]] <- arggs[[nn]]
  }
  res <- curveType()
  ff <- res[['formula']]
  params <- res[['params']]
  fit <- curveFitter(dat, ff, params, rho, numRefits, get.cov.only, ...)
  #fit <- list(fit = cF, ff = ff)


  ## Practice verbose (put at top)
  #if (FALSE) message(thenames)

  # if (is.null(fit)) {
  #   return(list(fit = fit, R2 = NA, fitCode = 6, ff = ff)
  # }
  #
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

# #
#  debugonce(bdotsFitter)
# tt <- bdotsFitter(newdat[[1]], curveType = doubleGauss,
#                   rho = rho, numRefits = numRefits,
#                   verbose = FALSE, splitVars = splitVars,
#                   datVarNames = datVarNames)
#
# rr <- tt$fit[[1]]
#
#
# ff <- rr$call
# zz <- ff[[2]]


