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

#bdotsFitter(dat = newdat[[1]], curveList = curveList, rho = 0.9)

## This uses old method, but still very good reference
# bdotsFitter <- function(dat, curveList, rho, numRefits = 0,
#                         verbose, thenames = NULL, get.cov.only = NULL, ...) {
#
#
#   ## Check the dat object to see if params exist
#   ## need to be a bit more careful actually if this is used both
#   # here and bdots refit
#
#   ## This will need to be handled differently, as dat may not contain
#   # y, time, etc
#   if(any(dat[, .N, by = .(time)]$N > 1)) {
#     warning("Some subjects have multiple observations for unique time. These will be averaged")
#     dat[, y := mean(y), by = .(time)]
#     dat <- unique(dat, by = c("time"))
#   }
#
#   ## We can avoid all of this if we have a function passed!
#   ## curve environment
#   ## Obviously will have to change if we change curvetype input
#   curveEnv <- makeCurveEnv(curveList)
#
#   # the names of the function to be called
#   estCurveFit <- switch(names(curveList),
#     doubleGauss = estDgaussCurve,
#     logistic = estLogisticCurve #,
#     #poly = estPolyCurve()
#   )
#
#   ## This is where my error lives and I don't know why
#   ## Named list of current environment
#   ## This has made testing a motherfucker
#   #arggs <- c(as.list(curveEnv), rho = 0.9)
#   #arggs$dat <- dat
#   arggs <- c(as.list(environment()), as.list(curveEnv), list(...))
#
#
#
#   #arggs <- list(dat = dat, curveType = curveType, concave = concave, rho = rho, refits = refits, verbose = FALSE)
#   fit <- do.call(estCurveFit, arggs)
#
#   #fit <- estDgaussCurve(dat, rho = 0.9, get.cov.only = FALSE, cor = TRUE, refits = 0, conc = TRUE)
#   #fit <- estLogisticCurve(dat, rho = 0.9, get.cov.only = FALSE, cor = TRUE, refits = 4)
#
#   ## Practice verbose (put at top)
#   #if (FALSE) message(thenames)
#
#   if (is.null(fit[['fit']])) {
#     return(list(fit = fit[["fit"]], R2 = NA, fitCode = 6, ff = fit[['ff']]))
#   }
#
#   SSE <- sum(resid(fit[['fit']])^2)
#   SSY <- sum((dat$y - mean(dat$y))^2)
#   R2 <- 1 - SSE/SSY
#
#   hasCor <- !is.null(fit[['fit']]$modelStruct$corStruct)
#   fitCode <- 3*(!hasCor) + 1*(R2 < 0.95)*(R2 > 0.8) + 2*(R2 < 0.8)
#
#   list(fit = fit[["fit"]], R2 = R2, fitCode = fitCode, ff = fit[['ff']])
# }

## See if this  works
bdotsFitter <- function(dat, curveType, rho, numRefits = 0,
                        verbose, thenames = NULL, get.cov.only = NULL,
                        params = NULL, ...) {


  if(any(dat[, .N, by = .(time)]$N > 1)) {
    warning("Some subjects have multiple observations for unique time. These will be averaged")
    dat[, y := mean(y), by = .(time)]
    dat <- unique(dat, by = c("time"))
  }
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
  cF <- curveFitter(dat, ff, params, rho, numRefits, get.cov.only, ...)
  fit <- list(fit = cF, ff = ff)


  ## Practice verbose (put at top)
  #if (FALSE) message(thenames)

  if (is.null(fit[['fit']])) {
    return(list(fit = fit[["fit"]], R2 = NA, fitCode = 6, ff = fit[['ff']]))
  }

  SSE <- sum(resid(fit[['fit']])^2)
  SSY <- sum((dat$y - mean(dat$y))^2)
  R2 <- 1 - SSE/SSY

  hasCor <- !is.null(fit[['fit']]$modelStruct$corStruct)
  fitCode <- 3*(!hasCor) + 1*(R2 < 0.95)*(R2 > 0.8) + 2*(R2 < 0.8)

  ## This is an UNREASONABLY large object. makes up most of size of bdobject
  if (hasCor) {
    attr(fit[['fit']]$modelStruct$corStruct, "factor") <- NULL
  }

  list(fit = fit[["fit"]], R2 = R2, fitCode = fitCode, ff = fit[['ff']])
}


