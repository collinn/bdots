## Actual bdots fitting function, not exported to user

## We could almost assume that this takes subjectDat <- dat[subject == w/e, ]

## test environment
#load(file = "~/packages/bdots/data/test_env.RData")
#dat <- tdat[subject == "2", ]



## bdotsFitter takes a dat of a single individual and their curve type
# It will be simple and robust, and any changes at the subject level
# can be done here, in one place

## Actually, this function CAN be exported. Once the person knows which thing
# failed, they can inspect it OR they can take names from that list and just pass
# their own subset portion of the dataset to do it themselves.
# The modularity of this is coming along nicely, I think
for(i in seq_along(newdat)) {
  bdotsFitter(newdat[[i]])
}


## This will ultimately not be passed to end user.
bdotsFitter <- function(dat, curveType, concave, cor, rho, refits = FALSE,
                        verbose, params = NULL, thenames = NULL, get.cov.only = NULL, ...) {



  ## This will need to be handled differently, as dat may not contain
  # y, time, etc
  if(any(dat[, .N, by = .(time)]$N > 1)) {
    warning("Some subjects have multiple observations for unique time. These will be averaged")
    dat[, y := mean(y), by = .(time)]
    dat <- unique(dat, by = c("time"))
  }

  ## Well, determine what curveType, but for now, just  gauss
  # which means I need to make estLogistic
  # also the order of these variables (and in curveFitter and elsewhere) is fucked
  ## fuck me, no, just found out which function we need and then
  # args <- functionArgs
  # do.call(estCurve, args)

  estCurveFit <- switch(curveType,
    doubleGauss = estDgaussCurve(),
    logistic = estLogisticCurve() #,
    #poly = estPolyCurve()
  )

  fit <- do.call(estCurveFit, fargs)

  #fit <- estDgaussCurve(dat, rho = 0.9, get.cov.only = FALSE, cor = TRUE, refits = 0, conc = TRUE)
  #fit <- estLogisticCurve(dat, rho = 0.9, get.cov.only = FALSE, cor = TRUE, refits = 4)




  ## Practice verbose (put at top)
  if (FALSE) message(thenames)

  if (is.null(fit['fit'])) {
    return(list(fit = fit[["fit"]], R2 = NA, fitCode = 6, ff = fit[['ff']]))
  }

  SSE <- sum(resid(fit[['fit']])^2)
  SSY <- sum((dat$y - mean(dat$y))^2)
  R2 <- 1 - SSE/SSY

  hasCor <- !is.null(fit[['fit']]$modelStruct$corStruct)
  fitCode <- 3*(!hasCor) + 1*(R2 < 0.95)*(R2 > 0.8) + 2*(R2 < 0.8)

  list(fit = fit[["fit"]], R2 = R2, fitCode = fitCode, ff = fit[['ff']])
}




