
## Top part of this file are methods for plot.bdotsObj
## Bottom part of file are methods for plot.bdotsBootObj

## Generic currently works. Hurray
## Need to do some sort of match.call here for plotfun
# option to print to file?
plot.bdotsObj <- function(bdObj, fitCode, gridSize = NULL, plotfun = "fits", ...) {
  if (plotfun == 'fits') {
    plotFits(bdObj, fitCode, gridSize, ...)
  } else if (plotfun == 'pars') {
    plotPars(bdObj, fitCode, gridSize, ...)
  } else {
    stop("Invalid plotfun type. See ?plot.bdotsObj")
  }
}

#plot(bdObj, plotfun = "pars")
#plot(bdObj, plotfun = "fits")

# plot.bdotsBootObj <- function(...)



# so  plot(class, args) is how it will all look in the end
# add option for them to set own grid
## Might be interesting to do by Group?
plotPars <- function(bdObj, ...) {
  cc <- coef(bdObj)
  cn <- colnames(cc)

  oldPar <- par()$mfrow
  on.exit(par(mfrow = oldPar))

  sz <- ceiling(sqrt(ncol(cc)))
  par(mfrow = c(sz, sz))

  for (i in seq_len(ncol(cc))) {
    hist(cc[, i], xlab = cn[i], main = cn[i])
  }
}

## This plots fitted vs obs curves
# include other things to subset here
# fitcode, group, etc or arg to add data'
# or to set par
plotFits <- function(bdObj, fitCode, gridSize = NULL, ...) {

  ## do subsetting here

  bdCall <- attr(bdObj, "call")
  splitVars <- c(bdCall$subject, eval(bdCall$group))
  y <- bdCall$y # what observed values are in X

  if (!is.character(y)) stop("Error 123")

  X <- attr(bdObj, "X")
  dfname <- deparse1(bdCall$data)
  if (is.null(X) & exists(dfname)) {
    X <- get(dfname)
  } else if (!exists(dfname)) {
    stop("Cannot find dataset used to construct bdObj, please pass as argument")
  }


  ## Assume that subsets have been made
  oldPar <- par()$mfrow
  on.exit(par(mfrow = oldPar))

  if (is.null(gridSize)) {
    gridSize <- min(2, sqrt(nrow(X)))
  }
  par(mfrow = c(gridSize, gridSize))

  Xs <- split(X, by = splitVars)

  ## Not sure if this is what I want yet
  ## Also not sure about legends
  time <- attr(bdObj, "time")

  # set legend placement (?)
  lgn <- switch(attr(bdObj, "curveType"),
         "logistic" = "bottomright",
         "doubleGauss" = "topright",
         "topright")

  # should also make sure that axes are all the same
  for (i in seq_len(nrow(bdObj))) {
    code <- as.integer(bdObj[i, ]$fitCode)
    r2 <- round(as.numeric(bdObj[i, ]$R2), 3)
    if (code == 6) next
    obs <- unlist(bdObj[i, splitVars, with = FALSE])
    obs2 <- paste(obs, collapse = ".")
    obsY <- Xs[[obs2]][[y]]
    fitY <- fitted.values(bdObj[i, ]$fit[[1]])
    title <- paste(paste0(obs, collapse = " "), "\n fitCode = ", code, ", R2 = ", r2)
    plot(x = time, y = obsY, lty = 2, lwd = 2, type = 'l',
         ylab = y, main = title, col = 'blue')
    lines(x = time, y = fitY, lty = 1, lwd =2, type = 'l')
    # perhaps change legend based on  plot type?
    legend(lgn, legend = c("Observed", "Fit"), lty = c(2, 1),
           lwd = c(2, 2), col = c('blue', 'black'))
  }

}

#####################################################

## For now, only focus on plotting the 'diff' portion
# but later add elements to look at inner/outer diff  separately

# plot.bdotsBootObj <- function(...)
curveList <- bdBootObj[['curveList']]
plotDiff <- function(bdBootObj, alpha = 0.05, ...) {
  diff <- bdBootObj[['curveList']][['diff']]
  mm <- makePlotCI(diff, alpha)
  ## add some shaading to this, that would be lit
  matplot(mm, lty = c(2, 1, 2), type = 'l',
          col = c('grey', 'black', 'gray'))
}

# Fuck yeah, this is cool
plotCompare <- function(bdBootObj, alpha = 0.05, ...) {

  oldPar <- par()$mfrow
  on.exit(par(mfrow = oldPar))
  par(mfrow = c(2, 1))

  ## Don't actually have anything fancy here
  cv1 <- bdBootObj[['curveList']][['LI.M']]
  cv2 <- bdBootObj[['curveList']][['LI.W']]

  mm1 <- makePlotCI(cv1, alpha)
  mm2 <- makePlotCI(cv2, alpha)

  matplot(cbind(mm1, mm2), lty = rep(c(2, 1, 2), 2),
          type = 'l', col = rep(c("blue", "red"), each = 3))

  diff <- bdBootObj[['curveList']][['LI.diff']]
  mm <- makePlotCI(diff, alpha)
  ## add some shaading to this, that would be lit
  matplot(mm, lty = c(2, 1, 2), type = 'l',
          col = c('grey', 'black', 'gray'))

}

cl <- diff
## This should take in an object of curveList i.e., 'diff', 'LI.M', etc.
makePlotCI <- function(cl, alpha = 0.05, ...) {
  n <- cl[['n']]
  tv <- qt(1 - alpha / 2, n - 1)
  fit <- cl[['fit']]
  sd <- cl[['sd']]
  mm <- matrix(NA, nrow = 3, ncol = length(fit))
  mm[1, ] <- fit - sd * tv
  mm[2, ] <- fit
  mm[3, ] <- fit + sd * tv
  t(mm)
}
























