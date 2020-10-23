#' Plot a bdotsFit object
#'
#' Plot individual fits or model fit parameters from an object of class
#' 'bdotsObj'. These functions are not very stable
#'
#' @param bdObj An object of class 'bdotsObj' returned from \code{bdotsFit}
#' @param fitCode Currently not used
#' @param gridSize Length one numeric indicating size of plot grid. Default is
#' 2x2. For right now, they are square
#' @param plotfun Plot either subject fits or model parameters with "fits" or "pars"
#' @param ... ignore for now (other args to plot.generic)
#'
#' @details Right now, these functions are very unstable and expected to change.
#' The largest current issue is with the placement of the legend, which cannot
#' be adjusted. If you are running into issues with seeing things correctly, try
#' making the "Plots" window in RStudio larger before running this function
#' @export
plot.bdotsObj <- function(bdObj, fitCode, gridSize = NULL, plotfun = "fits", ...) {
  ## Top part of this file are methods for plot.bdotsObj
  ## Bottom part of file are methods for plot.bdotsBootObj

  ## Generic currently works. Hurray
  ## Need to do some sort of match.call here for plotfun
  # option to print to file?
  if (plotfun == 'fits') {
    plotFits(bdObj, fitCode, gridSize, ...)
  } else if (plotfun == 'pars') {
    plotPars(bdObj, fitCode, gridSize, ...)
  } else {
    stop("Invalid plotfun type. See ?plot.bdotsObj")
  }
}

# plot <- function(x, ...) {
#   UseMethod("plot")
# }

# so  plot(class, args) is how it will all look in the end
# add option for them to set own grid
## Might be interesting to do by Group?

#' @importFrom graphics par hist
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

#' @importFrom graphics par lines legend
plotFits <- function(bdObj, fitCode, gridSize = NULL, ...) {

  ## This plots fitted vs obs curves
  # include other things to subset here
  # fitcode, group, etc or arg to add data'
  # or to set par

  ## do subsetting here

  bdCall <- attr(bdObj, "call")
  splitVars <- c(bdCall$subject, eval(bdCall$group))
  y <- bdCall$y # what observed values are in X

  if (!is.character(y)) stop("Error 123")

  X <- attr(bdObj, "X")$X
  dfname <- deparse1(bdCall$data)
  if (is.null(X) & exists(dfname)) {
    X <- get(dfname)
  } else if (is.null(X) & !exists(dfname)) {
    stop("Cannot find dataset used to construct bdObj, please pass as argument")
  }


  ## Assume that subsets have been made
  oldPar <- par()$mfrow
  on.exit(par(mfrow = oldPar))

  if (is.null(gridSize)) {
    gridSize <- 2
  } else if (gridSize == "refit") {
    par(mfrow = c(1, 2))
  } else {
    par(mfrow = c(gridSize, gridSize))
  }

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
    code <- bdObj[i, ]$fitCode
    r2 <- round(as.numeric(bdObj[i, ]$R2), 3)
    if (code == 6) next
    obs <- unlist(bdObj[i, splitVars, with = FALSE])
    obs2 <- paste(obs, collapse = ".")
    obsY <- Xs[[obs2]][[y]]
    fitY <- fitted.values(bdObj[i, ]$fit[[1]])

    ## Janky fix for update
    if (gridSize == "refit") {
      if (i == 1) {
        title <- paste("Original Fit", "\n fitCode = ", code, ", R2 = ", r2)
      }
      if (i == 2) {
        title <- paste("Updated Fit", "\n fitCode = ", code, ", R2 = ", r2)
      }
    } else {
      title <- paste(paste0(obs, collapse = " "), "\n fitCode = ", code, ", R2 = ", r2)
    }


    plot(x = time, y = obsY, lty = 2, lwd = 2, type = 'l',
         ylab = y, main = title, col = 'blue')
    lines(x = time, y = fitY, lty = 1, lwd =2, type = 'l')
    # perhaps change legend based on  plot type?
    legend(lgn, legend = c("Observed", "Model Fit"), lty = c(2, 1),
           lwd = c(2, 2), col = c('blue', 'black'))


  }

}


#' @importFrom graphics matplot
plotDiff <- function(bdBootObj, alpha = 0.05, ...) {
  diff <- bdBootObj[['curveList']][['diff']]
  mm <- makePlotCI(diff, alpha)
  ## add some shaading to this, that would be lit
  matplot(mm, lty = c(2, 1, 2), type = 'l',
          col = c('grey', 'black', 'gray'))
}

#' Plot for object of class bdotsBootObj
#'
#' Allows a number of different but also unstable option for plotting an object
#' of class bdotsBoot
#'
#' @param bdBootObj An object of class bdotsBootObj
#' @param alpha Significance level for plotting shaded regions
#' @param diffs Boolean to plot "difference of difference" curve
#' @param group Specify group to plot if difference of difference was used
#' @param ... ignore for now
#'
#'
#' @details Use with care
#' @export
plot.bdotsBootObj <- function(bdBootObj, alpha = 0.05, diffs = NULL, group = NULL, ...) {
  # Fuck yeah, this is cool
  # diffs indicate if we should also plot the difference in addition to curves
  # group indicates which group, in the case we do diffs of diff
  # should also have an argument to alpha adjust for diff of diff. Because we may
  # end up with LI.M, LI.W and LI.diff, but we don't have significant regions for them.
  # it should definitely be an argument since that computation is expensive.
  # should perhaps consider if there is a way to add new_alphastar to a specific
  # group instead of using memoization. Mmm. But can't do that inside of a function
  # could possibly do it if they do bdBootObj <- plot(bdBootObj)
  ## This is way too flowery. Just make hard conditions and treat as separte functions for now
  plotCompare(bdBootObj, alpha, diffs, group, ...)
}

#' @importFrom graphics par legend
plotCompare <- function(bdBootObj, alpha = 0.05, diffs = NULL, group = NULL, ...) {

  cl <- bdBootObj[['curveList']]
  dod <- bdBootObj[['dod']] # diff of diff
  cvGroups <- bdBootObj[['diffs']]

  if (dod & is.null(diffs)) {
    diffs <- TRUE
  }

  # diff of diff and select group (this is kosher)
  if (dod & !is.null(group)) {
    gpidx <- grep(group, names(cl))
    if (length(gpidx) == 0) stop("Invalid 'group' name") # may check if they accidentally passed innerDiff var
    cl <- cl[[gpidx]]
  }

  ## if not diff of diff and group selected, there can be no diffs
  if (!dod & !is.null(group) & !is.null(diffs)) {
    if (diffs == TRUE) warning("Single curve selected, diff curve not available")
    diffs <- FALSE
  }


  diffidx <- grep("diff", names(cl))

  ## does diff exist? if so, let's break it up
  # if diffs is null, default is TRUE
  if (length(diffidx) != 0) {
    diffCurve <- cl[[diffidx]]
    idx <- seq_along(names(cl))[-diffidx]
    cl <- cl[idx]
    if (is.null(diffs)) diffs <- TRUE
  } else {
    diffs <- FALSE
  }

  if (dod & is.null(group)) {
    cl <- lapply(cl, function(x) {
      x[3]
    })
    names(cl) <- NULL
    cl <- unlist(cl, recursive = FALSE)
  }

  ## Plot matrices
  cvMat <- lapply(cl, makePlotCI, alpha)
  nn <- length(cvMat)
  cvMat <- Reduce(cbind, cvMat)
  cvGroups <- bdBootObj[['diffs']]

  time <- attr(bdBootObj, "bdObjAttr")[["time"]]

  if (nn == 1) {
    plotcol <- 'blue'
  } else if (nn == 2) {
    plotcol <- c("blue", "red")
  }

  if (diffs) {
    oldPar <- par()$mfrow
    on.exit(par(mfrow = oldPar))
    par(mfrow = c(1, 2))
    diffMat <- makePlotCI(diffCurve, alpha)
  }

  matplot(x = time, cvMat, lty = rep(c(2, 1, 2), nn),
          type = 'l', col = rep(plotcol, each = 3))
  # set legend placement (?)
  lgn <- switch(attr(bdBootObj, "bdObjAttr")[['curveType']],
                "logistic" = "bottomright",
                "doubleGauss" = "topright",
                "topright")
  legend(lgn, legend = names(cl), lty = c(2, 1),
         lwd = c(2, 2), col = plotcol)

  if (is.null(group)) {
    sigTime <- bdBootObj[["sigTime"]]
    bucketPlot(sigTime, ylim = c(min(cvMat), max(cvMat)))
  }

  ## Need to compute sigTime if they select for computation
  if (diffs & is.null(group)) {
    matplot(x = time, diffMat, lty = c(2, 1, 2), type = 'l',
            col = c('grey', 'black', 'gray'))
    bucketPlot(sigTime, ylim = c(min(diffMat), max(diffMat)))
  } else if (diffs & !is.null(group)) {
    matplot(x = time, diffMat, lty = c(2, 1, 2), type = 'l',
            col = c('grey', 'black', 'gray'))
  }
}

#' @importFrom grDevices rgb
#' @importFrom graphics rect lines
bucketPlot <- function(sigTime, ylim = c(0, 0.9), ...) {
  ## Add option to change colors later
  if(!is.null(sigTime)) {
    yellow <- rgb(255, 255, 0, alpha = 70, maxColorValue = 255)
    gray <- "gray44"
    for(i in 1:nrow(sigTime)) {
      rect(sigTime[i,1], ylim[1], sigTime[i,2], ylim[2], col = yellow, border = NA)
      lines(c(sigTime[i,1], sigTime[i,1]), ylim, col = gray, lwd = 1)
      lines(c(sigTime[i,2], sigTime[i,2]), ylim, col = gray, lwd = 1)
    }
  }
}


makePlotCI <- function(cl, alpha = 0.05, ...) {
  ## This should take in an object of curveList i.e., 'diff', 'LI.M', etc.
  tv <- qt(1 - alpha / 2, cl[['n']] - 1)
  fit <- cl[['fit']]
  sd <- cl[['sd']]
  mm <- matrix(NA, ncol = 3, nrow = length(fit))
  mm[ ,1] <- fit - sd * tv
  mm[, 2] <- fit
  mm[, 3] <- fit + sd * tv
  mm
}















