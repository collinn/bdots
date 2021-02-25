# Reference for passing plot pars
# https://github.com/pbreheny/ncvreg/blob/master/R/plot-ncvreg.R
# lien 15

#' Plot a bdotsFit object
#'
#' Plot individual fits or model fit parameters from an object of class
#' 'bdotsObj'. These functions are not very stable
#'
#' @param x An object of class 'bdotsObj' returned from \code{bdotsFit}
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
#' @import ggplot2
#' @export
plot.bdotsObj <- function(x, fitCode = NULL, gridSize = NULL, plotfun = "fits", ...) {

  ## Generic currently works. Hurray
  ## Need to do some sort of match.call here for plotfun
  # option to print to file?
  if (plotfun == 'fits') {
    plotFits(x, gridSize, ...)
  } else if (plotfun == 'pars') {
    plotPars(x, gridSize, ...)
  } else {
    stop("Invalid plotfun type. See ?plot.bdotsObj")
  }
}

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
plotFits <- function(bdObj, gridSize = NULL, ...) {

  bdCall <- attr(bdObj, "call")
  splitVars <- c(bdCall$subject, eval(bdCall$group))
  y <- bdCall$y

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
    par(mfrow = c(gridSize, gridSize))
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
  # lgn <- switch(attr(bdObj, "curveType"),
  #        "logistic" = "bottomright",
  #        "doubleGauss" = "topright",
  #        "topright")
  lgn <- "topright" # for now

  # should also make sure that axes are all the same
  for (i in seq_len(nrow(bdObj))) {
    code <- bdObj[i, ]$fitCode
    r2 <- round(as.numeric(bdObj[i, ]$R2), 3)
    obs <- unlist(bdObj[i, splitVars, with = FALSE])
    obs2 <- paste(obs, collapse = ".")
    obsY <- Xs[[obs2]][[y]]
    if (code != 6)
      fitY <- fitted.values(bdObj[i, ]$fit[[1]])

    ## Janky fix for update. Should just make a separate for refits
    if (gridSize == "refit") {
      if (i == 1) {
        title <- paste("Original Fit", "\n fitCode = ", code, ", R2 = ", r2)
      } else {
        title <- paste("Updated Fit", "\n fitCode = ", code, ", R2 = ", r2)
      }
    } else {
      title <- paste(paste0(obs, collapse = " "), "\n fitCode = ", code, ", R2 = ", r2)
    }

    plot(x = time, y = obsY, lty = 1, lwd = 2, type = 'l',
         ylab = y, main = title, col = 'blue')

    if (code != 6) {
      lines(x = time, y = fitY, lty = 1, lwd = 2, type = 'l')
      # perhaps change legend based on  plot type?
      legend(lgn, legend = c("Observed", "Model Fit"), lty = c(1, 1),
             lwd = c(2, 2), col = c('blue', 'black'))
    } else {
      legend(lgn, legend = c("Observed"), lty = 1,
             lwd = 2, col = c('blue'))
    }


  }

}

## Needs title
#' @importFrom graphics matplot
plotDiff <- function(bdBootObj, alpha = 0.05, ...) {
  diff <- bdBootObj[['curveList']][['diff']]
  mm <- makePlotCI(diff, alpha)
  ## add some shaading to this, that would be lit
  matplot(mm, lty = c(2, 1, 2), type = 'l',
          col = c('grey', 'black', 'gray'))
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

## this is valid at least for group fits, not 100% sure of diff
makePlotCI <- function(cl, alpha = 0.05, ...) {
  ## This should take in an object of curveList i.e., 'diff', 'LI.M', etc.
  tv <- qt(1 - alpha / 2, cl[['n']] - 1) # this will not be correct for diff, as 'n' will be wrong.
  fit <- cl[['fit']]
  sd <- cl[['sd']]
  mm <- matrix(NA, ncol = 3, nrow = length(fit))
  mm[ ,1] <- fit - sd * tv # need to divide by sqrt(n) ?
  mm[, 2] <- fit
  mm[, 3] <- fit + sd * tv
  mm
}

#' Plot for object of class bdotsBootObj
#'
#' Allows a number of different but also unstable option for plotting an object
#' of class bdotsBoot
#'
#' @param x An object of class bdotsBootObj
#' @param alpha Significance level for plotting confidence intervals. To readjust
#' areas of significance will value different than alpha used in \code{bdotsBoot} is
#' computationally expensive and is currently not an option but will be offered soon.
#' @param plotDiffs Boolean to plot difference curve
#' @param group Specify group to plot if difference of difference was used. The
#' user can also subset the bdotsBootObj prior to plotting
#' @param ciBands Boolean indicating whether or not to include confidence intervals
#' around fitted curves (currently only option is TRUE)
#' @param ... ignore for now, but will eventually allow plot parameters
#'
#'
#' @details Use with care
#' @export
plot.bdotsBootObj <- function(x, alpha = NULL, plotDiffs = TRUE, group = NULL, ciBands = TRUE, ...) {
  # value used in original call
  alpha <- x$alpha

  if (x$dod & !is.null(group)) {
    x <- subset(x, group, alpha)
    plotInnerGroup(x, alpha, plotDiffs, ciBands, ...)
  } else if (x$dod) {
    plotOuterGroup(x, alpha, plotDiffs, ciBands, ...)
  } else {
    plotInnerGroup(x, alpha, plotDiffs, ciBands, ...)
  }
}

#' @importFrom graphics par legend
plotInnerGroup <- function(bdBootObj, alpha = 0.05, plotDiffs = TRUE, ciBands, ...) {
  cList <- bdBootObj[['curveList']]
  diffList <- cList[['diff']]
  cList <- cList[setdiff(names(cList), "diff")]
  cvMat <- lapply(cList, makePlotCI, alpha)
  nn <- length(cvMat)
  cvMat <- Reduce(cbind, cvMat)

  Time <- attr(bdBootObj, "bdObjAttr")[["time"]]

  ## Set pars
  plotcol <- c("steelblue", "tomato")
  bdCall <- attr(bdBootObj, "bdObjAttr")[['call']]
  yylab <- bdCall[['y']]
  xxlab <- bdCall[['time']]
  mmain <- names(bdBootObj$curveGroups)

  ## If plotDiff
  if (plotDiffs) {
    oldMfrow <- par()$mfrow
    on.exit(par(mfrow = oldMfrow))
    par(mfrow = c(1, 2))
    diffMat <- makePlotCI(diffList, alpha)
  }

  ## With pars, makes more sense to use do.call
  matplot(x = Time, cvMat, lty = rep(c(2,1,2), nn),
          type = 'l', col = rep(plotcol, each = 3),
          xlab = xxlab, ylab = yylab, main = mmain)

  # set legend placement (?) ideally xpd out of plot or below
  # lgn <- switch(attr(bdBootObj, "bdObjAttr")[['curveType']],
  #               "logistic" = "topleft",
  #               "doubleGauss" = "topleft",
  #               "topleft")
  lgn <- "topright" # for now

  ## This only holds if diff of diff not used. Need to handle other case
  if (is.null(NULL) & !is.null(bdBootObj[['sigTime']])) {
    sigTime <- bdBootObj[["sigTime"]]
    bucketPlot(sigTime, ylim = c(min(cvMat), max(cvMat)))
  }

  legend(lgn, legend = names(cList), lty = c(1, 1),
         lwd = c(2, 2), col = plotcol, xpd = TRUE)

  if (plotDiffs) {
    matplot(x = Time, diffMat, lty = c(2, 1, 2), type = 'l',
            col = c('grey', 'black', 'gray'),
            main = "difference of curves", xlab = xxlab,
            ylab = bquote(paste(Delta, .(yylab), sep = " ")))
    if (!is.null(bdBootObj[['sigTime']]))
      bucketPlot(sigTime, ylim = c(min(diffMat), max(diffMat)))
  }
}

## For plotting diff of diff case
plotOuterGroup <- function(bdBootObj, alpha = 0.05, plotDiffs = TRUE, ciBands, ...) {
  cList <- bdBootObj$curveList
  diffList <- cList[['diff']]
  cList <- cList[setdiff(names(cList), "diff")]
  cList <- unlist(cList, recursive = FALSE)
  cList <- cList[grep("diff$", names(cList))]
  tt <- bdBootObj
  cList[['diff']] <- diffList
  tt$curveList <- cList
  plotInnerGroup(tt, alpha, plotDiffs, ciBands, ...)
}




