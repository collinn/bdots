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
#' @details Right now, these functions are a bit unstable and expected to change.
#' The largest current issue is with the placement of the legend, which cannot
#' be adjusted. If you are running into issues with seeing things correctly, try
#' making the "Plots" window in RStudio larger before running this function
#'
#' @returns This will return a list of all of the plots rendered.
#'
#' @export
plot.bdotsObj <- function(x, fitCode = NULL, gridSize = NULL, plotfun = "fits", ...) {

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

  ## Custodial
  parameter <- value <- NULL

  cList <- coefList(bdObj)
  cList <- lapply(cList, function(x) {
    xx <- melt(as.data.table(x), measure.vars = colnames(x),
               variable.name = "parameter", value.name = "value")
  })

  plotList <- vector("list", length(cList))
  names(plotList) <- names(cList)
  for (i in seq_along(cList)) {
    plotList[[i]] <- ggplot(cList[[i]], aes(value, group = parameter, color = 'gray')) + geom_histogram(bins = 8) +
      facet_wrap(~parameter, scales = "free_x") +
      scale_color_manual(values = "white") + ggtitle(names(cList)[i]) +
      theme(legend.position = "none")
  }
  return(invisible(plotList))
}

#' @importFrom graphics par lines legend
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
# This whole thing is actually really awful
plotFits <- function(bdObj, gridSize = NULL, ...) {

  bdCall <- attr(bdObj, "call")
  splitVars <- c(bdCall$subject, eval(bdCall$group))
  yname <- bdCall$y
  tname <- bdCall$time

  if (!is.character(yname)) stop("Error 123")

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
  time <- attr(bdObj, "time")

  # should also make sure that axes are all the same
  plotList <- vector("list", length = nrow(bdObj))
  ## The inside of this loop itself could be a function
  for (i in seq_len(nrow(bdObj))) {
    code <- bdObj[i, ]$fitCode
    r2 <- round(as.numeric(bdObj[i, ]$R2), 3)
    obs <- unlist(bdObj[i, splitVars, with = FALSE])
    obs2 <- paste(obs, collapse = ".")
    obsY <- Xs[[obs2]][[yname]]
    if (code != 6) {
      fitY <- fitted.values(bdObj[i, ]$fit[[1]])
      df <- as.data.table(cbind(time, fitY, obsY))
      #df2 <- melt(df, id.vars = "time", measure.vars = c("obsY", "fitY"))
    } else {
      df <- as.data.table(cbind(time, obsY))
    }
    df2 <- melt(df, id.vars = "time")
    df2$lty <- "dashed"
    df2[variable != "obsY", ]$lty <- "solid"

    df2$clr <- "tomato"
    df2[variable != "obsY", ]$clr <- "steelblue"

    df2[, `Curves` := variable]
    df2$Curves <- "Fit"
    df2[variable == "obsY", ]$Curves <- "Observed"

    ## Janky fix for update. Should just make a separate for refits
    if (gridSize == "refit") {
      if (i == 1) {
        title <- paste("Original Fit", "\n fitCode = ", code, ", R2 = ", r2)
      } else {
        title <- paste("Updated Fit", "\n fitCode = ", code, ", R2 = ", r2)
      }
    } else {
      title <- paste0(paste0(obs, collapse = " "),
                      "\nfitCode = ", code, ", R2 = ", r2, collapse = "")
    }

    ## This can't be the way to do it
    y <- NULL; variable <- NULL; value <- NULL; clr <- NULL; lty <- NULL
    Curves <- NULL;

    colors <- c(
      "#56B4E9",
      "#D55E00"
    )

    ## This thing is kinda gross
    plotList[[i]] <- ggplot(df2, aes(time, value, color = Curves, linetype = Curves)) +
      geom_line(size = 1) + theme_bw() + ggtitle(title) + xlab(tname) + ylab(yname) +
      theme(legend.position = "bottom") +
      scale_color_manual(values = colors)
  }

  ## this is a dumb way to do this
  n <- nrow(bdObj)
  gz <- ifelse(gridSize == "refit", 2, 4)
  idxList <- vector("list", length = ceiling(n/gz))
  rr <- seq(1, n, by = gz)
  for (i in seq_along(rr)) {
    if (rr[i] + gz - 1 > n) {
      idxList[[i]] <- rr[i]:n
    } else {
      idxList[[i]] <- rr[i]:(rr[i]+gz-1)
    }
  }

  ## Preliminary refit plot
  if (gridSize == 1) {
    print(plotList[[1]])
  } else {
    nnrow <- ifelse(gz == 2, 1, 2)
    for (i in seq_along(idxList)) {
      do.call(grid.arrange, c(plotList[idxList[[i]]], nrow = nnrow, ncol = 2))
    }
  }
  return(invisible(plotList))
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

#### DELETE after confirming ggplot works
## @importFrom grDevices rgb
## @importFrom graphics rect lines
# bucketPlot <- function(sigTime, ylim = c(0, 0.9), ...) {
#   ## Add option to change colors later
#   if(!is.null(sigTime)) {
#     yellow <- rgb(255, 255, 0, alpha = 70, maxColorValue = 255)
#     gray <- "gray44"
#     for(i in 1:nrow(sigTime)) {
#       rect(sigTime[i,1], ylim[1], sigTime[i,2], ylim[2], col = yellow, border = NA)
#       lines(c(sigTime[i,1], sigTime[i,1]), ylim, col = gray, lwd = 1)
#       lines(c(sigTime[i,2], sigTime[i,2]), ylim, col = gray, lwd = 1)
#     }
#   }
# }

#### DELETE after confirming ggplot works
## this is valid at least for group fits, not 100% sure of diff
# makePlotCI <- function(cl, alpha = 0.05, ...) {
#   ## This should take in an object of curveList i.e., 'diff', 'LI.M', etc.
#   tv <- qt(1 - alpha / 2, cl[['n']] - 1) # this will not be correct for diff, as 'n' will be wrong.
#   fit <- cl[['fit']]
#   sd <- cl[['sd']]
#   mm <- matrix(NA, ncol = 3, nrow = length(fit))
#   mm[ ,1] <- fit - sd * tv # need to divide by sqrt(n) ?
#   mm[, 2] <- fit
#   mm[, 3] <- fit + sd * tv
#   mm
# }

makePlotCI <- function(cl, alpha = 0.05, ...) {
  tv <- stats::qt(1 - alpha / 2, cl[['n']] - 1) # this will not be correct for diff, as 'n' will be wrong.
  fit <- cl[['fit']]
  sd <- cl[['sd']]
  mm <- matrix(NA, ncol = 3, nrow = length(fit))
  mm[ ,1] <- fit - sd * tv # need to divide by sqrt(n) ?
  mm[, 2] <- fit
  mm[, 3] <- fit + sd * tv
  colnames(mm) <- c("seL", "y", "seU")
  data.table::as.data.table(mm)
}


#' Plot for object of class bdotsBootObj
#'
#' Allows a number of different but also unstable option for plotting an object
#' of class bdotsBoot
#'
#' @param x An object of class bdotsBootObj
#' @param alpha Significance level for plotting confidence intervals. To readjust
#' areas of significance will value different than alpha used in \code{bdotsBoot} is
#' computationally expensive and is currently not an option.
#' @param ciBands Boolean indicating whether or not to include confidence intervals
#' around fitted curves
#' @param plotDiffs Boolean to plot difference curve
#' @param group Specify group to plot if difference of difference was used. The
#' user can also subset the bdotsBootObj prior to plotting
#' @param ... ignore for now, but will eventually allow plot parameters
#'
#' @returns List of ggplot objects, which may be helpful if the margins are weird
#'
#' @details This plot function is also a bit unstable and is expected to change
#' @export
plot.bdotsBootObj <- function(x, alpha = NULL, ciBands = TRUE, plotDiffs = TRUE, group = NULL, ...) {
  # value used in original call
  alpha <- x$alpha

  if (x$dod & !is.null(group)) {
    x <- subset(x, group, adjustAlpha = alpha) # this is a method
    plotInnerGroup(x, alpha, ciBands, plotDiffs, ...)
  } else if (x$dod) {
    plotOuterGroup(x, alpha, ciBands, plotDiffs, ...)
  } else {
    plotInnerGroup(x, alpha, ciBands, plotDiffs, ...)
  }
}

#' @importFrom graphics par legend
# plotInnerGroup <- function(bdBootObj, alpha = 0.05, plotDiffs = TRUE, ciBands, ...) {
#   cList <- bdBootObj[['curveList']]
#   diffList <- cList[['diff']]
#   cList <- cList[setdiff(names(cList), "diff")]
#   cvMat <- lapply(cList, makePlotCI, alpha)
#   nn <- length(cvMat)
#   cvMat <- Reduce(cbind, cvMat)
#
#   Time <- attr(bdBootObj, "bdObjAttr")[["time"]]
#
#   ## Set pars
#   plotcol <- c("steelblue", "tomato")
#   bdCall <- attr(bdBootObj, "bdObjAttr")[['call']]
#   yylab <- bdCall[['y']]
#   xxlab <- bdCall[['time']]
#   mmain <- names(bdBootObj$curveGroups)
#
#   ## If plotDiff
#   if (plotDiffs) {
#     oldMfrow <- par()$mfrow
#     on.exit(par(mfrow = oldMfrow))
#     par(mfrow = c(1, 2))
#     diffMat <- makePlotCI(diffList, alpha)
#   }
#
#   ## With pars, makes more sense to use do.call
#   matplot(x = Time, cvMat, lty = rep(c(2,1,2), nn),
#           type = 'l', col = rep(plotcol, each = 3),
#           xlab = xxlab, ylab = yylab, main = mmain)
#
#   # set legend placement (?) ideally xpd out of plot or below
#   lgn <- switch(as.character(attr(bdBootObj, "bdObjAttr")[['curveType']]),
#                 "logistic" = "topleft",
#                 "doubleGauss" = "topleft",
#                 "topleft")
#   #lgn <- "topright" # for now
#
#   ## This only holds if diff of diff not used. Need to handle other case
#   if (is.null(NULL) & !is.null(bdBootObj[['sigTime']])) {
#     sigTime <- bdBootObj[["sigTime"]]
#     bucketPlot(sigTime, ylim = c(min(cvMat), max(cvMat)))
#   }
#
#   legend(lgn, legend = names(cList), lty = c(1, 1),
#          lwd = c(2, 2), col = plotcol, xpd = TRUE)
#
#   if (plotDiffs) {
#     matplot(x = Time, diffMat, lty = c(2, 1, 2), type = 'l',
#             col = c('grey', 'black', 'gray'),
#             main = "difference of curves", xlab = xxlab,
#             ylab = bquote(paste(Delta, .(yylab), sep = " ")))
#     if (!is.null(bdBootObj[['sigTime']]))
#       bucketPlot(sigTime, ylim = c(min(diffMat), max(diffMat)))
#   }
# }

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
  plotInnerGroup(tt, alpha, ciBands, plotDiffs, dod = TRUE, ...)
}


## dod for plot titles
plotInnerGroup <- function(bdBootObj, alpha = 0.05, ciBands = TRUE, plotDiffs = TRUE, dod = FALSE, ...) {

  ## Custodial tasks
  Groups <- V1 <- V2 <- parameter <- seL <- seU <- value <- y <- NULL

  ggplot2::theme_set(theme_bw())
  cList <- bdBootObj[['curveList']]
  diffList <- cList[['diff']]
  cList <- cList[setdiff(names(cList), "diff")]

  cvMat <- lapply(cList, makePlotCI, alpha)
  cvMat <- data.table::rbindlist(cvMat, idcol = 'Groups')

  ## Likely no chance this will affect a real group name, but need for legend
  cvMat[, Groups := gsub("\\.diff$", "", Groups)]

  Time <- attr(bdBootObj, "bdObjAttr")[["time"]]
  bdCall <- attr(bdBootObj, "bdObjAttr")[['call']]
  yylab <- bdCall[['y']]
  xxlab <- bdCall[['time']]
  if (!dod) {
    mmain <- "Bootstrapped Fits"
  } else {
    mmain <- "Bootstrapped Differences"
  }

  cvMat$Time <- rep(Time, times = length(cList))

  ## sigTime
  sigT <- data.table(bdBootObj$sigTime)

  p <- ggplot(cvMat, aes(Time, y)) +
    geom_line(aes(group = Groups, color = Groups), size = 1) +
    labs(x = xxlab, y = yylab, title = mmain) +
    theme(legend.margin = margin(l = 0, unit = "cm")) +
    geom_rect(data = sigT, mapping = aes(xmin = V1, xmax = V2, fill = "yellow"),
              ymin = -Inf, ymax = Inf, color = "gray",
              alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE) +
    scale_fill_manual(values = "#FFF33D")
  if (ciBands) {
    p <- p + geom_line(aes(y = seL, group = Groups, color = Groups), linetype = 'dashed') +
      geom_line(aes(y = seU, group = Groups, color = Groups), linetype = 'dashed')
  }

  ## Plot diffs if included (switched order because I don't like the legend in middle)
  if (plotDiffs) {
    diffMat <- makePlotCI(diffList, alpha)
    diffMat$Time <- Time

    if (!dod) {
      dmain <- "Difference of Curves"
    } else {
      dmain <- "Difference of Difference Curves"
    }

    dp <- ggplot(diffMat, aes(Time, y)) +
      geom_line(size = 1) + ylab(bquote(paste(Delta, .(yylab), sep = " "))) +
      labs(x = xxlab, title = dmain) +
      geom_rect(data = sigT, mapping = aes(xmin = V1, xmax = V2, fill = "yellow"),
                ymin = -Inf, ymax = Inf, color = "gray",
                alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE) +
      scale_fill_manual(values = "#FFF33D")

    if (ciBands) {
      dp <- dp + geom_line(aes(y = seL), linetype = 'dashed', color = 'gray') +
        geom_line(aes(y = seU), linetype = 'dashed', color = 'gray')
    }

    grid.arrange(dp, p, ncol = 2, respect = TRUE, clip = "on")
    grid.arrange(dp, p, ncol = 2, respect = FALSE)
    return(invisible(list(diffPlot = dp, bootPlot = p)))
  } else {
    print(p)
    return(invisible(list(diffPlot = NULL, bootPlot = p)))
  }

}
