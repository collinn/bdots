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
    plotList[[i]] <- ggplot(cList[[i]], aes(value, group = parameter, color = 'gray')) +
      geom_histogram(bins = 8) +
      facet_wrap(~parameter, scales = "free_x") +
      scale_color_manual(values = "white") +
      ggtitle(names(cList)[i], subtitle = paste0("page ", i, " of ", length(cList))) +
      theme(legend.position = "none") + xlab("Value") + ylab("")

    print(plotList[[i]])
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

  ## Eventually we will delete this
  X <- attr(bdObj, "X")$X
  # dfname <- deparse1(bdCall$data)
  # if (is.null(X) & exists(dfname)) {
  #   X <- get(dfname)
  # } else if (is.null(X) & !exists(dfname)) {
  #   stop("Cannot find dataset used to construct bdObj, please pass as argument")
  # }

  ## Eventually we will delete this logic
  gridSize <- ifelse(is.null(gridSize), 2, gridSize)

  if (nrow(bdObj) < 4 & gridSize != "refit") gridSize <- 1

  # if (is.null(gridSize)) {
  #   gridSize <- 2
  #   par(mfrow = c(gridSize, gridSize))
  # } else if (gridSize == "refit") {
  #   par(mfrow = c(1, 2))
  # } else {
  #   par(mfrow = c(gridSize, gridSize))
  # }

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

    } else {
      df <- as.data.table(cbind(time, obsY))
    }
    # df2 <- melt(df, id.vars = "time")
    # df2$lty <- "dashed"
    # df2[variable != "obsY", ]$lty <- "solid"
    #
    # df2$clr <- "tomato"
    # df2[variable != "obsY", ]$clr <- "steelblue"
    #
    # df2[, `Curves` := variable]
    # df2$Curves <- "Fit"
    # df2[variable == "obsY", ]$Curves <- "Observed"

    df2 <- melt(df, id.vars = "time", variable.name = "Curves")
    df2[, `:=`(lty = ifelse(Curves == "obsY", "dashed", "solid"),
               clr = ifelse(Curves == "obsY", "tomato", "steelblue"),
               Curves = ifelse(Curves == "obsY", "Observed", "Fit"))]

    ## Janky fix for update. Should just make a separate for refits
    if (gridSize == "refit") {
      if (i == 1) {
        title <- paste0("Original Fit", "\nfitCode = ", code, ", R2 = ", r2)
      } else {
        title <- paste0("Updated Fit", "\nfitCode = ", code, ", R2 = ", r2)
      }
    } else {
      title <- paste0(paste0(obs, collapse = " "),
                      "\nfitCode = ", code, ", R2 = ", r2, collapse = "")
    }

    ## This can't be the way to do it
    y <- NULL; variable <- NULL; value <- NULL; clr <- NULL; lty <- NULL
    Curves <- NULL;

    #colors <- c("#56B4E9","#D55E00")

    ## This thing is kinda gross
    plotList[[i]] <- ggplot(df2, aes(time, value, color = Curves, linetype = Curves)) +
      geom_line(size = 1) + theme_bw() + ggtitle(title) + xlab(tname) + ylab(yname) +
      theme(legend.position = "bottom") +
      scale_color_manual(values = c("#56B4E9","#D55E00"))
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

## Generate confidence intervals around fit values
#' @import data.table
makePlotCI <- function(cl, alpha = 0.05, ...) {
  tv <- stats::qt(1 - alpha / 2, cl[['n']] - 1) # this will not be correct for diff, as 'n' will be wrong.
  fit <- cl[['fit']]
  sd <- cl[['sd']]
  mm <- data.table(seL = fit - sd * tv,
                   y   = fit,
                   seU = fit + sd * tv)
}


#' Plot for object of class bdotsBootObj
#'
#' Allows a number of different but also unstable option for plotting an object
#' of class bdotsBoot
#'
#' @param x An object of class bdotsBootObj
#' @param alpha Significance level for plotting confidence intervals.
#' @param ciBands Boolean indicating whether or not to include confidence intervals
#' around fitted curves
#' @param plotDiffs Boolean to plot difference curve
#' @param group Specify group to plot if difference of difference was used. The
#' user can also subset the bdotsBootObj prior to plotting. Currently not used
#' @param ... ignore for now, but will eventually allow plot parameters
#'
#' @returns List of ggplot objects, which may be helpful if the margins are weird
#'
#' @details This plot function is also a bit unstable and is expected to change
#' @export
plot.bdotsBootObj <- function(x, alpha = NULL, ciBands = TRUE, plotDiffs = TRUE, group = NULL, ...) {

  # value used in original call
  if (is.null(alpha)) alpha <- x$alpha

  ## Not used for now
  group <- NULL

  if (x$dod & !is.null(group)) {
    x <- subset(x, group, adjustAlpha = alpha) # this is a method
    plotInnerGroup(x, alpha, ciBands, plotDiffs, ...)
  } else if (x$dod) {
    plotOuterGroup(x, alpha, ciBands, plotDiffs, ...)
  } else {
    plotInnerGroup(x, alpha, ciBands, plotDiffs, ...)
  }
}

## For plotting diff of diff case
plotOuterGroup <- function(bdBootObj, alpha = 0.05, ciBands, plotDiffs, ...) {
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
plotInnerGroup <- function(bdBootObj, alpha = 0.05, ciBands, plotDiffs, dod = FALSE, ...) {

  ## Custodial tasks
  ggplot2::theme_set(theme_bw())
  Groups <- V1 <- V2 <- parameter <- seL <- seU <- value <- y <- NULL
  Time <- attr(bdBootObj, "bdObjAttr")[["time"]]
  bdCall <- attr(bdBootObj, "bdObjAttr")[['call']]
  yylab <- bdCall[['y']]
  xxlab <- bdCall[['time']]
  mmain <- ifelse(!dod, "Bootstrapped Fits", "Bootstrapped Differences")

  ## Build our datasets for plots
  cList <- bdBootObj[['curveList']]
  diffList <- cList[['diff']]
  cList <- cList[setdiff(names(cList), "diff")]

  ## Likely no chance this will affect a real group name, but need for legend
  cvMat <- lapply(cList, makePlotCI, alpha)
  cvMat <- data.table::rbindlist(cvMat, idcol = 'Groups')
  cvMat[, Groups := gsub("\\.diff$", "", Groups)]
  cvMat$Time <- rep(Time, times = length(cList))

  ## sigTime
  sigT <- data.table(bdBootObj$sigTime)

  ## Construction of plot
  p <- ggplot(cvMat, aes(Time, y)) +
    geom_line(aes(group = Groups, color = Groups), size = 1) +
    labs(x = xxlab, y = yylab, title = mmain) +
    theme(legend.margin = margin(l = 0, unit = "cm"), legend.position = "bottom") +
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

    dmain <- ifelse(!dod, "Difference of Curves", "Difference of Difference Curves")

    dp <- ggplot(diffMat, aes(Time, y)) +
      geom_line(size = 1) + ylab(bquote(paste(Delta, .(yylab), sep = " "))) +
      labs(x = xxlab, title = dmain) +
      geom_rect(data = sigT, mapping = aes(xmin = V1, xmax = V2, fill = "yellow"),
                ymin = -Inf, ymax = Inf, color = "gray",
                alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE) +
      scale_fill_manual(values = "#FFF33D") +
      theme(plot.margin = unit(c(5.5,5.5,34,5.5), "points"))

    if (ciBands) {
      dp <- dp + geom_line(aes(y = seL), linetype = 'dashed', color = 'gray') +
        geom_line(aes(y = seU), linetype = 'dashed', color = 'gray')
    }

    grid.arrange(p, dp, ncol = 2, respect = FALSE)

  } else {
    dp <- NULL
    print(p)
  }

  return(invisible(list(bootPlot = p, diffPlot = dp)))
}
