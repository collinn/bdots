

#' Correlation with fixed value in bdots
#'
#' Find the correlation of a fixed value with the bdots fitted
#' curves at each time point
#'
#' @param bdObj Object of class `bdotsObj`
#' @param val Character string of fixed value for correlation in dataset from `bdotsFit`
#' @param ciBands Boolean for including confidence intervals
#' @param method Arguments for `cor` or `cor.test`. The default option us `method = "pearson"`
#'
#' @export
bdotsCorr <- function(bdObj, val, ciBands = FALSE, method = "pearson") {

  Group <- Correlation <- lower <- upper <- . <- NULL

  X <- attr(bdObj, "X")$X

  if (!(val %in% names(X))) stop("val must be a column of dataset used in bdotsFit")
  if (!is.numeric(X[[val]])) stop("val must be a numeric column for correlation")

  id_cols <- getIdentifierCols(bdObj)
  tt <- split(bdObj, by = id_cols[-1])
  corMats <- lapply(tt, getFitCorforGroups, val = val,
                    ciBands = ciBands, method = method)

  if (!ciBands) {
    zz <- as.data.table(corMats)
    zz[['time']] <- attr(bdObj, "time")
    rr <- melt(zz, measure.vars = names(zz)[names(zz)!="time"],
               value.name = "Correlation", variable.name = "Group")

    gg <- as.data.table(tstrsplit(rr[["Group"]], split= "\\."))
    colnames(gg) <- paste0("Group", 1:(ncol(gg)))

    rr[, Group := gsub("\\.", " ", Group)]
    rr <- cbind(rr[, .(time, Correlation, Group)], gg)
  } else {
    zz <- lapply(names(corMats), function(nn) {
      ww <- corMats[[nn]]
      ww <- as.data.table(t(ww))
      colnames(ww) <- c("Correlation", "lower", "upper")
      ww[['time']] <- attr(bdObj, "time")
      ww[['Group']] <- nn

      gg <- as.data.table(tstrsplit(ww[["Group"]], split= "\\."))
      colnames(gg) <- paste0("Group", 1:(ncol(gg)))
      ww[, Group := gsub("\\.", " ", Group)]
      ww <- cbind(ww[, .(time, Correlation, lower,  upper, Group)], gg)
      ww
    })
    rr <- data.table::rbindlist(zz)
  }
  structure(.Data = rr, class = c("bdotsCorrObj", "data.table", "data.frame"),
            ciBands = ciBands)
}

#' Get Fit Correlations
#'
#' Helper function for finding correlation of fixed value and
#' fitted values within group
#'
#' @param x A split object of class `bdObj` split by identifiers
#' @param val Fixed value from dataset
#' @param ciBands boolean for including cibands
#' @param method method for correlation function
getFitCorforGroups <- function(x, val, ciBands = FALSE, method = "pearson") {
  l <- length(attr(x, "time"))
  n <- nrow(x)
  m <- matrix(NA, nrow = n, ncol = l)
  ## each row is subject, column is time point
  for (i in seq_len(n)) {
    m[i, ] <- fitted.values(x[i, ]$fit[[1]])
  }
  dt <- getSubX(x)
  id_cols <- getIdentifierCols(x)

  dt <- unique(dt[, c(val, id_cols), with = FALSE])

  v <- dt[[val]]

  ## possibly replace with cor.test if we want confidence intervals
  if (ciBands) {
    cormat <- apply(m, 2, function(y) {
      qq <- cor.test(x = y, y = v, method = method)
      c(cor = qq$estimate, lower = qq$conf.int[1], upper = qq$conf.int[2])
    })
    colnames(cormat) <- attr(x, "time")
  } else {
    cormat <- apply(m, 2, function(y) cor(y, v, method = method))
    cormat <- setNames(cormat, attr(x, "time"))
  }
  cormat
}


#' Plots for bdotsCorr
#'
#' Plots correlation of fixed value wtih fitted curves over time
#'
#' @param x object of class `bdotsCorrObj`
#' @param ciBands boolean. Whether or not to include confidence intervals in plots.
#' Must have been selected in `bdotsCorr`
#' @param window A length 2 numeric vector with start and end points
#' for the  plotting window
#' @param ... super secret, don't use
#'
#' @export
plot.bdotsCorrObj <- function(x, ciBands = FALSE, window = NULL, ...) {

  ggplot2::theme_set(theme_bw())
  lower <- upper <- Correlation <- NULL

  if (!is.null(window)) {
    if (!is.numeric(window) | length(window) != 2) stop("Window must be length 2 numeric")
    wind <- seq.int(from = window[1], to = window[2], by = 1)
    x <- x[time %in% wind, ]
  }

  ci <- attr(x, "ciBands")
  if (ci & ciBands) {
    ggplot(x, aes(x = time, y = Correlation)) +
      geom_line() +
      geom_line(aes(x = time, y = lower), linetype='dashed') +
      geom_line(aes(x = time, y = upper), linetype='dashed') +
      facet_wrap(~Group)
  } else {
    ggplot(x, aes(x = time, y = Correlation)) +
      geom_line() +
      facet_wrap(~Group)
  }
}



