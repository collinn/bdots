

#' Write fits from \code{bdotsBoot} to csv file
#'
#' The function is used to write out columns for each group for which a curve
#' was bootstrapped
#'
#' @param bootObj An object of class \code{bdotsBootObj}
#' @param file file name to write out csv
#' @param alpha alpha level for upper/lower CI
#' @param ... Other arguments passed to \code{data.table::fread}
#'
#' @details This is potentially useful for constructing plots in a separate
#' application. There is an additional column, \code{Significant} indicating if
#' a particular time point was considered significant between the difference curves.
#' For difference of difference objects, this only indicates significance for the
#' outer difference.
#'
#' @export
writeCSV <- function(bootObj, file, alpha = 0.05, ...) {

  ## Get time column
  time <- attr(bootObj, "bdObjAttr")$time
  time <- matrix(time, ncol = 1)
  colnames(time) <- attr(bootObj, "bdObjAttr")$call$time

  ## Add significant indicator for differences
  sigMat <- rep(0, nrow(time))
  sigMat <- 0 * time
  colnames(sigMat) <- "Significant"
  st <- bootObj$sigTime
  for (i in seq_len(nrow(st))) {
    s1 <- which(time == st[i, 1])
    s2 <- which(time == st[i, 2])
    idx <- `:`(s1, s2)
    sigMat[idx, 1] <- 1
  }
  startMat <- cbind(time, sigMat)

  ## Return 'observed' values?
  # ...

  cl <- bootObj$curveList
  if (inherits(cl, "outerGroupCurveList")) {
    main_diff <- cl[['diff']] # don't like using the name here
    cl[['diff']] <- NULL
    cl <- unlist(cl, use.names = TRUE, recursive = FALSE)
    cl <- cl
    cl[['diff']] <- main_diff
  }

  fitMat <- lapply(cl, makePlotCI, alpha = alpha)
  fitMat <- Map(function(x, y) {
    nn <- paste(y, c("Lower CI", "Fit", "Upper CI"), sep = " - ")
    colnames(x) <- nn
    x[, c(2,1,3)]
  }, fitMat, names(cl))

  out <- Reduce(cbind, fitMat, init = startMat)


  out <- data.table::data.table(out)
  data.table::fwrite(out, file, ...)
}
