## Need summary functions for bdotsObj (from bdotsFit) and bdotsBootObj
## BUT here's the thing. We need a summary that constructs the information
# but then a print method specifically for that summary

## Why is the above true? Because the summary function *should* perform
# all of the summary methods that a user may want. Consequently, the use of
# this information is what should determine the value of the object returned
# by the summary function. How it is actually displayed, as visual information
# not to be further processed, is why we employ the print method as well

## bdotsObj
# Number of unique subjects
# Subjects in each group
# classification of fits
# info about parameters

summary.bdotsObj <- function(bdObj, ...) {
  bdCall <- attr(bdObj, "call")
  subj <- bdCall$subject
  grps <- eval(bdCall$group)
  crvType <- attr(bdObj, "curveType")
  formula <- deparse1(attr(bdObj, "formula"))
  time <- attr(bdObj, "time")
  timeRange <- range(time)

  # cheap workaround to reduce size when split
  X <- attr(bdObj, 'X')
  on.exit(attr(bdObj, 'X') <- X)
  attr(bdObj, 'X') <- NULL

  ## Length of this is number of grp permutations fit
  # names of are names of those group permutations
  # values are number of each fit code
  grpSummary <- lapply(split(bdObj, by = grps, drop = TRUE), function(x) {
    mm <- coef(x)
    parMean <- colMeans(mm)
    vv <- var(mm)
    fitCount <- table(x[['fitCode']])
    nobs <- length(unique(x[[subj]]))
    list(nobs = nobs, pars = parMean, varmat = vv, fitCount = fitCount)
  })
  mm <- coef(bdObj)
  parMean <- colMeans(mm)
  vv <- var(mm)
  fitCount <- table(bdObj[['fitCode']])
  nobs <- length(unique(bdObj[[subj]]))
  totalSummary <- list(nobs = nobs, pars = parMean,
                       varmat = vv, fitCount = fitCount)
  allSummary <- c(total = list(totalSummary), grpSummary)

  structure(.Data = list(curveType = crvType,
                         formula = formula,
                         groups = names(grpSummary),
                         ntime = length(time),
                         timeRange = timeRange,
                         summaries = allSummary),
            class = "bdotsSummary",
            call = bdCall)

}

#rr <- summary.bdotsObj(bdObj)

print.bdotsSummary <- function(x, ...) {
  cat("\nbdotsFit Summary\n\n")
  cat("Curve Type:", x$curveType, "\n")
  cat("Formula:", x$formula, "\n")
  cat("Time Range:", paste0("(", paste0(x$timeRange, collapse = ", "), ")"))
  cat(paste0(" [", x$ntime, " points]\n"))
  cnts <- x$summaries
  for(i in seq_along(cnts)) {
    cat(names(cnts)[i], "\n\n")
    cat("Num Obs: ", cnts[[i]][['nobs']], "\n")
    cat("Parameter Values: \n")
    cat(names(cnts[[i]][['pars']]), "\n", cnts[[i]][['pars']], "\n")
    printFitCount(cnts[i])
  }
}

## This takes an actual summaries value
printFitCount <- function(x) {
  x <- x[['fitCount']]
  printLine <- c(paste("AR1,       0.95 <= R2        --", x[1], "\n"),
                 paste("AR1,       0.80 < R2 <= 0.95 --", x[2], "\n"),
                 paste("AR1,       R2 < 0.8          --", x[3], "\n"),
                 paste("Non-AR1,   0.95 <= R2        --", x[4], "\n"),
                 paste("Non-AR1,   0.8 < R2 <= 0.95  --", x[5], "\n"),
                 paste("Non-AR1,   R2 < 0.8          --", x[6], "\n"),
                 paste("No Fit                       --", x[7], "\n"))
  cat("########################################\n")
  cat("############### FITS ###################\n")
  cat("########################################\n")
  for(i in seq_along(printLine)) {
    cat(printLine[i])
  }
}

## From original
cat("########################################\n")
cat("############### FITS ###################\n")
cat("########################################\n")
cat(paste("AR1,       0.95 <= R2        --", ar1.good, "\n"))
cat(paste("AR1,       0.80 < R2 <= 0.95 --", ar1.ok, "\n"))
cat(paste("AR1,       R2 < 0.8          --", ar1.bad, "\n"))
cat(paste("Non-AR1,   0.95 <= R2        --", nonar1.good, "\n"))
cat(paste("Non-AR1,   0.8 < R2 <= 0.95  --", nonar1.ok, "\n"))
cat(paste("Non-AR1,   R2 < 0.8          --", nonar1.bad, "\n"))
cat(paste("No Fit                       --", nofit, "\n"))
cat("########################################\n\n")
