#' Summary for bdotsObj
#'
#' Provides summary information for bdotsObj
#'
#' @param object An object of class bdotsObj
#' @param ... not used
#'
#' @return Returns an object of class "bdotsSummary". There is some summarized
#' information included if assigned to an object, i.e., `summ <- summary(bdObj)`
#' then `str(summ)`
#' @import stats
#' @import nlme
#' @export
summary.bdotsObj <- function(object, ...) {
  bdObj <- object
  bdCall <- attr(bdObj, "call")
  subj <- bdCall$subject
  grps <- eval(bdCall$group)
  crvType <- attr(bdObj, "curveType")
  formula <- deparse1(attr(bdObj, "formula"))
  time <- attr(bdObj, "time")
  timeRange <- range(time)
  groups <- attr(bdObj, "groups")

  # cheap workaround to reduce size when split
  X <- attr(bdObj, 'X')$X
  #attr(bdObj, 'X')$X <- NULL

  ## Length of this is number of grp permutations fit
  # names of are names of those group permutations
  # values are number of each fit code
  grpSummary <- lapply(split(bdObj, by = grps, drop = TRUE), function(x) {
    mm <- coef(x)
    parMean <- colMeans(mm, na.rm = TRUE)
    vv <- var(mm, na.rm = TRUE)
    fitCount <- table(factor(x[['fitCode']], levels = 0:6))
    nobs <- length(unique(x[[subj]]))
    list(nobs = nobs, pars = parMean, varmat = vv, fitCount = fitCount)
  })
  mm <- coef(bdObj)
  parMean <- colMeans(mm, na.rm = TRUE)
  vv <- var(mm, na.rm = TRUE)
  fitCount <- table(factor(bdObj[['fitCode']], levels = 0:6))
  nobs <- nrow(bdObj)
  totalSummary <- list(nobs = nobs, pars = parMean,
                       varmat = vv, fitCount = fitCount)
  allSummary <- c(total = list(totalSummary), grpSummary)
  allSummary <- rev(allSummary)

  structure(.Data = list(curveType = crvType,
                         formula = formula,
                         groups = groups,
                         ntime = length(time),
                         timeRange = timeRange,
                         summaries = allSummary),
            class = "bdotsSummary",
            call = bdCall)

}


#' Print bdotsObj Summary
#'
#' Stuff
#'
#' @param x object to be printed
#' @param ... not used
#'
#' @details That's pretty much it. This is a print method, so there is likely
#' not much need to call it directly
#' @export
print.bdotsSummary <- function(x, ...) {
  cat("\nbdotsFit Summary\n\n")
  cat("Curve Type:", x$curveType, "\n")
  cat("Formula:", x$formula, "\n")
  cat("Time Range:", paste0("(", paste0(x$timeRange, collapse = ", "), ")"))
  cat(paste0(" [", x$ntime, " points]\n"))
  grpNames <- c(makeGroupNameVal(x$groups), "All Fits")
  cnts <- x$summaries
  cat("\n\n")
  for(i in seq_along(cnts)) {
    if (i != 1) cat("\n\n")
    #cat(names(cnts)[i], "\n\n")
    cat(grpNames[[i]], "\n")
    cat("Num Obs: ", cnts[[i]][['nobs']], "\n")
    cat("Parameter Values: \n")
    #cat(names(cnts[[i]][['pars']]), "\n", cnts[[i]][['pars']], "\n")
    print(cnts[[i]][['pars']])
    printFitCount(cnts[[i]])
  }

  # return invisibly
  invisible(x)
}


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

##########################

#' Summary for bdotsBootObj
#'
#' Provides summary information for bdotsBootObj
#'
#' @param object An object of class bdotsObj
#' @param ... Ignored for now
#'
#' @return Returns an object of class "bdotsBootSummary". There is some summarized
#' information included if assigned to an object, i.e., `summ <- summary(bdBootObj)`
#' then `str(summ)`
#' @export
summary.bdotsBootObj <- function(object, ...) {
  bdBootObj <- object
  ## Header info
  bdCall <- attr(bdBootObj, "call")
  alphastar <- bdBootObj[['adj.alpha']]
  sigTime <- bdBootObj[['sigTime']]
  rho <- bdBootObj[['rho']]
  dod <- bdBootObj[['dod']]
  paired <- bdBootObj[['paired']]
  curveGroup <- bdBootObj[['curveGroups']]
  formula <- deparse1(attr(bdBootObj, "bdObjAttr")[['formula']])
  curveType <- attr(bdBootObj, "bdObjAttr")[['curveType']]
  time <- attr(bdBootObj, "bdObjAttr")[['time']]
  timeRange <- range(time)

  padj_method <- match.arg(attr(bdBootObj, "call")[['p.adj']],
                           c("oleson", stats::p.adjust.methods))


  ## group specific info
  diffs <- bdBootObj[['diffs']]
  outerDiff <- diffs[['outerDiff']]
  outDiffVals <- curveGroup[[outerDiff]]

  ## we will always collect outerdiff info
  # that would just include outerDiff and pair status (already have above)
  if (dod) {
    innerDiff <- diffs[['innerDiff']]
    inDiffVals <- curveGroup[[innerDiff]]
    innerList <- setNames(vector("list", length = 2), outDiffVals)
    ## Can safely replace with with lapply
    for(gp in outDiffVals) {
      ll <- names(bdBootObj[['curveList']][[gp]])
      diffidx <- grep("diff", ll)
      ip <- bdBootObj[['curveList']][[gp]][[diffidx]][['paired']]
      groups <- list(groups = diffs, vals = ll[-diffidx])
      innerList[[gp]] <- list(groups = groups, paired = ip)
    }
  } else {
    innerList <- NULL
  }

  structure(.Data = list(formula = formula,
                         alphastar = alphastar,
                         sigTime = sigTime,
                         rho = rho,
                         dod = dod,
                         diffs = diffs,
                         curveGroup = curveGroup,
                         paired = paired,
                         innerList = innerList,
                         timeRange = timeRange,
                         ntime = length(time),
                         curveType = curveType,
                         padj_method = padj_method),
            class = "bdotsBootSummary",
            call = bdCall)

}



#' Print bdotsBoot Summary
#'
#' That's pretty much it. This is a print method, so there is likely
#' not much need to call it directly
#' @param x generic name, but this will be an object of bdotsBootSummary
#' @param ... ignored for now
#'
#' @export
print.bdotsBootSummary <- function(x, ...) {
  cat("\nbdotsBoot Summary\n\n")
  cat("Curve Type:", x$curveType, "\n")
  cat("Formula:", x$formula, "\n")
  cat("Time Range:", paste0("(", paste0(x$timeRange, collapse = ", "), ")"))
  cat(paste0(" [", x$ntime, " points]\n\n"))

  dod <- x[['dod']]
  cat("Difference of difference:", dod, "\n")
  cat("Paired t-test:", x[['paired']], "\n")
  if (dod) {
    cat("Outer Difference:", x[['diffs']][['outerDiff']], "\n")
    cat("Inner Difference:", x[['diffs']][['innerDiff']], "\n")
  } else {
    cat("Difference:", x[['diffs']][['outerDiff']], "\n")
  }
  cat("\n")
  cat("Alpha adjust method:", x$padj_method, "\n")
  cat("Adjusted Alpha:", x[['alphastar']], "\n")
  cat("Significant Intervals at adjusted alpha:\n")
  print(x[['sigTime']])

  ## Return the summary invisibly
  invisible(x)
}


#
# summary <- function(x, ...)
#   UseMethod("summary")
#












