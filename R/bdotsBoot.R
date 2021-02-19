#' Create bootstrapped curves from bdotsObj
#'
#' Creates bootstrapped curves and performs alpha adjustment. Can perform
#' "difference of difference" for nested comparisons
#'
#' @param formula See details.
#' @param bdObj An object of class 'bdotsObj'
#' @param N.iter Number of iterations of bootstrap to draw
#' @param alpha Significance level
#' @param p.adj Adjustment to make to pvalues for significance. Will be able to
#' use anything from \code{p.adjust} function, but for now, just "oleson"
#' @param cores Number of cores to use in parallel. Default is zero, which
#' uses half of what is available.
#' @param ... not used
#'
#' @details The formula is the only tricky part of this. There will be a minor
#' update to how it works in the future.  The three parts we will examine here
#' are Groups, the LHS, and the RHS. This is incompelete below, but is detailed
#' in the original vignette sent out
#'
#' ## Groups
#'
#' To be confusing, we use "groups" in two senses here. The first use, which we
#' will denote with captial "G" Group, indicates the column name of the dataset
#' representing one of the variables by which we will distinguish categories of
#' observations. Lowercase "g" group will be the values in this column.
#'
#' ## LHS
#'
#' If we are comparing two groups, we will just use the name of the outcome, i.e.,
#' `y ~` .
#'
#' @return  Object of class 'bdotsBootObj'
#'
#' @import data.table
#' @export

bdotsBoot <- function(formula,
                      bdObj,
                      N.iter = 1000,
                      alpha = 0.05,
                      p.adj = "oleson",
                      cores = 0, ...) {

  if (cores < 1) cores <- detectCores()/2

  if (any(bdObj[['fitCode']] == 6)) {
    warning("Some observations had NULL gnls fits. These and their pairs will be removed")
    bdObj <- bdRemove(bdObj, fitCode = 6, removePairs = TRUE)
  }

  prs <- bootParser(formula, bdObj)
  bdObj <- bootGroupSubset(prs, bdObj)
  innerDiff <- prs[["innerDiff"]] # unnamed char vec
  outerDiff <- prs[["outerDiff"]] # unnamed char vec

  curveGrps <- setNames(prs[['subargs']], prs[['subnames']])

  makeCurveFun <- function(bdObj) {
    time <- attr(bdObj, "call")[['time']]
    f_bod <- attr(bdObj, "formula")[[3]]
    f_args <- c(colnames(coef(bdObj)), time)
    f_args <- setNames(as.pairlist(rep("", length(f_args))), f_args)
    eval(call("function", f_args, f_bod), parent.frame())
  }
  curveFun <- makeCurveFun(bdObj)

  ## This is either an innerDiff or outerDiff list (check class)
  curveList <- curveBooter(bdObj,
                           outerDiff = outerDiff,
                           innerDiff = innerDiff,
                           N.iter = N.iter,
                           curveFun = curveFun)
  ip <- curveList[['diff']][['paired']] # paired?

  res <- alphaAdjust(curveList, p.adj, alpha, cores)
  pval <- res[['pval']]
  rho <- res[['rho']]
  alphastar <- res[['alphastar']]
  adjpval <- res[['adjpval']]
  time <- attr(bdObj, "time")
  sigTime <- bucket(pval <= alphastar, time)
  dod <- ifelse(is.null(innerDiff), FALSE, TRUE)

  ## maybe consider keeping ancillary data in list, i.e., res.
  structure(class = "bdotsBootObj",
            .Data = list(curveList = curveList,
                         alpha = alpha,
                         adj.alpha = alphastar,
                         adj.pval = adjpval,
                         sigTime = sigTime,
                         rho = rho,
                         paired = ip,
                         # also this is duplicated in curveGroups, i.e., diffs == names(curveGroups)
                         # however, it does not  specify between inner and outer group, so idk
                         diffs = c("outerDiff" = outerDiff,
                                   "innerDiff" = innerDiff), # could replaces this with `prs`
                         curveGroups = curveGrps,
                         dod = dod,
                         curveFun = curveFun),
            call = match.call(),
            bdObjAttr = attributes(bdObj))
}








