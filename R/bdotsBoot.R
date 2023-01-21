#' Create bootstrapped curves from bdotsObj
#'
#' Creates bootstrapped curves and performs alpha adjustment. Can perform
#' "difference of difference" for nested comparisons
#'
#' @param formula See details.
#' @param bdObj An object of class 'bdotsObj'
#' @param Niter Number of iterations of bootstrap to draw
#' @param alpha Significance level
#' @param padj Adjustment to make to pvalues for significance. Will be able to
#' use anything from \code{p.adjust} function, but for now, just "oleson"
#' @param permutation Boolean indicating whether to use permutation testing rather
#' thank adjusting alpha to control FWER. WARNING: This option is very much in beta testing
#' and not recommended for general use. Also not available for paired tests or difference of difference
#' @param cores Number of cores to use in parallel. Default is zero, which
#' uses half of what is available.
#' @param ... not used
#'
#' @details The formula is the only tricky part of this. There will be a minor
#' update to how it works in the future. The three parts we will examine here
#' are Groups, the LHS, and the RHS. For all variable names, special characters
#' should be included with backticks, i.e., \code{`my-var`}
#'
#' ## Groups
#'
#' The Groups are the values input in \code{group} in the \code{bdotsFit} function,
#' which are columns of the dataset used. These will be denoted G_i
#' Within each group, we will designate the unique values within each group as v_j, ..., whereby
#' G_i(v_1, v_2) will designate unique two unique values within G_i. The possible
#' values of v_i will be implied by the group with which they are associated.
#'
#' For example, if we have groups \code{vehicle} and \code{color}, we could specify
#' that we are interested in all red cars and trucks with the expression
#' \code{vehicle(car, truck) + color(red)}.
#'
#' ## Formula
#'
#' ### Bootstrapped difference of curves
#'
#' This illustrates the case in which we are taking a simple bootstraped difference
#' between two curves within a single group
#'
#' If only one group was provided in \code{bdotsFit}, we can take the bootstrapped
#' difference between two values within the group with
#'
#' \code{y ~ Group1(val1, val2)}
#'
#' If more than two groups were provided, we must specify within which values of the
#' other groups we would like to compare the differences from Group1 in order to
#' uniquely identify the observations. This would be
#'
#' \code{y ~ Group1(val1, val2) + Group2(val1)}
#'
#' For example, bootstrapping the differences between cars and trucks when \code{color}
#' was provided as a second group, we would need \code{y ~ vehicle(car, truck) + color(red)}.
#'
#' ### Bootstrapped difference of difference curves
#'
#' This next portion illustrates the case in which we are interested in studying
#' the difference between the differences between two groups, which we will call
#' the innerGroup and the outerGroup following a nested container metaphor. Here,
#' we must use caution as the order of these differences matter. Using again the
#' vehicle example, we can describe this in two ways:
#'
#' \enumerate{
#'    \item We may be interested in comparing the difference between red trucks and cars (d_red) with
#'    the difference between blue trucks and cars (d_blue). In this case, we will be finding
#'    the difference between cars and trucks twice (one for blue, one for red). The
#'    vehicle type is the innerGroup, nested within the outerGroup, in this case, color.
#'    \item We may also be interested in comparing the difference between red trucks
#'    and blue trucks (d_truck) with the difference between red and blue cars (d_car).
#'    Here, innerGroup is the color and outerGroup is the vehicle
#' }
#'
#' As our primary object of interest here is not the difference in outcome itself, but the difference
#' of the outcome within two groups, the LHS of the formula is written
#' \code{diffs(y, Group1(val1, val2))}, where Group1 is the innerGroup. The RHS
#' is then used to specify the groups of which we want to take the outer difference of. The
#' syntax here is the same as above. Together, then, the formula looks like
#'
#' \code{diffs(y, Group1(val1, val2)) ~ Group2(val1, val2)}
#'
#' in the case in which only two grouping variables were provided to \code{bdotsFit}
#' and
#'
#' \code{diffs(y, Group1(val1, val2)) ~ Group2(val1, val2) + Group3(val1) + ...}
#'
#' is used to uniquely identify the sets of differences when three or more groups were provided.
#'
#' @return  Object of class 'bdotsBootObj'
#'
#' @examples
#' \dontrun{
#'
#' ## fit <- bdotsFit(cohort_unrelated, ...)
#'
#' boot1 <- bdotsBoot(formula = diffs(Fixations, LookType(Cohort, Unrelated_Cohort)) ~ Group(50, 65),
#'                    bdObj = fit,
#'                    N.iter = 1000,
#'                    alpha = 0.05,
#'                    p.adj = "oleson",
#'                    cores = 4)
#'
#' boot2 <- bdotsBoot(formula = Fixations ~ Group(50, 65) + LookType(Cohort),
#'                    bdObj = fit,
#'                    N.iter = 1000,
#'                    alpha = 0.05,
#'                    p.adj = "oleson",
#'                    cores = 4)
#' }
#'
#' @import data.table
#' @export
bdotsBoot <- function(formula,
                      bdObj,
                      Niter = 1000,
                      alpha = 0.05,
                      padj = "oleson",
                      permutation = FALSE,
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
  curveFun <- makeCurveFun(bdObj)

  ## Next, we want to get a bootstrapped distribution for each of the groups
  splitGroups <- split(bdObj, by = c(innerDiff, outerDiff)) # ok even if null


  ## This cannot stay here forever, its a super secret shortcut for just right now
  existsC <- exists("collinshortcut")
  takeshortcut <- ifelse(existsC, collinshortcut, FALSE)
  if (takeshortcut) {
    # just skip doing this
    curveList <- NULL
    ip <- NULL
  } else {
    # do the normal thing
    ## Make this parallel (maybe later)
    if (Sys.info()['sysname'] == "Darwin") {
      cl <- makePSOCKcluster(cores, setup_strategy = "sequential")
    } else {
      cl <- makePSOCKcluster(cores)
    }
    invisible(clusterEvalQ(cl, {library(bdots)}))
    # this needs to happen anyways, this is bootstrapped distributions
    ## WE NEED TO INDICATE IF PAIRED HERE
    groupDists <- parLapply(cl, splitGroups, getBootDist, b = Niter)

    stopCluster(cl)

    ## This is where we construct inner/outer groups
    # (ideally matching old bdots, at least for now)
    curveList <- createCurveList(groupDists, prs, splitGroups) # this is whats creates diff
    ip <- curveList[['diff']][['paired']]
  }



  # Determine first if we are doing difference of differences
  dod <- ifelse(is.null(innerDiff), FALSE, TRUE)

  ## Currently permutation doesnt exist for dod
  if (dod & permutation) {
    warning("Permutation testing does not yet work for difference of difference analysis. Switching to padj='oleson' instead")
    permutation <- FALSE
  }

  ## And here we determine significant regions
  if (permutation) {
    message("WARNING: permutation testing is work in progress and limited in scope")
    # do permutation
    res <- permTest(splitGroups, prs, alpha = alpha, P = Niter) # in permutation.R
    obsT <- res[['obst']]
    nullT <- res[['nullt']]

    # Do I want to return null distribution and vectors?
    # Do I compute any pvalues here? For now, just significant time
    time <- attr(bdObj, "time")
    sigTime <- bucket(res[["sigIdx"]], time)

    # Remove things we don't need (temporary)
    pval <- NULL
    rho <- NULL
    alphastar <- NULL
    adjpval <- NULL
  } else {
    res <- alphaAdjust(curveList, padj, alpha, cores)
    pval <- res[['pval']]
    rho <- res[['rho']]
    alphastar <- res[['alphastar']]
    adjpval <- res[['adjpval']]
    time <- attr(bdObj, "time")
    sigTime <- bucket(adjpval <= alpha, time)

    # Remove things from perm (this is temporary fix)
    obsT <- NULL
    nullT <- NULL
  }

  ## maybe consider keeping ancillary data in list, i.e., res.
  structure(class = "bdotsBootObj",
            .Data = list(curveList = curveList,
                         alpha = alpha,
                         adjalpha = alphastar,
                         adjpval = adjpval,
                         sigTime = sigTime,
                         rho = rho,
                         paired = ip,
                         diffs = c("outerDiff" = outerDiff,
                                   "innerDiff" = innerDiff), # could replaces this with `prs`
                         curveGroups = curveGrps,
                         dod = dod,
                         curveFun = curveFun),
            call = match.call(),
            bdObjAttr = attributes(bdObj))
}

