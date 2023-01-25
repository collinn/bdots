## This file contains:
# function to subset bdObj from parser
# function to create distribution for each group (bootstrap)
# function to create curveList (along with inner/outer functions)
# alpha adjustment function
# bucket function

##------------------------------------------------------------------------------

## Subset bdObj based on groups being compared
bootGroupSubset <- function(l, bdObj) {
  subnames <- l[["subnames"]]
  subargs  <- l[["subargs"]]
  resNames <- l[['resNames']]

  bd <- bdObj[, c(resNames, subnames), with = FALSE]

  ## This will iteratively subset itself
  for(i in seq_along(subnames)) {
    ss_vec <- bd[[subnames[i]]] %in% subargs[[i]]
    bd <- bd[ss_vec, ]  #subsets multiple times
  }

  ## I'm still going to keep order
  bd[, c(resNames[1], subnames, resNames[-1]), with = FALSE]
}


#' Create old bdots curve list
#'
#' @param x list of group distributions
#' @param prs named list of parsed call
#' @param splitGroups splitGroups
createCurveList <- function(x, prs, splitGroups) {

  ## Whats what
  innerDiff <- prs[["innerDiff"]]
  outerDiff <- prs[["outerDiff"]]
  curveGrps <- setNames(prs[['subargs']], prs[['subnames']])

  ## Ok, are we dod or not?
  dod <- !is.null(innerDiff)

  if (!dod) {
    ## Really, this should be called outer diff list since there is no inner group
    diffList <- makeInnerDiffList(x, splitGroups)
    res <- structure(.Data = setNames(c(x, list(diffList)), c(names(x), "diff")),
                     class = c("innerGroupCurveList", "groupCurveList"))
  } else {
    xnames <- vapply(strsplit(names(x), "\\."), `[[`, character(1), 2)

    ## we subset by innerDiff because we want like terms together
    gname <- paste0("^", curveGrps[[outerDiff]][1], "$")
    idx <- grepl(gname, xnames)

    ## List of soon-to-be innerLists
    xx <- setNames(list(x[idx], x[!idx]), c(xnames[idx][1], xnames[!idx][1]))

    ## Make innerDiffList
    diffList <- Map(makeInnerDiffList,
                    curveList = xx,
                    oP = list(splitGroups[idx], splitGroups[!idx]))

    outerlists <- Map(function(x, y) {
      nn <- vapply(strsplit(names(x), "\\."), `[[`, character(1), 1)
      structure(.Data = setNames(c(x, list(y)), c(nn, "diff")),
                class = c("innerGroupCurveList", "groupCurveList"))
    }, x = xx, y = diffList)

    obj <- rbindlist.bdObjList(splitGroups)
    obj <- split.bdotsObj(obj, by = outerDiff, drop = TRUE)
    diffList <- makeOuterDiffList(outerlists, obj)

    res <- structure(.Data = setNames(c(outerlists, list(diffList)),
                                      c(names(outerlists), "diff")),
                     class = c("outerGroupCurveList","groupCurveList"))
  }
  return(res)
}

##------------------------------------------------------------------------------
## Make diffList from curveList
makeInnerDiffList <- function(curveList, oP) {
  diffList <- Map(function(x, y) {
    y - x
  }, curveList[[1]], curveList[[2]])

  if (ip <- isPaired(oP)) {
    diffList$sd <- apply(diffList$curveMat, 2, sd) # this is correct
    diffList$n <- nrow(oP[[1]]) - 1L
  } else {
    diffList$sd <- nopairSD(curveList)
    diffList$n <- sum(vapply(oP, nrow, numeric(1))) - 2L
  }
  diffList$paired <- ip
  structure(.Data = diffList,
            class = c("bdInnerDiffList", "bdDiffList"))
}


## Join and take diff of two inner diffs
makeOuterDiffList <- function(res, obj) {
  res <- unlist(res, recursive = FALSE)
  idx <- grep("diff", names(res))
  if (length(idx) != 2L) stop("something weird in curveBooter. Contact author")

  ## diff of diff (length one list)
  diffList <- Map(function(x, y) {
    Map(function(a, b) {
      a - b
    }, x, y)
  }, res[idx[1]], res[idx[2]])

  ## Map returns a lenght 1 list
  diffList <- diffList[[1]]

  ## snap, we can
  if (ip <- isPaired(obj)) {
    diffList$sd <- apply(diffList$curveMat, 2, sd)
    diffList$n <- nrow(obj[[1]]) - 1L
  } else {
    diffList$sd <- nopairSD(res[idx])
    diffList$n <- sum(vapply(obj, nrow, numeric(1L))) - 2L
  }
  diffList$paired <- ip
  structure(.Data = diffList,
            class = c("bdOuterDiffList", "bdDiffList"))
}



##------------------------------------------------------------------------------

## alphaAdjust
## curveList - returned from curveBooter
## group - either Group name or Group value, i.e., LI = LI.M/LI.W or W = LI.W/TD.W
## For right now, group can only be group name (since I would need to recalculate diff for group value)
alphaAdjust <- function(curveList, p.adj = "oleson", alpha = 0.05, cores, group = NULL) {
  if (is.null(group)) {
    curve <- curveList[['diff']]
  } else {
    idx <- grep(group, names(curveList))
    if (length(idx) == 0) stop("Invalid group name")
    d_idx <- grep(paste0(group, "\\.diff"), names(curveList[[idx]]))
    curve <- curveList[[idx]][[d_idx]]
  }
  tstat <- curve[["fit"]]/curve[['sd']]
  pval <- 2 * (1 - pt(abs(tstat), df = curve[['n']]))

  p.adj <- match.arg(p.adj, c("oleson", stats::p.adjust.methods))
  rho <- ar1Solver(tstat)

  adjpval <- p_adjust(pval, p.adj, length(tstat), alpha,
                      df = curve[['n']], rho = rho, cores = cores)
  alphastar <- attr(adjpval, "alphastar")

  list(pval = pval, adjpval = adjpval, alphastar = alphastar, rho = rho)
}

## Determines time regions where difference significant
# sig are boolean values
bucket <- function(sig, time) {
  if (sum(sig) == 0) return(NULL)
  rr <- rle(sig)
  rr$idx <- cumsum(rr$lengths)
  mm <- rr$values*(rr$idx) + (1 - rr$values)*(rr$idx + 1)

  ## Alt condition if first region true
  if (rr$values[1] == 1) {
    mm <- c(1, mm)
  }
  # If ends on false
  if (rr$values[length(rr$values)] == 0) {
    mm <- mm[1:(length(mm) - 1)]
  }
  matrix(time[mm], ncol = 2, byrow = TRUE)
}





