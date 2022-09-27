## curveBooter
# this function is dangerously long
# and potentially complicated
# need to revist for potential simplication.

curveBooter <- function(Obj, outerDiff, innerDiff = NULL, N.iter, curveFun) {

  ## This only happens IF there is both inner/outer diff
  if (!is.null(innerDiff)) {
    obj <- split.bdotsObj(Obj, by = outerDiff, drop = TRUE)
    res <- lapply(obj, curveBooter, outerDiff = innerDiff,
                  N.iter = N.iter, curveFun = curveFun)

    diffList <- makeOuterDiffList(res, obj)

    return(structure(.Data = setNames(c(res, list(diffList)),
                               c(names(res), "diff")),
                     class = c("outerGroupCurveList","groupCurveList")))
  }

  ## Determine correlation matrix, if paired
  oP <- split.bdotsObj(Obj, by = outerDiff, drop = TRUE)
  if (ip <- isPaired(oP)) {
    cmat <- lapply(oP, coef)
    corMat <- do.call("cor", setNames(cmat, c("x", "y")))
  } else {
    corMat <- NULL
  }

  ## Bootstrap values
  if (!is.null(corMat)) {
    outDiffL <- split.bdotsObj(Obj, by = attr(Obj, "call")[['subject']], drop = TRUE)
    bootPars <- lapply(outDiffL, bdotsBooter, N.iter, corMat)
    meanMat <- parMatSplit(Reduce(`+`,  bootPars)/length(bootPars))
  } else {
    outDiffL <- lapply(oP, split.bdotsObj, by = attr(Obj, "call")[['subject']], drop = TRUE)
    meanMat <- lapply(outDiffL, function(x) {
      bootPars <- lapply(x, bdotsBooter, N.iter, corMat)
      meanMat <- Reduce(`+`,  bootPars)/length(bootPars)
    })
  }

  ## class bdCurveList
  curveList <- makeCurveList(meanMat, curveFun, oP)

  ## Class bdDiffList
  diffList <- makeInnerDiffList(curveList, oP)

  structure(.Data = setNames(c(curveList, list(diffList)),
                             c(unique(Obj[[outerDiff]]), "diff")),
            class = c("innerGroupCurveList", "groupCurveList"))
}

###------------------------------------------------

## This function is responsible for taking the mean parameter matrix
# from each group over N.iter iterations. Along with time and the function
# specifying the curve to be fit, this returns a length 2 list, one for each
# of the groups being fit

# takes list of meanMatrix for each group (from bdotsBooter),
# and a numeric vec with timeName attributes (from original call)
makeCurveList <- function(meanMat, curveFun, oP) {
  time <- attr(oP[[1]], "time")
  timeName <- attr(oP[[1]], "call")$time

  lapply(seq_along(meanMat), function(i) {
    mm <- meanMat[[i]]
    parNames <- colnames(mm)
    mmList <- lapply(split(mm, row(mm)), function(x) {
      x <- as.list(x)
      x[[timeName]] <- time
      setNames(x, c(parNames, timeName))
    })
    res <- lapply(mmList, function(x) {force(x); do.call(curveFun, x)})
    res <- matrix(unlist(res, use.names = FALSE), nrow = length(res), byrow = TRUE)
    curveFit <- colMeans(res) # each column is a time point
    curveSD <- apply(res, 2, sd)
    structure(.Data = list(fit = curveFit, sd = curveSD,
                           curveMat = res, parMat = mm,
                           n = nrow(oP[[i]])),
              class = "bdCurveList")
  })
}


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
