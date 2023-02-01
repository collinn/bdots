## curveBooter
# this function is dangerously long
# and potentially complicated
# need to revist for potential simplication.

# _sm designation means its for single mean assumptions

curveBooter_sm <- function(Obj, outerDiff, innerDiff = NULL, N.iter, curveFun) {

  if (!is.null(innerDiff)) {
    obj <- split.bdotsObj(Obj, by = outerDiff, drop = TRUE)
    res <- lapply(obj, curveBooter, outerDiff = innerDiff,
                  N.iter = N.iter, curveFun = curveFun)

    diffList <- makeOuterDiffList_sm(res, obj)

    return(structure(.Data = setNames(c(res, list(diffList)),
                               c(names(res), "diff")),
                     class = c("outerGroupCurveList","groupCuveList")))
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
  curveList <- makeCurveList_sm(meanMat, curveFun, oP)

  ## Class bdDiffList
  diffList <- makeInnerDiffList_sm(curveList, oP)

  structure(.Data = setNames(c(curveList, list(diffList)),
                             c(unique(Obj[[outerDiff]]), "diff")),
            class = c("innerGroupCurveList", "groupCuveList"))
}

###------------------------------------------------

## This function is responsible for taking the mean parameter matrix
# from each group over N.iter iterations. Along with time and the function
# specifying the curve to be fit, this returns a length 2 list, one for each
# of the groups being fit

# takes list of meanMatrix for each group (from bdotsBooter),
# and a numeric vec with timeName attributes (from original call)
makeCurveList_sm <- function(meanMat, curveFun, oP) {
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
makeInnerDiffList_sm <- function(curveList, oP) {
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
makeOuterDiffList_sm <- function(res, obj) {
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

## bdotsBooter
# Takes subset dat with iter and corMat, returns
# More specifically, it takes a bdotsObj
# N.iter x npars matrix of random draws

# Notes
## bdotsBooter
# This will have at most two curves, for when things need
# to be bivariate normal. It will also have an argument for
# correlation coefficient. It does not need to know diffGroup or fitGroup
# it only needs to know correlated or not
## The bivariate normal matrix is going to always be problematic
# because var for base1, base2, ht are nearly 0 (widely different scales than mu, sig1, sig2)
# just look at kappa(sig11)
# for now, I will leave it. I  will try to come up with my own solution before
# giving this to jake
## very much same issue with logistic (cross is HUGE)
#' @import mvtnorm
bdotsBooter <- function(bdo, N.iter, corMat = NULL) {

  ## for now
  if (nrow(bdo) > 2) stop("something weird in bdotsBooter")

  mm <- coef(bdo)
  ## Can only be one or two (will this be always true?)
  if (!is.null(corMat)) {
    sig11 <- getVarMat(bdo[1, ])
    sig22 <- getVarMat(bdo[2, ])
    sig12 <- 0 * corMat %*% sqrt((diag(sig11)) %*% t(diag(sig22)))
    sig <- cbind(rbind(sig11, sig12), rbind(t(sig12), sig22))
    # ee <- eigen(sig)$values
    # ll <- min(abs(ee))/max(abs(ee))
    # sig <- Matrix::nearPD(sig, keepDiag = TRUE, eig.tol = ll*1e-2, maxit = 1e7)$mat
    pars <- mvtnorm::rmvnorm(N.iter, mean = c(t(mm)), sigma = sig)
  } else {
    sig <- getVarMat(bdo)
    mm <- coef(bdo)
    pars <- mvtnorm::rmvnorm(N.iter, mm, sigma = sig)
  }
  colnames(pars) <- rep(colnames(mm), ncol(pars)/ncol(mm))
  pars
}
