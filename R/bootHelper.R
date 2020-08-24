
### Functions used in bdotsBoot ###

## bdotsBooter
# Takes subset dat with iter and corMat, returns
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

bdotsBooter <- function(dat, N.iter, corMat = NULL) {

  ## for now
  if (nrow(dat) > 2) stop("something weird in bdotsBooter")

  mm <- coef(dat)
  ## Can only be one or two (will this be always true?)
  if (!is.null(corMat)) {
    sig11 <- getVarMat(dat[1, ])
    sig22 <- getVarMat(dat[2, ])
    sig12 <- 0 * corMat %*% sqrt((diag(sig11)) %*% t(diag(sig22)))
    sig <- cbind(rbind(sig11, sig12), rbind(t(sig12), sig22))
    # ee <- eigen(sig)$values
    # ll <- min(abs(ee))/max(abs(ee))
    # sig <- Matrix::nearPD(sig, keepDiag = TRUE, eig.tol = ll*1e-2, maxit = 1e7)$mat
    pars <- mvtnorm::rmvnorm(N.iter, mean = c(t(mm)), sigma = sig)
  } else {
    sig <- getVarMat(dat)
    mm <- coef(dat)
    pars <- rmvnorm(N.iter, mm, sigma = sig)
  }
  colnames(pars) <- rep(colnames(mm), ncol(pars)/ncol(mm))
  pars
}


## getVarMat
# takes subset data with single observation
# returns covariance matrix of fit parameters
getVarMat <- function(dat) {
  if(nrow(dat) != 1) stop("only for single row of bdotsObj")
  dat$fit[[1]]$varBeta
}



## Extract coef from  bdotsObj
## uh, this doesn't adress null
# Ah, mother fucker, that's fitcode 6!
# this can potentially be made cleaner
#### Can't replace fit[[1]] since it's unnamed length 1 list. Could name it, I guess
## Needs to be made generic, address fitcode mentioned above.
# this needs to be renamed coef.bdObj
coef.bdotsObj <- function(dat) {
  #if (!inherits(dat, "bdotsObj")) stop('need bdotsObj')
  nnfit_v <- which(vapply(dat$fit, function(x) !is.null(x$fit), logical(1))) #dat$fitCode != 6 (change here and somewhere else I remember)
  if (!length(nnfit_v)) {
    warning("No models contain valid coefficients")
    # return(NULL)
  }
  mm <- matrix(NA, nrow = nrow(dat), ncol = length(cc <- coef(dat[nnfit_v[1], ]$fit[[1]])))
  colnames(mm) <- names(cc)
  for (i in seq_along(1:nrow(mm))) {
    if (dat[i, ]$fitCode != 6) mm[i, ] <- coef(dat[i, ]$fit[[1]])
  }
  mm
}


## curveBooter
# this function is dangerously long
# and potentially complicated
# need to revist for potential simplication.

## Returns length 2 nested list (outdated) ((super outdated))
# 1. curve1
#   i. curveMat
#   ii. parMat
# 2. curve2
#   i. curveMat
#   ii. parMat
## Make make sense to change inner/outer diff
# The 'n' for diff is associated with the t-statistic, not actual count
## This function is a bit of a bohemeth
curveBooter <- function(Obj, outerDiff, innerDiff = NULL, N.iter, curveFun) {

  if(!is.null(innerDiff)) {
    obj <- split(Obj, by = outerDiff, drop = TRUE)
    res <- lapply(obj, curveBooter, outerDiff = innerDiff,
                  N.iter = N.iter, curveFun = curveFun)
    res <- unlist(res, recursive = FALSE)

    idx <- grep("diff", names(res))
    if (length(idx) != 2) stop("something weird in curveBooter. Contact author")

    ## diff of diff (length one list)
    diffList <- Map(function(x, y) {
      Map(function(a, b) {
        a - b
      }, x, y)
    }, res[idx[1]], res[idx[2]])

    diffList <- diffList[[1]]

    ## snap, we can
    if (ip <- isPaired(obj)) {
      diffList$sd <- apply(diffList$curveMat, 2, sd)
      diffList$n <- nrow(obj[[1]]) - 1
    } else {
      diffList$sd <- nopairSD(res[idx])
      diffList$n <- sum(vapply(obj, nrow, numeric(1))) - 2
    }
    diffList$paired <- ip

    ## Maybe first do something to combine these?
    ## specfically, we should only return a single diff

    return(setNames(c(res, list(diffList)), c(names(res), "diff")))
  }
  oP <- split(Obj, by = outerDiff, drop = TRUE)
  if (ip <- isPaired(oP)) {
    cm <- lapply(oP, coef)
    corMat <- do.call("cor", setNames(cm, c("x", "y")))
  } else {
    corMat <- NULL
  }
  ## Should investigate how these are different
  # oh, on the split
  if (!is.null(corMat)) {
    outDiffL <- split(Obj, by = "Subject", drop = TRUE)
    bootPars <- lapply(outDiffL, bdotsBooter, N.iter, corMat)
    meanMat <- parMatSplit(Reduce(`+`,  bootPars)/length(bootPars))
  } else {
    outDiffL <- lapply(oP, split, by = "Subject")
    meanMat <- lapply(outDiffL, function(x) {
      bootPars <- lapply(x,  bdotsBooter, N.iter, corMat)
      meanMat <- Reduce(`+`,  bootPars)/length(bootPars)
    })
  }

  time <- attr(Obj, "time")
  curveList <- lapply(seq_along(meanMat), function(i) {
    mm <- meanMat[[i]]
    parNames <- colnames(mm)
    mmList <- lapply(split(mm, row(mm)), function(x) {
      x <- as.list(x)
      x$time <- time
      setNames(x, c(parNames, "time"))
    })
    ## Note, use this  to get mean and sd for the curves
    res <- lapply(mmList, function(x) {force(x); do.call(curveFun, x)})
    res <- matrix(unlist(res, use.names = FALSE), nrow = length(res), byrow = TRUE)
    curveFit <- colMeans(res)
    curveSD <- apply(res, 2, sd)
    list(fit = curveFit, sd = curveSD, curveMat = res, parMat = mm, n = nrow(oP[[i]]))
  })

  ## This is a bit sloppy
  nn <- names(curveList[[1]])
  diffList <- structure(vector("list", length = length(nn) + 1),
                        names = c(nn, "paired"))

  ## Could probably do the double map here like above
  # maybe see which is faster?
  tmp <- unlist(curveList, recursive = FALSE, use.names = FALSE)
  for (i in seq_along(nn)) {
    diffList[[i]] <- tmp[[i]] - tmp[[i + length(nn)]]
  }

  if (ip) {
    diffList$sd <- apply(diffList$curveMat, 2, sd)
    diffList$n <- nrow(oP[[1]]) - 1
  } else {
    diffList$sd <- nopairSD(curveList)
    diffList$n <- sum(vapply(oP, nrow, numeric(1))) - 2
  }
  diffList$paired <- ip


  ## Let's return all that above, and do the diff stuff here (wasted computation if not needed, I guess, but it's only matrix calc)
  setNames(c(curveList, list(diffList)), c(unique(Obj[[outerDiff]]), "diff"))
}




