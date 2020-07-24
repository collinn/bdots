
outerDiff <- "Group"; innerDiff <- NULL; bdObj2 <- bdObj[TrialType == "M", ]
outerDiff <- "TrialType"; innerDiff <- NULL; bdObj2 <- bdObj[Group == "LI", ]


## Start with simple case. innerDiff false (so we are doing TD_W - LI_W or TD_W - TD_M)
bdObj2 <- bdObj[TrialType == "M", ]
#  bdObj2 <- bdObj[Group == "LI", ] # need to try  both
outDiffL <- split(bdObj[TrialType == "M", ], by = outerDiff, drop = TRUE)
#outDiffL <- split(bdObj[Group == "LI", ], by = "Group", drop = TRUE)

## BELOW
# assume bdObj and bdObj2 are already subset correctly

if(!is.null(innerDiff)) {
  outDiffL <- split(bdObj, by = outerDiff, drop = TRUE)
  # stuff
} else { #is.null(innerDiff)

  ## So it appears that for is.null(innerDiff)
  # 1.) need to split bdObj by outer group to see if correlated
  # 2.) need to split bdObj by subject for fits
  # 3.) This is in contrast to !is.null(innerDiff)
  # 4.) I really just need a better way to determine if correlated
  # 5.) Actually, I think it's fine - outerDiff will ALWAYS be length 2. My tests here are wrong

  ## Determine if paired
  ## Could make isPaired more robust?
  outDiffL_P <- split(bdObj2, by = outerDiff, drop = TRUE)
  if (isPaired(outDiffL_P)) {
      cm <- lapply(outDiffL_P, coef.bdots)
      corMat <- do.call(cor, setNames(cm, c("x", "y")))
  } else {
      corMat <- NULL
  }

  ## Dealing with paired case
  ## curveList is length 2
  if(!is.null(corMat)) {
    ## Now make nested list of subjects
    # If !is.null(corMat), then there will be two subjects in each
    # component of this list. Otherwise, only 1
    outDiffL <- split(bdObj2, by = "Subject", drop = TRUE)

    ## This has to be N.iter x (pars * x) since corMat
    bootPars <- lapply(outDiffL, bdotsBooter, N.iter, corMat)
    meanMat <- Reduce(`+`,  bootPars)/length(bootPars)

    curveList <- lapply(parMatSplit(meanMat), function(mm) {
      parNames <- colnames(mm)
      mmList <- lapply(split(mm, row(mm)), function(x) {
        x <- as.list(x)
        x$time <- time
        setNames(x, c(parNames, "time"))
      })
      res <- lapply(mmList, function(x) {force(x); do.call(curveFun, x)})
      res <- matrix(unlist(res, use.names = FALSE), nrow = length(res), byrow = TRUE)
    })
  } else {
    ## Here, we have list of lists
    outDiffL <- split(bdObj2, by = c(outerDiff), drop = FALSE)
    outDiffL <- lapply(outDiffL, split, by = "Subject")
    meanMat <- lapply(outDiffL, function(x) {
      bootPars <- lapply(x,  bdotsBooter, N.iter, corMat)
      meanMat <- Reduce(`+`,  bootPars)/length(bootPars)
    })
    curveList <- lapply(meanMat, function(mm) {
      parNames <- colnames(mm)
      mmList <- lapply(split(mm, row(mm)), function(x) {
        x <- as.list(x)
        x$time <- time
        setNames(x, c(parNames, "time"))
      })
      res <- lapply(mmList, function(x) {force(x); do.call(curveFun, x)})
      res <- matrix(unlist(res, use.names = FALSE), nrow = length(res), byrow = TRUE)
    })
  }

}






