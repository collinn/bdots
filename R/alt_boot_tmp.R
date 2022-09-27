

bdotsBoot2 <- function(formula,
                      bdObj,
                      Niter = 1000,
                      alpha = 0.05,
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
  
  x <- splitGroups[[1]]
  time <- attr(bdObj, "time") # right now, just use the union of times
  # time <- seq(min(time), max(time), length.out = 1000)
  groupDists <- lapply(splitGroups, getBootDist, b = Niter) # don't wanna call this Niter
  
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

## Create distribution for group
x <- splitGroups[[1]]
getBootDist <- function(x, b = 1000) {
  
  ## First, things I need from x
  curveFun <- makeCurveFun(x)
  time <- attr(x, "time")
  timeName <- attr(x, "call")$time
  parNames <- colnames(coef(x))
  
  ## Then we get the bootstrapped parameters as a list
  bsPars <- function(x) {
    idx <- sample(seq_len(nrow(x)), replace = TRUE)
    xn <- x[idx, ]
    xn$splitvar <- seq_len(nrow(x))
    xns <- split(xn, by = "splitvar")
    xpar <- vapply(xns, function(z) {
      rmvnorm(1, coef(z), vcov(z$fit[[1]]))
    }, numeric(4))
    rowMeans(xpar)
  }
  mm <- replicate(b, bsPars(x), simplify = TRUE) |> t()
  parList <- lapply(split(mm, row(mm)), function(y) {
    y <- as.list(y)
    y[[timeName]] <- time
    setNames(y, c(parNames, timeName))
  })
  
  ## Finally, function time
  res <- lapply(parList, function(y) {force(y); do.call(curveFun, y)})
  res <- matrix(unlist(res, use.names = FALSE), nrow = length(res), byrow = TRUE)
}
