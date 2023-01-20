

## Permutation test for bdots
# not available for difference of difference
# right now going to be a lot of piecewise logic as i figure out what im doing
#' Permutation testing for bboot
#' 
#' Provides alternative for controlling FWER
#' 
#' @param x this is splitGroups (each item of list bdotsObj)
#' @param prs from parser, indicates inner/outer groups
#' @param alpha alpha for fwer
#' @param P number of permutations
permTest <- function(x, prs, alpha, P) {
  ip <- isPaired(x)
  if (ip) warning("paired permutation testing current not available. Running unpaired")
  ip <- FALSE # paired testing currently not implemented
  dod <- !is.null(prs$innerDiff)
  if (dod) stop("Permutation not available yet for difference of differences")
  
  if (!ip & !dod) {
    pgroups <- prs$outerDiff
    x <- bdots:::rbindlist.bdObjList(x)
    n <- nrow(x)
    permmat <- replicate(P, sample(seq_len(n), n))
    
    ## Get the t stats for observed
    tvec <- getT(x, seq_len(n), group = pgroups, whole = TRUE)
    
    ## Get max tstat across permutations (null distribution)
    tnull <- apply(permmat, 2, function(y) {
      getT(x, y, group = pgroups)
    })
    qq <- quantile(tnull, probs = 1-alpha/2)
    sigIdx <- tvec > qq
    return(list(obst = tvec, nullt = tnull, sigIdx = sigIdx))
  }
}

#' Generate t-distribution for permutation
#' 
#' Creates vector of t-statistics for given permutation returning
#' either entire vector or max values
#' 
#' @param x bdObj
#' @param idx permutation to use
#' @param group group that we are permuting against
#' @param whole return vector of T stats or just the max
getT <- function(x, idx, group, whole = FALSE) {
  x[[group]] <- x[[group]][idx]
  
  ## Stuff I need to  call the function (gross)
  timeName <- attr(x, "call")$time
  TIME <- attributes(x)$time
  ff <- makeCurveFun(x)
  
  fit_s <- split(x, by = group)
  mvl <- lapply(fit_s, function(y) {
    cc <- coef(y)
    cl <- apply(cc, 1, function(z) {
      z <- as.list(z)
      z[[timeName]] <- TIME
      do.call(ff, z)
    })
    mm <- rowMeans(cl)
    vv <- apply(cl, 1, var)
    vvn <- vv/nrow(cc)
    list(mean = mm, nvar = vvn)
  })
  
  ## Guess I don't need a function for this
  x <- mvl[[1]]; y <- mvl[[2]]
  xm <- x$mean; xv <- x$nvar
  ym <- y$mean; yv <- y$nvar
  Tt <- abs(xm-ym) / sqrt(yv + xv)
  
  ifelse(whole, return(Tt), return(max(Tt)))
}
