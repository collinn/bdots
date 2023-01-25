

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
permTest <- function(x, prs, alpha, P, cores = detectCores()-1L) {
  ip <- isPaired(x)
  # if (ip) warning("paired permutation testing current not available. Running unpaired")
  #ip <- FALSE # paired testing currently not implemented
  dod <- !is.null(prs$innerDiff)
  if (dod) stop("Permutation not available yet for difference of differences")

  ## Need these for all
  pgroups <- prs$outerDiff
  x <- bdots:::rbindlist.bdObjList(x)
  n <- nrow(x)

  ## Get the t stats for observed
  tvec <- getT(x, seq_len(n), group = pgroups, whole = TRUE)

  ## I guess let's do parallel its faster
  if (Sys.info()['sysname'] == "Darwin") {
    cl <- makePSOCKcluster(cores, setup_strategy = "sequential")
  } else {
    cl <- makePSOCKcluster(cores)
  }
  invisible(clusterEvalQ(cl, {library(bdots)}))

  ## Now compute null distribution
  if (!ip) {
    permmat <- replicate(P, sample(seq_len(n), n))

    ## Get max tstat across permutations (null distribution)
    tnull <- parApply(cl, permmat, 2, function(y) {
      getT(x, y, group = pgroups)
    })
    qq <- quantile(tnull, probs = 1-alpha/2)
    sigIdx <- tvec > qq
  } else  if (ip) {
    ## First half of permmat
    permmat <- replicate(P, sample(c(TRUE, FALSE), n/2, replace = TRUE))
    permmat <- rbind(permmat, !permmat)
    tnull <- parApply(cl, permmat, 2, function(y) {
      bidx <- bool2idx(y)
      getT(x, bidx, group = pgroups)
    })
    qq <- quantile(tnull, probs = 1-alpha/2)
    sigIdx <- tvec > qq
  }

  stopCluster(cl)

  return(list(obst = tvec, nullt = tnull, sigIdx = sigIdx))

}


#' Need to turn  T/F into indices and I think I have an idea how
#' @param b vector of boolean
bool2idx <- function(b) {
  n <- length(b)
  nidx <- seq_len(n)
  idx <- c(nidx[b], nidx[!b])
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

#' Get null distribution of permutation
#'
#' Get null distribution of permutation
#'
#' @param x bdObj
#' @param P number of permutations
getTDist <- function(x, P = 1000) {
  n <- nrow(x)
  permmat <- replicate(P, sample(seq_len(n), n))
  tvec <- getT(x, seq_len(n), whole = TRUE)
  tnull <- apply(permmat, 2, function(y) {
    getT(x, y)
  })
}
