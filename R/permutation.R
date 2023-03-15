

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
  # these all implicitly assumed sorted by subject so maybe should do that in bdotsFit
  ip <- isPaired(x)
  dod <- !is.null(prs$innerDiff)
  if (dod) {
    msg <- "Permutation not available for difference of differences."
    msg <- paste(msg, "Please contact the package author if this is something",
    "you're interested in having implemented.")
    stop(msg)
  }

  ## I guess let's do parallel its faster
  if (Sys.info()['sysname'] == "Darwin") {
    cl <- makePSOCKcluster(cores, setup_strategy = "sequential")
  } else {
    cl <- makePSOCKcluster(cores)
  }
  invisible(clusterEvalQ(cl, {library(bdots)}))

  ## Need these for all
  pgroups <- prs$outerDiff
  x <- bdots:::rbindlist.bdObjList(x)
  n <- nrow(x)

  ## Get the t stats for observed
  tvec <- getT(x, seq_len(n), group = pgroups, whole = TRUE, addVar = FALSE, ip)

  ## Start with case no DOD
  if (!dod) {
    ## Now compute null distribution
    if (!ip) {
      permmat <- replicate(P, sample(seq_len(n), n))

      ## Get max tstat across permutations (null distribution)
      tnull <- parApply(cl, permmat, 2, function(y) {
        getT(x, y, group = pgroups)
      })
      qq <- quantile(tnull, probs = 1-alpha)
      sigIdx <- tvec > qq
    } else if (ip) {
      ## First half of permmat
      permmat <- replicate(P, sample(c(TRUE, FALSE), n/2, replace = TRUE))
      permmat <- rbind(permmat, !permmat)
      tnull <- parApply(cl, permmat, 2, function(y) {
        bidx <- bool2idx(y)
        getT(x, bidx, group = pgroups, whole = FALSE, addVar = FALSE, ip = ip)
      })
      qq <- quantile(tnull, probs = 1-alpha)
      sigIdx <- tvec > qq
    }
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
#' @param ip boolean indicating if paired
getT <- function(x, idx, group, whole = FALSE, addVar = TRUE, ip = FALSE) {
  x[[group]] <- x[[group]][idx]

  ## Stuff I need to  call the function (gross)
  timeName <- attr(x, "call")$time
  TIME <- attributes(x)$time
  ff <- makeCurveFun(x)

  fit_s <- split(x, by = group)

  ## Do i sample with the added variation?
  if (addVar == FALSE) {
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
      list(mean = mm, nvar = vvn, nn = nrow(cc), curveList = cl)
    })
  } else {
    mvl <- lapply(fit_s, function(y) {
      ## Weird that this becomes a list
      cc <- apply(y, 1, function(z) {
        rmvnorm(1, coef(z$fit), vcov(z$fit))
      }) |> t()
      #cc <- coef(y)
      cl <- apply(cc, 1, function(z) {
        z <- as.list(z)
        z[[timeName]] <- TIME
        do.call(ff, z)
      })
      mm <- rowMeans(cl)
      vv <- apply(cl, 1, var)
      vvn <- vv/nrow(cc)
      list(mean = mm, nvar = vvn, nn = nrow(cc), curveList = cl)
    })
  }


  ## Guess I don't need a function for this
  x <- mvl[[1]]; y <- mvl[[2]]
  xm <- x$mean; xv <- x$nvar
  ym <- y$mean; yv <- y$nvar

  ## If paired, this will be same as y$nn
  clx <- x$curveList; cly <- y$curveList

  if (!ip) {
    Tt <- abs(xm-ym) / sqrt(yv + xv)
  } else {
    dif <- clx - cly
    dd <- rowMeans(dif)
    vv <- apply(dif, 1, var)
    vv <- vv / x$nn
    Tt <- abs(dd) / sqrt(vv)
  }

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


#### I'm leaving this all commented out
# Haven't figured out degrees of freedom and various other tidbits for t statistic
# and not mission critical to dissertation

# ## Step 1. What kind of pair we got
# ps <- isDODpaired(x, prs)
# P <- 10
#
# ## Step 2. Combine into paired up inner groups
# xx <- rbindlist(x)
# xx <- split(xx, by = prs[["outerDiff"]])
# group <- prs[['innerDiff']]
#
# x <- xx[[1]]
# idx <- sample(1:20)
# getInnerPermuteMean <- function(x, idx, group, ip) {
#   newvec <-  x[[group]][idx]
#   set(x, j = group, value = newvec)
#
#   timeName <- attr(x, "call")$time
#   TIME <- attributes(x)$time
#   ff <- makeCurveFun(x)
#
#   fit_s <- split(x, by = group)
#   mvl <- lapply(fit_s, function(y) {
#     cc <- coef(y)
#     cl <- apply(cc, 1, function(z) {
#       z <- as.list(z)
#       z[[timeName]] <- TIME
#       do.call(ff, z)
#     })
#     #rowMeans(cl)
#   })
#
#   if (ip) {
#     dm <- Reduce(`-`, mvl)
#     ## nrow(x)/2 since paired
#     n <- nrow(x)/2
#     dv <- apply(dm, 1, var) / n
#     ds <- sqrt(dv)
#   } else {
#     dm <- lapply(mvl, rowMeans)
#     dm <- Reduce(`-`, dm)
#   }
#
#   dm <- Reduce(`-`, mvl) # this won't work if they are different sizes so rowmeans first
#   dm
# }
#
# ## Step 3. Get indices for permutation, idxa and idxb
# # a. if paired in group find half permute then !flip (since boolean)
# # b. if paired in outer make idxa == idxb
# # c. if not paired at all generate four random idx based on n subjects
# ## I will do this by creating a function that gives inner diff index
# # if (b) just use once, it will return conditional on (a) or (b)
# innerIndex <- function(P, n, pair) {
#   if (pair) {
#     permmat <- replicate(P, sample(c(TRUE, FALSE), n/2, replace = TRUE))
#     permmat <- rbind(permmat, !permmat)
#     permmat <- apply(permmat, 2, bool2idx)
#   } else {
#     permmat <- replicate(P, sample(seq_len(n)))
#   }
#   permmat
# }
#
# x <- splitGroups
# ps <- isDODpaired(x, prs)
# xx <- rbindlist(x)
# xx <- split(xx, by = prs[["outerDiff"]])
# nv <- vapply(xx, nrow, numeric(1))
# P <- 20
# ip <- ps['ip']
# op <- ps['op']
#
# ## if both paired
# if (sum(ps) == 2) {
#   permA <- innerIndex(P, nv[1], pair = TRUE)
#   permB <- permA
# } else if (sum(ps) == 0) {
#   ## Neither paired
#   permA <- innerIndex(P, nv[1], pair = FALSE)
#   permB <- innerIndex(P, nv[2], pair = FALSE)
# } else if (ip) {
#   # both paired but outer different
#   permA <- innerIndex(P, nv[1], pair = TRUE)
#   permB <- innerIndex(P, nv[2], pair = TRUE)
# } else if (op) {
#   permA <- innerIndex(P, nv[1], pair = FALSE)
#   permB <- permA
# }
#
# ## Step 4. Compute the inner mean vectors based on this
# # I guess do this in two steps at first
# clusterExport(cl, c("getInnerPermuteMV", "group", "xx", "ip"))
# clusterEvalQ(cl, library(bdots))
# clusterEvalQ(cl, devtools::load_all("~/packages/bdots"))
#
# ## Could also do this with Map but does not lend itself to parallelization like this does
# diffmatA <- parApply(cl, permA, 2, function(y) {
#   getInnerPermuteMV(xx[[1]], y, group, ip)
# })
# diffmatB <- parApply(cl, permB, 2, function(y) {
#   getInnerPermuteMV(xx[[2]], y, group, ip)
# })
# # I now need a t statisitic for each row which means i also need some sort of sd
# diffmat <- diffmatA - diffmatB
#
