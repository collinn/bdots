#' Create bootstrapped distribution of groups
#'
#' Function to create bootstrapped distribution of groups depending
#' on paired status and difference of difference
#'
#' @param x list of bdObj
#' @param prs list from boot parser
#' @param b number of bootstraps
#' @param cores cores
createGroupDists <- function(x, prs, b, cores) {

  ## Now tend to parallel needs
  if (Sys.info()['sysname'] == "Darwin") {
    cl <- makePSOCKcluster(cores, setup_strategy = "sequential")
  } else {
    cl <- makePSOCKcluster(cores)
  }
  invisible(clusterEvalQ(cl, {library(bdots)}))

  ## Is this a difference of difference?
  dod <- !is.null(prs$innerDiff)

  if (!dod) {
    ## This means only length 2 so checking if paired is simple
    ip <- isPaired(x)
    if (!ip) {
      groupDists <- parLapply(cl, x, getBootDistUnpaired, b) # <- the old way, nice
      #stopCluster(cl)
      #return(groupDists) # done
    } else {
      # if it is paired
      groupDists <- getBootDistPaired(cl, x, b)
      #stopCluster(cl)
      #return(groupDists)
    }
  } else { # if difference of difference
    ## Find out if inner/outer paired
    ps <- isDODpaired(x, prs)

    ## Based on paired status
    if (sum(ps) == 2) {
      ## both paired
      groupDists <- getBootDistPaired(cl, x, b)
    } else if (sum(ps) == 0) {
      ## Neither paired
      groupDists <- parLapply(cl, x, getBootDistUnpaired, b)
    } else {
      ## This is the case if either inner or outer not but other is
      nx <- strsplit(names(x), split = "\\.")

      if (ps['ip']) {
        # paired inner group will share the same outer group in the split object
        ogn <- prs$subargs[[2]] # I hate this but this is values for outerDiff which also doesn't even make sense
        idx <- vapply(nx, function(y) ogn[1] %in% y, logical(1))
      } else if (ps['op']) {
        ign <- prs$subargs[[1]]
        idx <- vapply(nx, function(y) ign[1] %in% y, logical(1))
      }
      gd1 <- getBootDistPaired(cl, x[idx], b)
      gd2 <- getBootDistPaired(cl, x[!idx], b)
      groupDists <- c(gd1, gd2)
    }
  }
  stopCluster(cl)
  return(groupDists)
}

#' Create distribution for paired groups
#'
#' @param cl cluster for parallel
#' @param x A list of bdObj
#' @param b Number of bootstraps
getBootDistPaired <- function(cl, x, b) {

  curveFun <- makeCurveFun(x[[1]])
  time <- attr(x[[1]], "time")
  timeName <- attr(x[[1]], "call")$time
  parNames <- colnames(coef(x[[1]]))

  ## Stolen from other function
  bsPars <- function(x, idx) {
    #idx <- sample(seq_len(nrow(x)), replace = TRUE)
    xn <- x[idx, ]
    pp <- ncol(coef(xn))
    xn$splitvar <- seq_len(nrow(x))
    xns <- split(xn, by = "splitvar")
    xpar <- vapply(xns, function(z) {
      rmvnorm(1, coef(z), vcov(z$fit[[1]]))
    }, numeric(pp))
    rowMeans(xpar)
  }

  # getIdxThenPars <- function(y) {
  #   idx <- sample(seq_len(nrow(y[[1]])), replace = TRUE)
  #   pars1 <- t(bsPars(y[[1]], idx))
  #   pars2 <- t(bsPars(y[[2]], idx))
  #   list(pars1 = pars1, pars2 = pars2)
  # }

  getIdxThenPars <- function(y) {
    idx <- sample(seq_len(nrow(y[[1]])), replace = TRUE)
    pars <- lapply(y, bsPars, idx)
  }

  #clusterEvalQ(cl, library(mvtnorm))
  #clusterExport(cl, varlist = c("bsPars", "getIdxThenPars", "x"))
  clusterExport(cl, varlist = ls(), envir = environment())

  ## Ok get coefs from each group
  mm <- parSapply(cl, integer(b), function(...) getIdxThenPars(x),
                  USE.NAMES = FALSE)

  ## Really easier to just do this in long steps even if gross
  turnToDist <- function(m) {
    parList <- lapply(m, function(y) {
      y <- as.list(y)
      y[[timeName]] <- time
      setNames(y, c(parNames, timeName))
    })
    res <- lapply(parList, function(y) {force(y); do.call(curveFun, y)})
    res <- matrix(unlist(res, use.names = FALSE), nrow = length(res), byrow = TRUE)

    ## Ok, let's end by making this a bdCurveList to match other bdots
    fit <- colMeans(res)
    sd <- apply(res, 2, sd)

    pp <- matrix(unlist(m), nrow = length(m), byrow = TRUE)

    structure(.Data = list(fit = fit, sd = sd,
                           curveMat = res, parMat = pp,
                           n = length(m)),
              class = "bdCurveList")
  }

  d <- lapply(split(mm, row(mm)), turnToDist)
  res <- setNames(d, names(x))

  # d1 <- turnToDist(mm[1, ])
  # d2 <- turnToDist(mm[2, ])
  # res <- setNames(list(d1, d2), names(x))
  return(res)

}


#' Create Distribution for Groups
#'
#' @param x A subset object from bdObj
#' @param b Number of bootstraps
getBootDistUnpaired <- function(x, b = 1000) {

  ## First, things I need from x
  curveFun <- makeCurveFun(x)
  time <- attr(x, "time")
  timeName <- attr(x, "call")$time
  parNames <- colnames(coef(x))

  ## Then we get the bootstrapped parameters as a list
  bsPars <- function(x) {
    idx <- sample(seq_len(nrow(x)), replace = TRUE)
    xn <- x[idx, ]
    pp <- ncol(coef(xn))
    xn$splitvar <- seq_len(nrow(x))
    xns <- split(xn, by = "splitvar")
    xpar <- vapply(xns, function(z) {
      rmvnorm(1, coef(z), vcov(z$fit[[1]]))
    }, numeric(pp))
    rowMeans(xpar)
  }
  mm <- replicate(b, bsPars(x), simplify = TRUE)
  mm <- t(mm)
  parList <- lapply(split(mm, row(mm)), function(y) {
    y <- as.list(y)
    y[[timeName]] <- time
    setNames(y, c(parNames, timeName))
  })

  ## Finally, function time
  res <- lapply(parList, function(y) {force(y); do.call(curveFun, y)})
  res <- matrix(unlist(res, use.names = FALSE), nrow = length(res), byrow = TRUE)

  ## Ok, let's end by making this a bdCurveList to match other bdots
  fit <- colMeans(res)
  sd <- apply(res, 2, sd)

  structure(.Data = list(fit = fit, sd = sd,
                         curveMat = res, parMat = mm,
                         n = nrow(x)),
            class = "bdCurveList")
}


## Need a function that determines if inner or
# outer or both or neither of 4 groups is paired
#' Difference of difference pairing
#'
#' Function to determine if list of 4 groups indicating difference of
#' difference is paired or not
#'
#' @param x list of bdObjs
#' @param prs parsed boot objects
#'
#' @description Took the easy, inefficient way for now
isDODpaired <- function(x, prs) {
  if (length(x) != 4) stop("something wrong in isDODpaired")
  id <- prs[["innerDiff"]]
  od <- prs[["outerDiff"]]

  xx <- rbindlist(x)

  ## Check outer (this is single test, all against all)
  xo <- split(xx, by = od)
  odPaired <- isPaired(xo)

  ## To check inner paired, need to look at both inner groups
  ign <- xx[[id]] |> unique()
  nx <- strsplit(names(x), split = "\\.")
  idx <- vapply(nx, function(y) ign[1] %in% y, logical(1))

  ip1 <- isPaired(x[idx])
  ip2 <- isPaired(x[!idx])
  idPaired <- ip1 & ip2

  return(c(ip = idPaired, op = odPaired))
}
