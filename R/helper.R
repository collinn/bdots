
## Subset a bdotsBootObj based on group
#' Subset a nested group bdotsBoot objects
#'
#' @param x An object returned from \code{bdotsBoot}
#' @param group A group to subset. Must be an outer group
#' @param adjustAlpha currently not used. Will give option to recompute adjusted alpha
#' @param ... Not used
#'
#' @details This function is used to subset a bdotsBootObject that was fit to compute
#' the difference of differences. This allows the user to subset out the outer group
#' in the comparison for plotting and investigation
#'
#' @export
subset.bdotsBootObj <- function(x, group, adjustAlpha = NULL, ...) {
  bdBootObj <- x #  need to just rename and not be lazy here
  if (!bdBootObj$dod)
    stop("No inner group to subset")

  if (!(group %in% names(bdBootObj$curveList)))
    stop("Invalid group for subset")

  ## Take only old group
  bdBootObj$curveList <- bdBootObj$curveList[[group]]
  bdBootObj$diffs <- setNames(bdBootObj$diffs['innerDiff'], 'outerDiff')
  bdBootObj$dod <- FALSE
  bdBootObj$curveGroups <- bdBootObj$curveGroups[bdBootObj$diffs]

  ## See if recomputed alpha, for now, w/e. Also, this feels gross
  bdBootObj$sigTime <- NULL
  bdBootObj$adj.alpha <- NULL
  bdBootObj$rho <- NULL
  bdBootObj$adj.pval <- NULL
  bdBootObj$paired <- bdBootObj$curveList[['diff']][['paired']]
  bdBootObj
}

## getVarMat
# takes subset data with single observation
# returns covariance matrix of fit parameters
getVarMat <- function(dat) {
  if(nrow(dat) != 1) stop("only for single row of bdotsObj")
  dat$fit[[1]]$varBeta
}

#' Extract bdotsFit Moedel Coefficients
#'
#' Returns coefficient matrix for bdotsFit object
#'
#' @param object A bdotsObj
#' @param ... not used
# @importFrom stats coef
#' @export
coef.bdotsObj <- function(object, ...) {
  #if (!inherits(dat, "bdotsObj")) stop('need bdotsObj')
  nnfit_v <- which(vapply(object$fit, function(x) !is.null(x$fit), logical(1))) #dat$fitCode != 6 (change here and somewhere else I remember)
  if (!length(nnfit_v)) {
    warning("No models contain valid coefficients")
    # return(NULL)
  }
  mm <- matrix(NA, nrow = nrow(object), ncol = length(cc <- coef(object[nnfit_v[1], ]$fit[[1]])))
  colnames(mm) <- names(cc)
  for (i in seq_along(1:nrow(mm))) {
    if (object[i, ]$fitCode != 6) mm[i, ] <- coef(object[i, ]$fit[[1]])
  }
  mm
}



## Make split retain bdotsObj class
# Need to also split data attribute
#' Split object of class bdotsObj
#'
#' Analogous to other splitting functions, but retains necessary attributes
#' across the split object. As of now, it can only be unsplit with bdots::rbindlist
#'
#' @param x Object of class bdotsObj
#' @param f For consistency with generic, but is not used
#' @param drop logical. Default FALSE will not drop empty list elements caused
#' by factor levels not referred by that factor. Analagous to data.table::split
#' @param by Character vector of column names on which to split. Usually will
#' be Subject or one of the fitted groups
#' @param ... not used
#'
#' @import data.table
#' @export
split.bdotsObj <- function(x, f, drop = FALSE, by,...) {
  oldAttr <- attributes(x)
  class(x) <- c("data.table", "data.frame")
  res <- lapply(split(x, by = by, drop = drop, ...), function(y) {
    attributes(y) <- oldAttr
    y
  })
  structure(.Data = res, class = c("bdObjList"))
}
# split.bdotsObj <- function(bdo, by, ...) {
#   oldAttr <- attributes(bdo)
#   res <- lapply(data.table:::split.data.table(bdo, by = by, ...), function(x) {
#     attributes(x) <- oldAttr
#     x
#   })
#   structure(.Data = res, class = c("bdObjList"))
# }

## Don't export for now because fuck S3 generic matching
# @export
rbindlist <- function(x, ...) {
  UseMethod("rbindlist")
}

#' @importFrom data.table rbindlist
rbindlist.default <- function(x, ...) {
  data.table::rbindlist(x, ...)
}

#' rbindlist for bdotsObjects
#'
#' Similar to data.table::rbindlist, but preserves botsObjects attributes
#'
#' @param bdo bdotsObject
#' @param ... for compatability with data.table
#'
#'
#' @export
rbindlist.bdObjList <- function(x, ...) {
  oldAttr <- attributes(x[[1]])
  class(x) <- "list"
  x <- rbindlist(x, ...)
  attributes(x) <- oldAttr
  x
}


## Used in bdotsBoot
isPaired <- function(l) { # this only works for lists of bdObj
  subject <- attr(l[[1]], "call")[['subject']]
  reduceEq <- function(x, y) if (identical(x[[subject]], y[[subject]])) x else FALSE
  tmp <- Reduce(reduceEq, l) # how did I do this?
  if(identical(tmp, FALSE)) FALSE else TRUE
}

##  parMatSplit
# takes an N x (2p) matrix and returns list of (N x p), (N x p)

parMatSplit <- function(x) {
  n <- ncol(x)/2
  m1 <- x[, 1:n]
  m2 <- x[, (n+1):(2*n)]
  list(m1, m2)
}
###########################################
#####
## https://en.wikipedia.org/wiki/Student%27s_t-test
#####
### There are a bunch that are possible

## Get standard deviation for t statistic
# in unpaired test
## this t-stat is not even close to correct
# nopairSD <- function(l) {
#   if(length(l) != 2) stop("contact author with 123")
#   s <- lapply(l, function(x) {
#     vv <- x[['sd']]^2
#     n <- x[['n']]
#     (n-1) * vv
#   })
#   s <- Reduce(`+`, s) * sqrt(1/l[[1]]$n + 1/l[[2]]$n)
# }
## I think this is correct, but gives strange fits sometimes
nopairSD <- function(l) {
  if(length(l) != 2) stop("contact author with 1239")
  s <- lapply(l, function(x) {
    vv <- x[['sd']]^2
    n <- x[['n']]
    (n-1) * n * vv
  })
  s <- Reduce(`+`, s) / (l[[1]]$n + l[[2]]$n - 2) * (1/l[[1]]$n + 1/l[[2]]$n)
  s <- sqrt(s)
}

## this ONLY returns the sd for the t stat, not the multiplier
## For now, making it the entire denominator. we will check plot and see if it makes sense
nopairSD2 <- function(l) {
  sd_ratio <- Reduce(function(x, y) {x$sd/y$sd}, l)
  ## Proportion of sd ratio within bounds should be some val, lets say 0.5
  var_sim <- mean(sd_ratio > 1/2 & sd_ratio < 2) > 0.5
  if (var_sim) {
    s <- lapply(l, function(x) {
      vv <- x[['sd']]^2
      n <- x[['n']]
      (n - 1) * vv
    })
    s <- Reduce(`+`, s) / (l[[1]]$n + l[[2]]$n - 2) * (l[[1]]$n + l[[2]]$n)
    s <- sqrt(s)
    dof <- l[[1]]$n + l[[2]]$n - 2
  } else {
    s <- lapply(l, function(x) {
      vv <- x[['sd']]^2
      n <- x[['n']]
      vv/n
    })
    dof_denom <- Reduce(`+`, Map(function(x, y) {x^2 / (y[['n']] - 1)}, x = s, y = l))
    s <- Reduce(`+`, s)
    dof <- s^2 / dof_denom
    s <- sqrt(s)
  }
  return(list(sd = s, dof = dof))
}

###################################
## Stolen from purrrrrrrr
vec_depth <- function(x) {
  if (is.null(x)) {
    0L
  } else if (is.atomic(x)) {
    1L
  } else if (is.list(x)) {
    depths <- vapply(x, vec_depth, numeric(1))
    1L + max(depths, 0L)
  } else {
    stop("'x' must be a vector")
  }
}


## Can't be entirely sure how this generalizes yet
# since my cases will always be 2 x (2 x (2 x ...( 2 x 2)))
## Motherfucker, this does work
# x <- "x"; y <- "y"; z <- "z"
# ll <- list()
# ll[[1]] <- list(x, y)
# ll[[2]] <- list(y, list(z, z, x))
# ll[[3]] <- list(z, list(x, x), list(y, z, list(x , x)))
## Also a trillion times faster with rlang::squash. Sad
## Make sure it definitely returns a list
unzipList <- function(l) {
  cc <- vapply(l, is.list, logical(1L))
  res <- lapply(l[cc], unzipList)
  res <- c(l[!cc], unlist(res, recursive = FALSE))
  if (!is.list(res)) res <- list(res)
  res
}


## Comes up quite a bit, we want to associate
# group names with group values for summaries and plots
# function takes a list, first element names, second is values
# may or may not be split
# this list item should be named as such
# l <- list(groups = c("Group", "TrialType"),
#           vals = c("LI.M", "LI.W", "TD.M", "TD.W"))
# l <- list(groups = c("Group"), vals = c("LI", "TD"))
makeGroupNameVal <- function(l) {
  groups <- l[['groups']]
  vals <- l[['vals']]
  vals <- strsplit(vals, "\\.")
  nn <- lapply(vals, function(x) {
    paste0(groups, c(": "), x)
  })
  nn
}

## From hadds
## Names of dots
dots <- function(...) {
  eval(substitute(alist(...)))
}

compact <- function(x) Filter(Negate(is.null), x)


## Capture function arg, build function inside
curve2Fun <- function(curve) {
  arggs <- as.list(curve)[-1]
  cname <- deparse1(as.list(curve)[[1]])
  cbody <- body(get(cname))
  cformal <- formals(get(cname))
  for(nn in names(arggs)) {
    cformal[[nn]] <- arggs[[nn]]
  }
  idx <- which(names(cformal) == "...")
  if  (length(idx)) {
    ss <- seq_along(cformal)[-idx]
    cformal <- cformal[c(ss, idx)]
  }
  w <- function() {}
  body(w) <- cbody
  formals(w) <- cformal
  environment(w) <- new.env()
  w
}

## Verify this works like above, need to add part to adjust for `...`
# call2Fun <- function(x) {
#   arggs <- as.list(x)[-1] # calls can be treated as lists
#   cname <- deparse1(as.list(x)[[1]])
#   w <- get(cname)
#   for (nn in names(arggs)) {
#     formals(w)[[nn]] <- arggs[[nn]]
#   }
#   w
# }


## Pull from attributes names of vars to split
# data by observation (subject, group)
getSplitVars <- function(bdObj) {
  bdCall <- attr(bdObj, "call")
  nn <- c(eval(bdCall[['subject']]), eval(bdCall[['group']]))
  nn
}

## Correctly subsets dataset for observation
getSubX <- function(bdo) {
  X <- setDT(attr(bdo, "X")$X)
  nn <- getSplitVars(bdo)
  Xnames <- do.call(paste, X[, nn, with = FALSE])
  bdNames <- do.call(paste, bdo[, nn, with = FALSE])
  x_idx <- Xnames %in% bdNames
  X[x_idx, ]
}
















