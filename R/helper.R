


## getVarMat
# takes subset data with single observation
# returns covariance matrix of fit parameters
getVarMat <- function(dat) {
  if(nrow(dat) != 1) stop("only for single row of bdotsObj")
  dat$fit[[1]]$varBeta
}

# coef <- function(x, ...) {
#   UseMethod("coef")
# }


## Extract coef from  bdotsObj
## uh, this doesn't adress null
# Ah, mother fucker, that's fitcode 6!
# this can potentially be made cleaner
#### Can't replace fit[[1]] since it's unnamed length 1 list. Could name it, I guess
## Needs to be made generic, address fitcode mentioned above.
# this needs to be renamed coef.bdObj
#' Extract bdotsFit Moedel Coefficients
#'
#' Returns coefficient matrix for bdotsFit object
#'
#' @param object A bdotsObj
#' @param ... not used
#' @import stats
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
#' @import data.table
split.bdotsObj <- function(bdo, by, ...) {
  oldAttr <- attributes(bdo)
  class(bdo) <- c("data.table", "data.frame")
  res <- lapply(split(bdo, by = by, ...), function(x) {
    attributes(x) <- oldAttr
    x
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


rbindlist <- function(x, ...) {
  UseMethod("rbindlist")
}

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
rbindlist.bdObjList <- function(bdo, ...) {
  oldAttr <- attributes(bdo[[1]])
  class(bdo) <- "list"
  bdo <- rbindlist(bdo, ...)
  attributes(bdo) <- oldAttr
  bdo
}


## Used in bdotsBoot
## Probably having subject in string is wrong, but ok for now
isPaired <- function(l) { # this only works for lists of bdObj
    reduceEq <- function(x, y) if(identical(x[["Subject"]], y[["Subject"]])) x else FALSE
    tmp <- Reduce(reduceEq, l)
    # !identical(tmp, FALSE) <- see if this works
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
  if(length(l) != 2) stop("contact author with 123")
  s <- lapply(l, function(x) {
    vv <- x[['sd']]^2
    n <- x[['n']]
    (n-1) * n * vv
  })
  s <- Reduce(`+`, s) / (l[[1]]$n + l[[2]]$n - 2) * (1/l[[1]]$n + 1/l[[2]]$n)
  s <- sqrt(s)
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

#' Adjust P-values for Multiple Comparisons
#'
#' Identical to stats::p.adjust, but includes \code{method = "oleson"}
#'
#' @param p numeric vector of p-values (possibly with NAs).
#' @param method correction method, a character string. Can be any of the methods in
#' p.adjust.methods, with the additional value \code{method = "oleson"}
#' @param n number of comparisons, must be at least \code{length(p)}; only set this
#' (to non-default) when you know what you are doing!
#' @param alpha adjustment to be made with method oleson
#' @param df degrees of freedom, if using \code{method = "oleson"}
#' @param rho correlation coefficient, if using \code{method = "oleson"}
#' @param cores number of cores for use in parallel, only valid for
#' \code{method = "oleson"}. Default is zero, using half of the available cores
#'
#' @details I will revist this function later to make sure the arguments are consistent
#' with how they may be used in the wild. Also need to include details from p.adjust in stats
#' I can add oleson there with there references. Need example too, I think. ALSO - it's sometimes
#' helpful to know alphastar. Should I return a similar value for the other methods? Do they have
#' an equivalent? I know BH does, but I'm not positive the rest are linear adjustments.
#' ALSO::: can ar1Solver determine rho from pvalue, or is tvalue needed?
#'
#' @import stats
#'
#' @export
p.adjust <- function(p, method = "oleson", n = length(p), alpha = 0.05,
                     df, rho, cores = 0) {
  if (method == "oleson") {
    if (cores < 1) cores <- detectCores()/2
    if (missing(rho)) {
      message('rho not assigned with method "olseon". Using rho <- 0.9')
      rho <- 0.9
    }
    if (missing(df)) {
      stop('Require value for df when using method "oleson"')
    }

    alphastar <- findModifiedAlpha(rho, n, df, cores)
    k <- alphastar/alpha
    adjpval <- p/k
    attr(adjpval, "alphastar") <- alphastar
    attr(adjpval, "rho") <- rho
    return(adjpval)
  } else {
    if (missing(method)) {
      method <- stats::p.adjust.methods
    }
    stats::p.adjust(p, method, n)
  }
}















