## This document should only contain interally used helper functions
# not anything that is eventually exported

## ----------

## Two functions to bind parameters to unique identifiers
#' Create \code{data.table} with \code{bdotsObj} parameters
#'
#' Creates an object of class \code{data.table} that matches
#' parameter values for each observation. This can then be
#' passed to the \code{bdotsRefit} function
#'
#' @param bdObj An object returned from \code{bdotsFit} or \code{bdotsRefit}
#'
#' @return A \code{data.table} matching parameter values to observations
#'
#' @examples
#' \dontrun{
#' fit <- bdotsFit(data = cohort_unrelated,
#'                 subject = "Subject",
#'                 time = "Time",
#'                 y = "Fixations",
#'                 group = c("Group", "LookType"),
#'                 curveType = doubleGauss(concave = TRUE),
#'                 cor = TRUE,
#'                 numRefits = 2,
#'                 cores = 0,
#'                 verbose = FALSE)
#' parDT <- coefWriteout(fit)
#' }
#' @export
coefWriteout <- function(bdObj) {
  cmat <- coef(bdObj)
  idcols <- getIdentifierCols(bdObj)
  idcols <- bdObj[, idcols, with = FALSE]
  res <- cbind(idcols, cmat) # attributes not preserved when writing out so don't add
}

getIdentifierCols <- function(bdo) {
  sub <- attr(bdo, "call")$subject
  grps <- eval(attr(bdo, "call")$group)
  c(sub, grps)
}



## getVarMat
# takes subset data with single observation
# returns covariance matrix of fit parameters
getVarMat <- function(dat) {
  if(nrow(dat) != 1) stop("only for single row of bdotsObj")
  dat$fit[[1]]$varBeta
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

## Taken from purrrrrrrr
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

## A trillion times faster with rlang::squash. Sad
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


## Pull from attributes names of vars to split
# data by observation (subject, group)
# this is exactly getIdentifierCols. Get rid of one of these. Probably this one
# since the other is more generally named
getSplitVars <- function(bdObj) {
  bdCall <- attr(bdObj, "call")
  nn <- c(eval(bdCall[['subject']]), eval(bdCall[['group']]))
  nn
}

## Correctly subsets dataset for observation
getSubX <- function(bdo) {
  X <- setDT(attr(bdo, "X")$X)
  nn <- getIdentifierCols(bdo)
  Xnames <- do.call(paste, X[, nn, with = FALSE])
  bdNames <- do.call(paste, bdo[, nn, with = FALSE])
  x_idx <- Xnames %in% bdNames
  X[x_idx, ]
}

## Create curve function from formula
# used in bdotsBoot and bdUpdate_NULL
makeCurveFun <- function(bdObj) {
  time <- attr(bdObj, "call")[['time']]
  f_bod <- attr(bdObj, "formula")[[3]]
  parnames <- attributes(attr(bdObj, "formula"))[['parnames']]
  f_args <- c(parnames, time)
  f_args <- setNames(as.pairlist(rep("", length(f_args))), f_args)
  eval(call("function", f_args, f_bod), parent.frame())
}


## Capture function arg, build function inside, bdotsFit
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


## Standard deviation stuff. Not 1000% sure, but I think ok for now
###########################################
#####
## https://en.wikipedia.org/wiki/Student%27s_t-test
#####
### There are a bunch that are possible
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
# nopairSD2 <- function(l) {
#   sd_ratio <- Reduce(function(x, y) {x$sd/y$sd}, l)
#   ## Proportion of sd ratio within bounds should be some val, lets say 0.5
#   var_sim <- mean(sd_ratio > 1/2 & sd_ratio < 2) > 0.5
#   if (var_sim) {
#     s <- lapply(l, function(x) {
#       vv <- x[['sd']]^2
#       n <- x[['n']]
#       (n - 1) * vv
#     })
#     s <- Reduce(`+`, s) / (l[[1]]$n + l[[2]]$n - 2) * (l[[1]]$n + l[[2]]$n)
#     s <- sqrt(s)
#     dof <- l[[1]]$n + l[[2]]$n - 2
#   } else {
#     s <- lapply(l, function(x) {
#       vv <- x[['sd']]^2
#       n <- x[['n']]
#       vv/n
#     })
#     dof_denom <- Reduce(`+`, Map(function(x, y) {x^2 / (y[['n']] - 1)}, x = s, y = l))
#     s <- Reduce(`+`, s)
#     dof <- s^2 / dof_denom
#     s <- sqrt(s)
#   }
#   return(list(sd = s, dof = dof))
# }
