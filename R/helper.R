


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

## Make split retain bdotsObj class
# Need to also split data attribute
split.bdotsObj <- function(bdo, by, ...) {
  oldAttr <- attributes(bdo)
  res <- lapply(data.table:::split.data.table(bdo, by = by, ...), function(x) {
    attributes(x) <- oldAttr
    x
  })
  structure(.Data = res, class = c("bdObjList"))
}


## Otherwise, rbindlist is not a generic
rbindlist <- function(x, ...) {
  UseMethod("rbindlist")
}

rbindlist.list <- function(x, ...) {
  data.table::rbindlist(x, ...)
}

rbindlist.bdObjList <- function(bdo, ...) {
  oldAttr <- attributes(bdo[[1]])
  class(bdo) <- "list"
  bdo <- rbindlist(bdo)
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
  ss <- seq_along(cformal)[-idx]
  cformal <- cformal[c(ss, idx)]
  w <- function() {}
  body(w) <- cbody
  formals(w) <- cformal
  w
}

## Pull from attributes names of vars to split
# data by observation (subject, group)
getSplitVars <- function(bdo) {
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

## Takes value from bdotsFitter
# returns an appropriate DT
## WE SHOULDN'T NEED THIS
# fitObj <- res[[1]]
# bdFit2DT <- function(fitObjs) {
#   fitList <- lapply(fitObjs, function(fitObj) {
#     fn <- fitObj$fitName
#     dat <- as.data.table(matrix(c(fn, c("fit", "R2", "AR1", "fitCode")),
#                                 ncol = length(fn) + 4))
#     names(dat) <- c(names(fn), c("fit", "R2", "AR1", "fitCode"))
#     dat$fit <- I(list(fitObj['fit']))
#     dat$R2 <- fitObj[['R2']]
#     dat$AR1 <- (fitObj[['fitCode']] < 3)
#     dat$fitCode <- fitObj[['fitCode']]
#     dat
#   })
#   fitList <- rbindlist(fitList)
#   fitList[, fitCode := factor(fitCode, levels = 0:6)]
#   fitList
# }






