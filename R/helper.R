
## Make split retain bdotsObj class
# Need to also split data attribute
split.bdotsObj <- function(x, by, ...) {
    oldClass <- class(x)
    #grps <- eval(attr(bdObj, "call")[['group']])
    class(x) <- c("data.table", "data.frame")
    res <- lapply(split(x, by = by, ...), function(y) {
        class(y) <- oldClass
        y
    })
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

## Get standard deviation for t statistic
# in unpaired test

nopairSD <- function(l) {
  if(length(l) != 2) stop("contact author with 123")
  s <- lapply(l, function(x) {
    vv <- x[['sd']]^2
    n <- x[['n']]
    (n-1) * vv
  })
  s <- Reduce(`+`, s) * sqrt(1/l[[1]]$n + 1/l[[2]]$n)
}

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











