
## Make split retain bdotsObj class

split.bdotsObj <- function(x, ...) {
    oldClass <- class(x)
    class(x) <- c("data.table", "data.frame")
    res <- lapply(split(x, ...), function(y) {
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
















