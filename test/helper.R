
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
isPaired <- function(l) { # this only works for lists
    reduceEq <- function(x, y) if(identical(x[["Subject"]], y[["Subject"]])) x else FALSE
    tmp <- Reduce(reduceEq, l)
    if(identical(tmp, FALSE)) FALSE else TRUE
}
