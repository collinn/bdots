## Methods for bdots objects (not including plots or summary)
## also contains p_adjust

## ----------

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

#' Extract bdotsFit Moedel Coefficients
#'
#' Returns coefficient matrix for bdotsFit object
#'
#' @param object A bdotsObj
#' @param ... not used
#'
#' @return Returns matrix of model coefficients for observations in \code{object}
#'
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

## I need to make sure these make sense
# specifically, I need the documentation to render for either when ?rbindlist called
#' @export
rbindlist <- function(x) UseMethod("rbindlist")

#' ## Don't export for now because fuck S3 generic matching
#' rbindlist <- function(x, ...) {
#'   UseMethod("rbindlist")
#' }
#'
#' @importFrom data.table rbindlist
#' @export
rbindlist.default <- function(x, ...) {
  data.table::rbindlist(x, ...)
}
#'
#' @importFrom data.table rbindlist
#' @export
rbindlist.list <- function(x, ...) {
  data.table::rbindlist(x, ...)
}

## Not 100% sure I should include this
#' rbindlist for bdotsObjects
#'
#' Similar to data.table::rbindlist, but preserves botsObjects attributes
#'
#' @param x bdotsObject
#' @param ... for compatability with data.table
#'
#' @export
rbindlist.bdObjList <- function(x, ...) {
  oldAttr <- attributes(x[[1]])
  #class(x) <- "list"
  x <- data.table::rbindlist(x, ...)
  attributes(x) <- oldAttr
  x
}


#' Adjust P-values for Multiple Comparisons
#'
#' Identical to \code{stats::p.adjust}, but includes \code{method = "oleson"}
#'
#' @param p numeric vector of p-values (possibly with NAs).
#' @param method correction method, a character string. Can be any of the methods in
#' p.adjust.methods, with the additional value \code{method = "oleson"}
#' @param n number of comparisons, must be at least \code{length(p)}; only set this
#' (to non-default) when you know what you are doing!
#' @param alpha adjustment to be made with method oleson
#' @param df degrees of freedom, if using \code{method = "oleson"}
#' @param rho AR1 correlation coefficient, if using \code{method = "oleson"}
#' @param cores number of cores for use in parallel, only valid for
#' \code{method = "oleson"}. Default is zero, using half of the available cores
#'
#' @details This function works identically to the function \code{p.adjust}, with
#' the additional option to use \code{method = "oleson"}. For this option, user
#' must include a value for \code{df}, \code{alpha}. If \code{method = "oleson"} and
#' no value is given for \code{rho}, 0.9 will be used. To compute a value for \code{rho}
#' from t-statistics, use \code{ar1Solver}.
#'
#'
#' @return Returns a vector of adjusted p-values just as in \code{p.adjust}, but
#' with additional attributes for alphastar and rho.
#'
#' @seealso \code{\link[bdots]{ar1Solver}}
#' @import stats
#'
#' @export
p_adjust <- function(p, method = "oleson", n = length(p),
                     alpha = 0.05, df, rho, cores = 0) {

  method <- match.arg(method, c("oleson", stats::p.adjust.methods))

  if (method == "oleson") {
    if (missing(df)) stop('Require value for df when using method "oleson"')
    if (cores < 1) cores <- detectCores()/2
    if (missing(rho)) {
      message('rho not assigned with method "olseon". Using rho = 0.9')
      rho <- 0.9
    }

    alphastar <- findModifiedAlpha(rho, n, df, alpha, cores = cores)
    k <- alphastar/alpha
    adjpval <- pmin(p/k, 1)
    attr(adjpval, "alphastar") <- alphastar
    attr(adjpval, "rho") <- rho
  } else {
    adjpval <- stats::p.adjust(p, method, n)
    if (method == "bonferroni") {
      attr(adjpval, "alphastar") <- alpha/n
    } else {
      attr(adjpval, "alphastar") <- NA
    }
    attr(adjpval, "rho") <- NA
  }
  return(adjpval)
}

#' Print `bdotsBootObj`
#'
#' Generic for printing `bdotsBootObj`
#'
#' @param x An object of class `bdotsBootObj`
#' @param ... Top secret alpha one code red
#'
#' @description Prints argument. Really, just the summary function
#'
#' @export
print.bdotsBootObj <- function(x, ...) {
  y <- summary(x)
  print(y)
}
