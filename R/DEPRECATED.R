
#' Fit nlme curves to grouped observations
#'
#' DEPRECATED. Use \code{bfit} instead
#'
#' @param data Dataset used
#' @param subject Column name of dataset containing subject identifiers
#' @param time Column name containing time variable
#' @param y Column name containing outcome of interest
#' @param group Character vector containing column names of groups. Can be
#' greater than one
#' @param curveType See details/vignette
#' @param cores number of cores. Default is \code{0}, indicating half cores available
#' @param cor autocorrelation TRUE/FALSE
#' @param ar Value indicates estimate for autocorrelation. A value of zero indicates to fit without AR(1) assumption
#' @param numRefits Integer indicating number of attempts to fit an observation
#' if the first attempt fails
#' @param verbose currently not used
#' @param ... Secret
#'
#' @return Object of class 'bdotsObj', inherits from data.table
#'
#' @details This function is being deprecated. Use \code{bfit} instead
#'
#' @examples
#' \dontrun{
#' res <- bdotsFit(data = cohort_unrelated,
#'                 subject = "Subject",
#'                 time = "Time",
#'                 y = "Fixations",
#'                 group = c("Group", "LookType"),
#'                 curveType = doubleGauss(concave = TRUE),
#'                 numRefits = 2,
#'                 cores = 0,
#'                 verbose = FALSE)
#' }
#'
#' @import data.table
#' @import parallel
#' @importFrom utils object.size
#' @export
bdotsFit <- function(data, # dataset
                     subject, # subjects
                     time, # column for time
                     y, # response vector
                     group, # groups for subjects
                     curveType = doubleGauss(concave = TRUE),
                     cor = TRUE,
                     ar = 0, # autocorrelation?
                     numRefits = 0,
                     cores = 0, # cores to use, 0 == 50% of available
                     verbose = FALSE,
                     ...) {
  tt <- match.call(expand.dots = TRUE)
  .Deprecated("bfit")
  do.call("bfit", as.list(tt)[-1])
}


#' Create bootstrapped curves from bdotsObj
#'
#' DEPRECATED. Use \code{bboot} instead
#'
#' @param formula See details.
#' @param bdObj An object of class 'bdotsObj'
#' @param Niter Number of iterations of bootstrap to draw
#' @param alpha Significance level
#' @param padj Adjustment to make to pvalues for significance. Will be able to
#' use anything from \code{p.adjust} function, but for now, just "oleson"
#' @param permutation Boolean indicating whether to use permutation testing rather
#' thank adjusting alpha to control FWER. WARNING: This option is very much in beta testing
#' and not recommended for general use. Also not available for paired tests or difference of difference
#' @param cores Number of cores to use in parallel. Default is zero, which
#' uses half of what is available.
#' @param skipDist do not use
#' @param singleMeans definitely do not use
#' @param permAddVar Boolean. Add observed variability for perms?
#' @param ... not used
#'
#' @details Deprecated. Use \code{bboot} instead
#' @import data.table
#' @export
bdotsBoot <- function(formula,
                      bdObj,
                      Niter = 1000,
                      alpha = 0.05,
                      padj = "oleson",
                      permutation = FALSE, skipDist = FALSE,
                      singleMeans = FALSE, permAddVar = TRUE,
                      cores = 0, ...) {
  tt <- match.call(expand.dots = TRUE)
  .Deprecated("bboot")
  do.call("bboot", as.list(tt)[-1])
}


#' Refit Observations Returned from bdotsFit
#'
#' DEPRECATED. Use \code{brefit} instead
#'
#' @param bdObj An object of class 'bdotsObj' returned from \code{bdotsFit}
#' @param fitCode A length one integer indicating observations to refit. See Details
#' @param subset Either an expression that evaluates to a logical used to subset the \code{bdObj},
#' (using \code{data.table} syntax) or a numeric vector of indices to subset. Default is \code{NULL}.
#' When not \code{NULL}, any arguments to \code{fitCode} are ignored.
#' @param quickRefit Boolean indicating if a quick refit should be used. If TRUE,
#' rather than prompting the user for adjustments for each observation, \code{bdotsReft}
#' will jitter the parameters of all observations indicated by \code{fitCode} and attempt
#' to refit. Between the original and the refitted curve, this will place priority on
#' the higher \code{fitCode}. If these are equal, R2 will take precedence. Otherwise,
#' the original fit will be kept.
#' @param numRefits Integer indicating the number of refit attempts after jittering
#' parameters, either with quickRefit or when done individually
#' @param paramDT A \code{data.table} or \code{data.frame} that matches the what is
#' returned by \code{coefWriteout(bdObj)}. That is, it should have columns
#' uniquely identifying observations with subjects and groups, as well as named
#' columns for the paramters. NA parameters are OK. Can also be a subset of the original rows.
#' Note, if this argument is not \code{NULL}, the remaining arguments will be ignored.
#' @param ... not used
#'
#' @return Returns bdObj with updated fits
#'
#' @details DEPRECATED. Use \code{brefit} instead
#'
#' @import data.table
#' @export
bdotsRefit <- function(bdObj, fitCode = 1L, subset = NULL, quickRefit = FALSE,
                       numRefits = 2L, paramDT = NULL, ...) {
  tt <- match.call(expand.dots = TRUE)
  .Deprecated("brefit")
  do.call("brefit", as.list(tt)[-1])
}


#' Correlation with fixed value in bdots
#'
#' DEPRECATED. Use \code{bcorr} instead.
#'
#' @param bdObj Object of class `bdotsObj`
#' @param val Character string of fixed value for correlation in dataset from `bdotsFit`
#' @param ciBands Boolean for including confidence intervals
#' @param method Arguments for `cor` or `cor.test`. The default option us `method = "pearson"`
#'
#' @export
bdotsCorr <- function(bdObj, val, ciBands = FALSE, method = "pearson") {
  tt <- match.call(expand.dots = TRUE)
  .Deprecated("bcorr")
  do.call("bcorr", as.list(tt)[-1])
}
