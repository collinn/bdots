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
#' an equivalent? I know BH does, but the rest are not linear adjustments.
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
      message('rho not assigned with method "olseon". Using rho = 0.9')
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
  } else {
    if (missing(method)) {
      method <- stats::p.adjust.methods
    }
    adjpval <- stats::p.adjust(p, method, n)
    attr(adjpval, "alphastar") <- NA
    attr(adjpval, "rho") <- NA
  }
  return(adjpval)
}
