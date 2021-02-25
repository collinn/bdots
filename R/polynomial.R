#' Polynomial curve function for nlme
#'
#' Logistic function used in fitting nlme curve for observations
#'
#' @param dat subject data to be used
#' @param y outcome variable
#' @param time time variable
#' @param degree degree of polynomial
#' @param raw Boolean, use raw polynomials?
#' @param params \code{NULL} unless user wants to specify starting parameters for gnls
#' @param ... just in case
#'
#' @details It's recommended that one uses raw polynomials for this function for
#' numerical stability. As inference is not performed on the parameters themselves,
#' this should have minimial consequences
#'
#' @details \code{y ~ mini + (peak - mini) / (1 + exp(4 * slope * (cross - (time)) / (peak - mini)))}
#' @export
polynomial <- function(dat, y, time, degree, raw = TRUE, params = NULL, ...) {

  polyPars <- function(dat, y, time, degree, raw) {
    pp <- lm(dat[[y]] ~ poly(dat[[time]], degree = degree, raw = raw))
    setNames(coef(pp), c(paste0("beta", seq(degree + 1L))))
  }

  if (is.null(params)) {
    params <- polyPars(dat, y, time, degree, raw)
  } else {
    if (length(params) != degree + 1L) stop("poly requires that number of params matches degree + 1 (for intercept)")
    if (!all(names(params) %in% paste0("beta", seq(degree + 1L)))) {
      stop("polynomial parameters must be labeled beta1, ..., beta[degree + 1]")
    }
  }

  time_names <- paste0("I(Time^", seq(degree + 1L) - 1L , ")")

  ff <- paste(names(params), time_names, sep = "*", collapse = "+")
  ff <- str2lang(ff)
  y <- str2lang(y)
  ff <- bquote(.(y) ~ .(ff))
  attr(ff, "parnames") <- names(params)
  return(list(formula = ff, params = params))
}

