#' Logistic curve function for nlme
#'
#' Logistic function used in fitting nlme curve for observations
#'
#' @param dat subject data to be used
#' @param y outcome variable
#' @param time time variable
#' @param params \code{NULL} unless user wants to specify starting parameters for gnls
#' @param ... just in case
#'
#' @details \code{y ~ mini + (peak - mini) / (1 + exp(4 * slope * (cross - (time)) / (peak - mini)))}
#' @export
logistic <- function(dat, y, time, params = NULL, ...) {

  # if (getFormulaOnly) {
  #   y <- str2lang(y)
  #   time <- str2lang(time)
  #   ff <- bquote(.(y) ~ mini + (peak - mini) / (1 + exp(4 * slope * (cross - (.(time))) / (peak - mini))))
  #   return(ff)
  # }

  ## Make functions subroutines so that
  # 1) pass less shit into parallel
  # 2) users will need to make inside anyways
  logisticPars <- function(dat, y, time, ...) {
    time <- dat[[time]]
    y <- dat[[y]]

    mini <- min(y)
    peak <- max(y)
    r <- (peak - mini)
    cross <- time[which.min(abs(0.5*r - y))]

    # slope
    q75 <- .75 * r + mini
    q25 <- .25 * r + mini
    time75 <- time[which.min(abs(q75 - y))]
    time25 <- time[which.min(abs(q25 - y))]
    slope <- (q75 - q25) / (time75 - time25)

    return(c(mini = mini, peak = peak, slope = slope, cross = cross))
  }

  if (is.null(params)) {
    params <- logisticPars(dat, y, time)
  } else {
    if (length(params) != 4) stop("logistic requires 4 parameters be specified for refitting")
    if (!all(names(params) %in% c("mini", "peak", "slope", "cross"))) {
      stop("logistic parameters for refitting must be correctly labeled")
    }
  }
  y <- str2lang(y)
  time <- str2lang(time)
  ff <- bquote(.(y) ~ mini + (peak - mini) / (1 + exp(4 * slope * (cross - (.(time))) / (peak - mini))))
  attr(ff, "parnames") <- names(params)
  return(list(formula = ff, params = params))
}
