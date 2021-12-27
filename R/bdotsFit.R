#' Fit nlme curves to grouped observations
#'
#' Creates observation level curves to use in bdotsBoot
#'
#' @param data Dataset used
#' @param subject Column name of dataset containing subject identifiers
#' @param time Column name containing time variable
#' @param y Column name containing outcome of interest
#' @param group Character vector containing column names of groups. Can be
#' greater than one
#' @param curveType See details/vignette
#' @param cores number of cores. Default is \code{0}, indicating half cores available
#' @param cor Boolean. Autocorrelation?
#' @param numRefits Integer indicating number of attempts to fit an observation
#' if the first attempt fails
#' @param verbose currently not used
#' @param ... Secret
#'
#' @return Object of class 'bdotsObj', inherits from data.table
#'
#' @details This is step one of the three step bdots process. Things should be
#' more or less straight forward. The only tricky part involves curveType. For now
#' know that one can use doubleGauss(concave = TRUE/FALSE) or logistic(). Should
#' be passed in as a call. See the vignette on customizing this
#'
#' @examples
#' \dontrun{
#' res <- bdotsFit(data = cohort_unrelated,
#'                 subject = "Subject",
#'                 time = "Time",
#'                 y = "Fixations",
#'                 group = c("Group", "LookType"),
#'                 curveType = doubleGauss(concave = TRUE),
#'                 cor = TRUE,
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
                     cor = TRUE, # autocorrelation?
                     numRefits = 0,
                     cores = 0, # cores to use, 0 == 50% of available
                     verbose = FALSE,
                     ...) {

  if (cores < 1) cores <- detectCores()/2
  curveType <- substitute(curveType)
  curveName <- curveType[[1]]
  curveType <- curve2Fun(curveType)

  ## Variable names on the dataset
  datVarNames <- c(y = y, subject = subject, time = time, group = group)
  haveVars <- datVarNames %in% names(data)
  if (!all(datVarNames %in% names(data))) {
    stopMsg <- paste0('"', datVarNames[!haveVars],'" is not a column of the dataset')
    stop(stopMsg)
  }

  # for removing rho (need to adjust in case it's in (...))
  # Should verify that we don't have rho set and cor = FALSE
  if (!exists("rho")) {
    rho <- ifelse(cor, 0.9, 0)
  } else {
    if (cor & (rho >= 1 | rho < 0)) {
      warning("cor set to TRUE with invalid rho. Setting rho to 0.9")
      rho <- 0.9
    }
  }

  ## Factors are bad, m'kay?
  dat <- setDT(data)
  dat[, (group) := lapply(.SD, as.character), .SDcols = group]
  dat[, (subject) := lapply(.SD, as.character), .SDcols = subject]

  ## Let's only keep the columns we need (have not tested this yet)
  dat <- dat[, c(y, time, subject, group), with = FALSE]


  timetest <- split(dat, by = group, drop = TRUE)
  timetest <- lapply(timetest, function(x) unique(x[[time]]))
  timeSame <- identical(Reduce(intersect, timetest, init = timetest[[1]]),
                        Reduce(union, timetest, init = timetest[[1]]))
  if (!timeSame) stop("Observed times are different between groups")

  ## This should work inside function. Let's check
  # if this happens, need to modify X for plots to work
  if (any(dat[, .N, by = c(subject, time, group)]$N > 1)) {
    warning("Some subjects have multiple observations for unique time. These will be averaged")
    yval <- deparse(substitute(y))
    dat[, substitute(y) := mean(get(y)), by = c(subject, time, group)]
    dat <- unique(dat, by = c(subject, time, group))
  }

  splitVars <- c(subject, group)
  newdat <- split(dat, by = splitVars, drop = TRUE)

  if (Sys.info()['sysname'] == "Darwin") {
    cl <- makePSOCKcluster(cores, setup_strategy = "sequential")
  } else {
    cl <- makePSOCKcluster(cores)
  }

  invisible(clusterEvalQ(cl, {library(bdots)}))

  ## Only for error checking, as devtools::load_all doesn't help
  ## with parallel
  # res <- vector("list", length(newdat))
  # for (i in seq_along(newdat)) {
  #   res[[i]] <- bdotsFitter(newdat[[i]], curveType, rho,
  #                           numRefits, verbose = FALSE,
  #                           splitVars = splitVars, datVarNames = datVarNames)
  # }

  res <- parLapply(cl, newdat, bdotsFitter,
                   curveType = curveType,
                   rho = rho, numRefits = numRefits,
                   verbose = FALSE,
                   splitVars = splitVars,
                   datVarNames = datVarNames)
  stopCluster(cl)

  ## Remove entries that had zero outcome variance
  ## NOTE:: need to make option to remove all entries from subject
  # to retain paired t-test
  nullEntries <- vapply(res, is.null, logical(1))
  if (sum(nullEntries) != 0) {
    nn <- names(newdat)[which(nullEntries)]
    message(paste0("Observations with zero outcome variance were removed (n = ", length(nn), "):"))
    for (i in seq_along(nn)) {
      message(nn[i])
    }
    res <- compact(res)
  }

  ff <- attr(res[[1]], "formula")
  fitList <- rbindlist(res, fill = TRUE)

  ## If too large, should get "name" of data
  ## and option to call it from global env
  # if (is.null(returnX)) {
  #   sz <- object.size(dat)
  #   if (sz < 1e8L) X <- dat
  # } else if (returnX) {
  #   X <- dat
  # } else {
  #   X <- dat # for now
  # }

  X_env <- new.env(parent = emptyenv())
  X_env$X <- dat

  ## Janky for now, but we want a groupName List
  vals <- do.call(function(...) paste(..., sep = "."),
                  unique(fitList[, group, with = FALSE]))
  groups <- list(groups = group,
                 vals = vals)

  res <- structure(class = c("bdotsObj", "data.table", "data.frame"),
                   .Data = fitList,
                   formula = ff,
                   curveType = curveName,
                   call = match.call(),
                   time = timetest[[1]],
                   rho = rho,
                   groups = groups,
                   X = X_env)
}






























