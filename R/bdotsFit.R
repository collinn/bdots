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
#' if the first attemp fails
#' @param verbose currently not used
#' @param returnX Boolean. Return data with bdObj? Ignore for now. There will be
#' a way around this later, but for now, if this is FALSE, things won't work
#' @param ... Secret
#'
#' @return Object of class 'bdotsObj', inherits from data.table
#'
#' @details This is step one of the three step bdots process. Things should be
#' more or less straight forward. The only tricky part involves curveType. For now
#' know that one can use doubleGauss(concave = TRUE/FALSE) or logistic(). Should
#' be passed in as a function. See the vignette on customizing this
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
                     returnX = NULL,
                     ...) {

  if (cores < 1) cores <- detectCores()/2
  curveType <- substitute(curveType)
  #curveName <- gsub("\\(|\\)", "", deparse1(curveType))
  curveName <- curveType[[1]]
  curveType <- curve2Fun(curveType)

  ## Variable names on the dataset
  datVarNames <- c(y = y, subject = subject, time = time, group = group)
  if (!all(datVarNames %in% names(data))) {
    stop("need more specific error. Either subject, time, y, or group is not in dataset")
  }

  # for removing rho (need to adjust in case it's in (...))
  # Should verify that we don't have rho set and cor = FALSE
  if (!exists("rho")) {
    rho <- ifelse(cor, 0.9, 0)
  } else {
    if(cor & (rho >= 1 | rho < 0)) {
      warning("cor set to TRUE with invalid rho. Setting rho to 0.9")
      rho <- 0.9
    }
  }

  ## Factors eff things up
  dat <- setDT(data)
  dat[, (group) := lapply(.SD, as.character), .SDcols = group]
  dat[, (subject) := lapply(.SD, as.character), .SDcols = subject]

  timetest <- split(dat, by = group, drop = TRUE)
  timetest <- lapply(timetest, function(x) unique(x[[time]]))
  timeSame <- identical(Reduce(intersect, timetest, init = timetest[[1]]),
                        Reduce(union, timetest, init = timetest[[1]]))
  if (!timeSame) stop("Yo, these times are different between groups, and until collin fixes it, that's going to make bdotsBoot wrong-ish")

  ## This shuld work inside function. Let's check
  if (any(dat[, .N, by = c(subject, time, group)]$N > 1)) {
    warning("Some subjects have multiple observations for unique time. These will be averaged")
    yval <- deparse(substitute(y))
    dat[, substitute(y) := mean(get(y)), by = c(subject, time, group)]
    dat <- unique(dat, by = c(subject, time, group))
  }

  splitVars <- c(subject, group)
  newdat <- split(dat, by = splitVars, drop = TRUE)

 # ## if(.platform$OStype == windows)
 cl <- makePSOCKcluster(cores)

 invisible(clusterEvalQ(cl, {library(nlme); library(bdots)}))

 res <- parLapply(cl, newdat, bdotsFitter,
                  curveType = curveType,
                  rho = rho, numRefits = numRefits,
                  verbose = FALSE,
                  splitVars = splitVars,
                  datVarNames = datVarNames)
 stopCluster(cl)

 ff <- attr(res[[1]], "formula")
 fitList <- rbindlist(res, fill = TRUE)

 ## If too large, should get "name" of data
 # and option to call it from global env
 if (is.null(returnX)) {
   sz <- object.size(data)
   if (sz < 1e8L) X <- data
 } else if (returnX) {
   X <- data
 } else {
   X <- NULL # for now
 }

 X_env <- new.env(parent = emptyenv())
 X_env$X <- X

 ## Janky for now, but we want a groupName List
 ## TOO janky. Fix this ASAP
 tmp <- unique(fitList[, group, with = FALSE])
 tmp <- names(split(tmp, by = group))
 groups <- list(groups = group,
                vals = tmp)

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






























