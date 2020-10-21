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
  curveName <- gsub("\\(|\\)", "", deparse1(curveType))
  curveType <- curve2Fun(curveType) # this function has an environment with a lot of stuff in it

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


  ## We can quickly check that time is kosher, at least within group divides
  ## Would need to verify this is we planned on doing paired stuff?
  ## Here's the thing, with some maneuvering, I could verify how many
  ## unique sets of time we have, then attach only that as list in attributes?
  ## For now, let's just PRAY that they are all equal. We just need it to drop
  # into bdotsBoot.
  timetest <- split(dat, by = group, drop = TRUE)
  timetest <- lapply(timetest, function(x) unique(x[[time]]))
  timeSame <- identical(Reduce(intersect, timetest, init = timetest[[1]]),
                        Reduce(union, timetest, init = timetest[[1]]))
  if (!timeSame) warning("Yo, these times are different between groups, and until collin fixes it, that's going to make bdotsBoot wrong-ish")

  ## This shuld work inside function. Let's check
  if (any(dat[, .N, by = c(subject, time, group)]$N > 1)) {
    warning("Some subjects have multiple observations for unique time. These will be averaged")
    yval <- deparse(substitute(y))
    dat[, substitute(y) := mean(get(y)), by = c(subject, time, group)]
    dat <- unique(dat, by = c(subject, time, group))
  }



 # #  ## if(.platform$OStype == windows)
 # #  ## This allows output to be print to console, possibly not possible in windows
 cl <- makePSOCKcluster(cores)
 ## Ideally, I could clusterEvalQ bdots
 #clusterExport(cl, c("curveFitter", "makeCurveEnv", "dots", "compact"), envir = parent.frame())
 invisible(clusterEvalQ(cl, {library(nlme); library(bdots)}))
 #invisible(clusterEvalQ(cl, library(bdots))) # someday

 splitVars <- c(subject, group)
 newdat <- split(dat, by = splitVars, drop = TRUE) # these needs to not be in string
 # res <- parLapply(cl, newdat, bdotsFitter,
 #                  curveList = curveList,
 #                  rho = rho, numRefits = numRefits,
 #                  verbose = FALSE)
 res <- parLapply(cl, newdat, bdotsFitter,
                  curveType = curveType,
                  rho = rho, numRefits = numRefits,
                  verbose = FALSE,
                  splitVars = splitVars,
                  datVarNames = datVarNames)
 ## This allows us to pass names for verbose
 #system.time(res <- clusterMap(cl, bdotsFitter, newdat, thenames = names(newdat)))
 stopCluster(cl)


 ## See, I don't like this bc hypothetically could be length 0
 # alternative is making logistic() more complicated
 ff <- attr(res[[1]], "formula")
 fitList <- rbindlist(res, fill = TRUE)

 ## I maybe don't like this here. But we can add it on
 # in the summary if necessary
 #fitList[, fitCode := factor(fitCode, levels = 0:6)]

 ## Dude, just store that covariate table they want in a list as well
 # and that can be it's own class if need be
 ## Neat idea - environment to data.table to work with objects like this
 # that are actually environments, but behave like data.table
 # useful!
 ## At some point, need to change the way this returns fit entry as "AsIs" object, nested list
 ## YO, mother-effer. We are doing this to each element of a list with
 # information from that list. In other words, this whole function
 # right here, everything that happens, we can throw that shit
 # into bdotsFitter. It has no place here. Goodness gracious. That's
 # why it'll be important to have the curvefunction pass something that
 # gives us names for the parameters!
 ## Ok, sure, but that messes with renaming c(subject, group)
 # which I don't want to  have to pass to bdotsFitter. I'll leave as is for now.
 ## So, to-do here:
 # remove matrix as attachment to it
 # consider putting it in attributes? This would
 # assist with subsetting and letting them adjust it
 # manually. We will put more thought in this later

 ## This could be a trillion times faster (maybe) with loop and data.table::set
 # fitList <- lapply(names(newdat), function(x) {
 #   result <- res[[x]] # list of length 3
 #   x <- strsplit(x, "\\.") # list of by variables for newdat
 #
 #   dat1 <- as.data.table(matrix(x[[1]], ncol = length(x[[1]])))
 #   names(dat1) <-  c(subject, group)
 #
 #
 #   #set(dat1, j = "fit", value = result[['fit']])
 #   dat1$fit <- I(list(result['fit']))
 #   dat1$R2 <- result[['R2']]
 #   dat1$AR1 <- (result[['fitCode']] < 3)
 #   dat1$fitCode <- result[['fitCode']]
 #   dat1
 # })
 # fitList <- rbindlist(fitList)
 #



 ## Return data matrix as well?
 ## Two things here
 # don't make X AND data, even if at the end, unneccesary size
 # if size too large, "get" name of data, and store that
 # then `getSubX` can either give the actual data set
 # or it can find it in the globalenv and use it there.
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

 ## Janky for now, but we want a gropuName List
 ## TOO janky. Fix this ASAP
 tmp <- unique(fitList[, group, with = FALSE])
 tmp <- names(split(tmp, by = group))
 groups <- list(groups = group,
                vals = tmp)

 ## Class bdots object

 ## 1) curveType
 ## 2) formula
 ## 3) fitList
   # i) columns for identifiers (subject, group)
   # ii) gnls model fit
   # iii) R2, AR1 for visual reference
   # iv) fitCode for refitting

 ## To refit
 # 1) Punch your own gnls object into fit (on your own for this one,
      # can use formula from bdots object)
 # 2) This process
  # i) fit <- bdotsFit(...)
  # ii) examine fits w/ plots or w/e
     # I) refits for all fitCode > n
     # II) specify specific entry to refit
  # iii) Perform refits
    # I) do it interactively (with a wizard!)
    # II) input your own starting parameters (as matrix)
       # a) must be same dimension as refits doing

 ## Returned object is just a data.table with attributes
 # this will make it significantly easier to work with/substitute
 # it will also keep formula/curveType with subsets
 ## Add number of parameters so that we know we are plotting last p columns

 ## We can extract this SUPER easily since it's named with
 # mod <- attr(res, "call")$curveType to get logistic(), poly(n), etc.,
 ## one thing to keep in mind with putting parameters in attributes
 # is retaining them upon subset
 ## Shit, yo, just make my own `[.bdotsObj` that inherits from data.table
 # preserves the matrix
 # same thing to be on lookout for when merging
 # rbind.bdotsObj == stop("don't do this, use merge")
 # shit, yo, could also use `makeActiveBinding()`
 # at some point, just consider R4 objects <- hard no on this
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






























