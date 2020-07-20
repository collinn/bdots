library(data.table)
library(magrittr)
library(parallel)

## currently, only actually using
# data.table
# parallel
# nlme




# pair verbose with message

## Need to not make my own names for this (i.e., y, subject, time)



## Yo, if we end up having to delete subjects, how do we do a paired t test?
## concave is doubleGauss onlly. Surely we can do better
bdotsFit <- function(data, # dataset
                     subject, # subjects
                     time, # column for time
                     y, # response vector
                     group, # groups for subjects
                     curveType = doubleGauss(concave = TRUE),
                     cor = TRUE, # autocorrelation?
                     refits = 0,
                     cores = 0, # cores to use, 0 == 50% of available
                     verbose = FALSE,
                     ...) {

  if (cores < 1) cores <- detectCores()/2


  # for removing rho (need to adjust in case it's in (...))
  rho <- ifelse(cor, 0.9, 0)


  ## Cheat around DT reference, conditionally set key for subset
  ## Can possibly avoid this set, as was done below
  # and get rid of magrittr
  # maybe get rid of data table :(
  ## For doubleGauss
  ## Names actually might not matter if I don't make bdotsFitter public
  dat <- data.table()
  dat$subject <- data[[subject]]
  dat$time <- data[[time]] %>% as.numeric()
  dat$y <- data[[y]] %>% as.numeric()
  #group <- c("Group", "LookType")

  ## Set group variables in data.tableb+
  for(gg in seq_along(group)) {
    set(dat, j = group[gg], value = data[[group[gg]]])
  }
#
# # for logistic
#     data(ci)
#     ci <- as.data.table(ci)
#     ci <- ci[LookType == "Target", ]
#     group <- "protocol"
#     dat <- data.table()
#     dat$subject <- ci[[subject]]
#     dat$time <- ci[[time]] %>% as.numeric()
#     dat$y <- ci[[y]] %>% as.numeric()
#     for(gg in seq_along(group)) {
#       set(dat, j = group[gg], value = ci[[group[gg]]])
#     }



  # Named list with args
  # but look in parser.R to see how to put this in environment (makeCurveEnv)
  curveList <- curveParser(substitute(curveType))
  #curveList <- curveParser(quote(doubleGauss(concave = TRUE)))
  #curveList <- curveParser(quote(logistic()))

  ## Set key
  #setkeyv(dat, c("subject", group))

 # #  ## if(.platform$OStype == windows)
 # #  ## This allows output to be print to console, possibly not possible in windows
 # # #cl <- makePSOCKcluster(4, outfile = "") # prefer to not have output here
 cl <- makePSOCKcluster(4)
 ## Ideally, I could clusterEvalQ bdots
 clusterExport(cl, c("curveFitter", "estDgaussCurve", "dgaussPars",
                     "estLogisticCurve", "logisticPars", "makeCurveEnv"), envir = parent.frame())
 #clusterEvalQ(cl, library(data.table))
 invisible(clusterEvalQ(cl, library(nlme)))
 #invisible(clusterEvalQ(cl, library(bdots))) # someday

 newdat <- split(dat, by = c("subject", group), drop = TRUE) # these needs to not be in string
 ## Would like to pass names into this
 #system.time(res <- parLapply(cl, newdat, bdotsFitter, curveType = "doubleGauss", concave = TRUE, rho = rho, refits = refits, verbose = verbose)) # this also needs arguments (does it?)
 system.time(res <- parLapply(cl, newdat, bdotsFitter,
                              curveList = curveList,
                              rho = rho, refits = refits,
                              verbose = verbose))
 ## This allows us to pass names for verbose
 #system.time(res <- clusterMap(cl, bdotsFitter, newdat, thenames = names(newdat))) # this also needs arguments (does it?)
 stopCluster(cl)

  #### TEST #####
  # curveEnv <- curveParser2(substitute(curveType))
  # cl <- makePSOCKcluster(4)
  # # if rho, refits, verbose are exported, I don't need to pass them to function
  # clusterExport(cl, c("rho", "refits", "verbose",
  #                     "curveFitter", "estDgaussCurve", "dgaussPars",
  #                     "estLogisticCurve", "logisticPars"), envir = parent.frame())
  # invisible(clusterEvalQ(cl, library(nlme)))
  # invisible(clusterExport(cl, ls(envir = curveEnv), envir = curveEnv))
  #
  # newdat <- split(dat, by = c("subject", group), drop = TRUE) # these needs to not be in string
  # ## Would like to pass names into this
  # system.time(res <- parLapply(cl, newdat, bdotsFitter)) # this also needs arguments (does it?)
  # ## This allows us to pass names for verbose
  # #system.time(res <- clusterMap(cl, bdotsFitter, newdat, thenames = names(newdat))) # this also needs arguments (does it?)
  # stopCluster(cl)

  ##### END TEST

 ## Some fits may be null, so we want to make sure we
 ## initialize their data.table below correctly
 # First, let's find a non-null fit
 #tt <- which(vapply(replicate(10, NULL), function(x) !is.null(x), logical(1))) # <- length 0 vector
 nnfit_v <- which(vapply(res, function(x) !is.null(x$fit), logical(1)))
 if(!length(nnfit_v)) stop("No models successfully fit") # here, which should still return something, but stop for now

 ## Null fit template for parameters from first fully fit model
 mat_null <- t(coef(res[[nnfit_v[1]]]$fit))
 mat_null[] <- NA


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
 ## So, to-do here:
 # remove matrix as attachment to it
 # consider putting it in attributes? This would
 # assist with subsetting and letting them adjust it
 # manually. We will put more thought in this later
 fitList <- lapply(names(newdat), function(x) {
   result <- res[[x]] # list of length 3
   x <- strsplit(x, "\\.") # list of by variables for newdat

   dat1 <- as.list(matrix(x[[1]], ncol = length(x[[1]])))
   names(dat1) <-  c(subject, group)


   #set(dat1, j = "fit", value = result[['fit']])
   #dat1$fit <- I(list(result['fit']))
   dat1$R2 <- result[['R2']]
   dat1$AR1 <- (result[['fitCode']] < 3)
   dat1$fitCode <- result[['fitCode']]
   #set(dat1, i = 1, j = "time", time)
   if(!is.null(result[['fit']])) {
     dat1$coef <- list(t(coef(result[['fit']])))
   } else {
     dat1$coef <- list(mat_null)
   }

   # as.data.table works, but do I REALLY want to
   # make these a part of my output? They are visual noise that
   # could be easily fixed with coef(bdotsObj)
   test <- structure(dat1, class = c("data.table", "data.frame"),
                     row.names = integer(0),
                     names = colnames(dat1))
   ## I bet there is a slick way to do this
   test2 <- data.table(colnames(dat1) = dat1)

   # Not sexy, but these aren't large things being copied either
   # Otherwise, need to do more complicated manuevering to determine
   # how to make par columns for NULL fits to be rbindlisted below
   cbind(dat1, dt_par)
 })
 fitList <- rbindlist(fitList)
 ff <- res[[1]][['ff']]

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
 # at some point, just consider R4 objects
 res <- structure(class = c("bdotsObj", "data.table", "data.frame"),
                  .Data = fitList,
                  formula = ff,
                  curveType = curveList,
                  call = match.call(),
                  npars = ncol(dt_null))

}






























