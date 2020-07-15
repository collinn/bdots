library(data.table)
library(magrittr)
library(parallel)

## currently, only actually using
# data.table
# parallel
# nlme




# pair verbose with message

## If rho null but cor = TRUE, couldn't we guess at rho?

## Need to not make my own names for this (i.e., y, subject, time)

## Add argument for minimally accepted R2. For example, if a
# model fits with AR1 = TRUE, but R2 =


## Yo, if we end up having to delete subjects, how do we do a paired t test?
bdotsFit <- function(data, # dataset
                     subject, # subjects
                     time, # column for time
                     y, # response vector
                     group, # groups for subjects
                     curveType = c("doubleGauss"), # logistic, doubleGauss, etc. maybe have length match responseGroup length?
                     concave = NULL, # doubleGauss only concavity
                     cor = TRUE, # autocorrelation?
                     rho = 0.9, # autocor value
                     refits = 0,
                     cores = 0, # cores to use, 0 == 50% of available
                     verbose = FALSE) {

  if(cores < 1) cores <- detectCores()/2


  ## Cheat around DT reference, conditionally set key for subset
  ## Can possibly avoid this set, as was done below
  # and get rid of magrittr
  # maybe get rid of data table :(
  ## For doubleGauss
  dat <- data.table()
  dat$subject <- data[[subject]]
  dat$time <- data[[time]] %>% as.numeric()
  dat$y <- data[[y]] %>% as.numeric()
  group <- c("Group", "LookType")

  ## Set group variables in data.tableb+
  for(gg in seq_along(group)) {
    set(dat, j = group[gg], value = data[[group[gg]]])
  }

  ## For logistic
  data(ci)
  ci <- as.data.table(ci)
  ci <- ci[LookType == "Target", ]
  group <- "protocol"
  dat <- data.table()
  dat$subject <- ci[[subject]]
  dat$time <- ci[[time]] %>% as.numeric()
  dat$y <- ci[[y]] %>% as.numeric()
  for(gg in seq_along(group)) {
    set(dat, j = group[gg], value = ci[[group[gg]]])
  }



  ## Set key
  setkeyv(dat, c("subject", group))

  ## if(.platform$OStype == windows)
  ## This allows output to be print to console, possibly not possible in windows
 cl <- makePSOCKcluster(4, outfile = "") # prefer to not have output here
 clusterExport(cl, c("curveType", "concave", "cor", "rho", "refits", "verbose",
                     "curveFitter", "estDgaussCurve", "dgaussPars", "gnls", "corAR1",
                     "estLogisticCurve", "logisticPars"), envir = parent.frame())
 #clusterEvalQ(cl, library(data.table))
 invisible(clusterEvalQ(cl, library(nlme)))
 newdat <- split(dat, by = c("subject", group), drop = TRUE) # these needs to not be in string
 ## Would like to pass names into this
 system.time(res <- parLapply(cl, newdat, bdotsFitter)) # this also needs arguments (does it?)
 ## This allows us to pass names for verbose
 #system.time(res <- clusterMap(cl, bdotsFitter, newdat, thenames = names(newdat))) # this also needs arguments (does it?)
 stopCluster(cl)


 ## At some point, need to change the way this returns fit entry as "AsIs" object, nested lis t
 fitList <- lapply(names(newdat), function(x) {
   result <- res[[x]] # list of length 3
   x <- strsplit(x, "\\.") # list of by variables for newdat

   dat1 <- as.data.table(matrix(x[[1]], ncol = length(x[[1]])))
   names(dat1) <-  c(subject, group)

   # I don't really like that this has to be like this
   # it's a length 1 list where result is a list of length 3, but w/e for now
   # since assigning class, could also make method `[` or `[[`
   dat1$fit <- I(list(result['fit'])) #<- this fixes it so that object is gnls inside list of length 1


   # These got named weird things from function, so use number. ew
   # Should also include starting parameters here (only used for refitting or end user)
   # that way, if user just throws in own model, they don't need to include extra nonsense
   # if they do nothing, bdotsRefit will just keep jittering the start values
   dat1$R2 <- result['R2']
   dat1$AR1 <- (result['fitCode'] < 3)
   dat1$fitCode <- result['fitCode']
   dat1
 })
 fitList <- rbindlist(fitList)


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

   res <- structure(class = c("bdotsObj", "data.table", "data.frame"),
                    list(curveType = curveType,
                         formula = res[[1]]$ff,
                         fits = fitList))


}






























