library(data.table)
library(magrittr)
library(parallel)

## currently, only actually using
# data.table
# parallel
# nlme
# we could almost remove data.table....
# since we don't have any operations  w/n groups
# and the dat[, list(f...), by = , .SDcols = ] isn't parallelized
# However, it will be nice for subsetting things like
# result[fitCode == xyz, this thing]
# or plot(result[lkajsdfl, alkdsjf])
## Keeping dependencies down is dope. Eff you stringr

## possible additions
# refits - number of times to jitter initial parameters

## for each fit, give numeric value, i.e., 1 - AR1 R2 > 0.95, 2 - AR1 R2 > 0.85, etc.
## this will make refitting SO much easier
## it will make dropping subjects easier too!
## Definitely going to be using data table for this

# pair verbose with message

## If rho null but cor = TRUE, couldn't we guess at rho?

## Need to not make my own names for this (i.e., y, subject, time)

## Add argument for minimally accepted R2. For example, if a
# model fits with AR1 = TRUE, but R2 =

## Go through and replace jitter with something else

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


 ## Unfortunatley, this is faster x3 (and far less sexy xInf)
  # ugh, about same speed with data.frame.
 cl <- makePSOCKcluster(4)
 clusterExport(cl, c("curveType", "concave", "cor", "rho", "refits", "verbose",
                     "curveFitter", "estDgaussCurve", "dgaussPars", "gnls", "corAR1",
                     "estLogisticCurve", "logisticPars"), envir = parent.frame())
 #clusterEvalQ(cl, library(data.table))
 invisible(clusterEvalQ(cl, library(nlme)))
 newdat <- split(dat, by = c("subject", group), drop = TRUE) # these needs to not be in string
 system.time(res <- parLapply(cl, newdat, bdotsFitter)) # this also needs arguments (does it?)

 stopCluster(cl)

 #tt <- strsplit(names(newdat), "\\.")
 # this gives me same as before (good), but still need to
 # do processing on 'result'
 #x <- names(newdat)[1]
 ## At some point, need to change the way this returns fit entry as "AsIs" object, nested lis t
 tt <- lapply(names(newdat), function(x) {
   result <- res[[x]] # list of length 3
   x <- strsplit(x, "\\.") # list of by variables for newdat

   dat1 <- as.data.table(matrix(x[[1]], ncol = length(x[[1]])))
   #dat1 <- as.data.frame(matrix(x[[1]], ncol = length(x[[1]])))
   names(dat1) <-  c(subject, group)

   # I don't really like that this has to be like this
   # it's a length 1 list where result is a list of length 3, but w/e for now
   # since assigning class, could also make method `[` or `[[`
   dat1$fit <- I(list(result[1])) #<- this fixes it so that object is gnls inside list of length 1


   # These got named weird things from function, so use number. ew
   # Should also include starting parameters here (only used for refitting or end user)
   # that way, if user just throws in own model, they don't need to include extra nonsense
   # if they do nothing, bdotsRefit will just keep jittering the start values
   dat1$R2 <- result['R2']
   dat1$AR1 <- (result['fitCode'] < 3)
   dat1$fitCode <- result['fitCode']
   dat1
 })
 tt <- rbindlist(tt)
 rr <- tt[1, ]
 ## This alone is a few seconds (2.097 vs 0.021)
 ## ok, w/e, find having this as data.frame, but it prints out wrong
 #system.time(tt1 <- Reduce(rbind, tt))
 #system.time(tt2 <- rbindlist(tt))

 ## What else should be returned? anything in sub level stuff?
 # perhaps a curve type, so list(curvetype, tt1), for example? what else?
 # think of summary functions

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



  return(dat)

}






























