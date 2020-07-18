## This is the parent function file, but for now, we are going to focus on a single subject
# what will likely happen (or maybe makes sense), in parent function, do
# 1) subset based on diff(y, group1(n1, n2)) ~ group2(m1, m2) + group3(w1) (or variant)
# so that at most, each subject may have 2 curves, and we then look between groups. At any rate
# the subject level boot function will determine PAIR and DIFFS but what is present in the
# subsetted dataset
# 2) (these are train of thought) - bdotsBoot(bdotObj, ) { datSub <- bdotObj[sub = 1, ]; bdotBootSub(datSub)}
# 3) here, first, then, is the bdotBootSub (but better named), bdotsBooter?

#### This is FALSE ####
#### But I'll leave it so I know why ####
## if sdat has 4 rows, its diffs and pair
# if sdat has 2 rows, have to check diffGroup/fitGroup
# if one row, nothing to be done, just boot that bitch
#############################################
## In truth, we could have background noise/no background noise (within subject)
# then response type cohort/rhyme
# Ah, but it doesn't make sense to have that be the case and also have subjects be learning disabled, not disabled
# will have to think on this more
## Any group subsets that are length one (diff(y, cond(M, W)) ~ g1(n1, n2) + g2(m1)) will be subset
# before being passed on, so that g2(m2), for example, won't be present

## This is what formula to bdotsBoot might look like
#bdotsParser(ff <- diffs(y, condition(M,W)) ~ group(DLD, TD) + g2(n1))
# Above formula means - within group(DLD) and group(TD), find condition(W) - condition(M).
# call this difference DLD_condDiff and TD_condDiff. Then, take DLD_condDiff - TD_condDiff <- difference of differences
# if DIFFS = FALSE, it would just be y ~ group(DLD, TD), meaning bootstrap DLD and TD, then take differences


## Single argument here, expected to have
# k rows, for the k permutations of groups, w/ unique subject type

## Shit, yo, the results of this can just be appended onto the bdotsObj in the associated rows

## Outline
# bdotsBoot
# 1. parse formula and make appropriate subsets to bdotsOjb
# 2. All we need to retain are `compareGroup` and `diffGroup`, if not NULL
# 3. If diffGroup, we next need to decide if drawing from multivariate
#  i. split bdotsObj into diffgroup1 and diffgroup2
#  ii. See if subject identifiers are unique or repeated
#  iii. If unique, bootstraps not correlated, but yes paired t-test
#  iv. If repeated, bootstraps ARE correlated (bivariate normal), no pair t-test
# 4. If diffGroup NULL, check if IDs are repeated
#   i. If yes, draw bivariate normal, paired t-test
#   ii. If no, standard draw, no paired t-test
# Get p-value adjustments
# find significant regions
# finish

## bdotsBooter
# This will have at most two curves, for when things need
# to be bivariate normal. It will also have an argument for
# correlation coefficient. It does not need to know diffGroup or fitGroup
# it only needs to know correlated or not



# diff group can be LookType
# Note that splitting removes bdotsObj class
# that's OK here - bdotsBooter will use gnls object for values
reslist <- split(res2, by = "LookType")
rl1 <- reslist[[1]]
mm <- matrix(NA, nrow = nrow(rl1), ncol = length(coef(rl1[1, ]$fit[[1]])))
for(i in seq_along(1:nrow(mm))) { # oops this is for all groups
  mm[i, ] <- coef(rl1[i, ]$fit[[1]])
}

mg1 <- mm[c(TRUE, FALSE), ]
mg2 <- mm[c(FALSE, TRUE), ]
corMat <- cor(mg1, mg2)

## This is what a typical entry may look like (this one paired)
dat <- reslist[[1]][Subject == 1, ]
class(dat) <- c("bdotsObj", "data.table", "data.frame")
## Extract coef from  bdotsObj
## uh, this doesn't adress null
# Ah, mother fucker, that's fitcode 6!
coef.bdots <- function(dat) {
  if (!inherits(dat, "bdotsObj")) stop('need bdotsObj')
  nnfit_v <- which(vapply(dat$fit, function(x) !is.null(x$fit), logical(1)))
  if (!length(nnfit_v)) stop("No models contain valid coefficients")
  mm <- matrix(NA, nrow = nrow(dat), ncol = length(cc <- coef(dat[nnfit_v[1], ]$fit[[1]])))
  colnames(mm) <- names(cc)
  for (i in seq_along(1:nrow(mm))) {
    if (dat[i, ]$fitCode != 6) mm[i, ] <- coef(dat[i, ]$fit[[1]])
  }
  mm
}

getVarMat <- function(dat) {
  if(nrow(dat) != 1) stop("only for single row of bdotsObj")
  dat$fit[[1]]$varBeta
}

## We can assume that the paired situation has been figured out previously
# here, all we need to do are fit the curves

## The bivariate normal matrix is going to always be problematic
# because var for base1, base2, ht are nearly 0 (widely different scales than mu, sig1, sig2)
# just look at kappa(sig11)
# for now, I will leave it. I  will try to come up with my own solution before
# giving this to jake
## very much same issue with logistic (cross is HUGE)


bdotsBooter <- function(dat, N.iter, corMat = NULL) {

  ## for now
  if (nrow(dat) > 2) stop("something weird in bdotsBooter")

  mm <- coef.bdots(dat)
  ## Can only be one or two (will this be always true?)
  if (!is.null(corMat)) {
    sig11 <- getVarMat(dat[1, ])
    sig22 <- getVarMat(dat[2, ])
    sig12 <- 0 * corMat %*% sqrt((diag(sig11)) %*% t(diag(sig22)))
    sig <- cbind(rbind(sig11, sig12), rbind(t(sig12), sig22))
    # ee <- eigen(sig)$values
    # ll <- min(abs(ee))/max(abs(ee))
    # tt <- nearPD(sig, keepDiag = TRUE, eig.tol = ll*1e-2, maxit = 1e7)
    pars <- rmvnorm(N.iter, mean = c(t(mm)), sigma = sig)
  } else {
    sig <- getVarMat(dat)
    mm <- coef.bdots(dat)
    pars <- rmvnorm(N.iter, mm, sigma = sig)
  }
  colnames(pars) <- rep(colnames(mm), ncol(pars)/ncol(mm))
  pars
}

## Another point to note somewhere
# a great use of magrittr is at the console
# you generate an output, think hmm, I wanna do
# something with that. Just up in the keyboard, then %>% f()



## Outline
# bdotsBoot
# 1. parse formula and make appropriate subsets to bdotsOjb
# 2. All we need to retain are `compareGroup` and `diffGroup`, if not NULL
# 3. If diffGroup, we next need to decide if drawing from multivariate
#  i. split bdotsObj into diffgroup1 and diffgroup2
#  ii. See if subject identifiers are unique or repeated
#  iii. If unique, bootstraps not correlated, but yes paired t-test
#  iv. If repeated, bootstraps ARE correlated (bivariate normal), no pair t-test
# 4. If diffGroup NULL, check if IDs are repeated
#   i. If yes, draw bivariate normal, paired t-test
#   ii. If no, standard draw, no paired t-test
# Get p-value adjustments
# find significant regions
# finish

### Somewhere, I need to verify order is the same for all of these once they are
# split up, i.e., split by Group, want Cohort always in first row, unrelated in second

## Once this is done, I need to confirm with bob that i have diffs right (fuck this is hard)
# one thing i could do is run the analysis and verify it matches his with what i have
## yo, he goes over this in Cognition paper. Nice
## could make an additional "tests" arguments that could find spots, time, or param
## I need to work on the wording for this/names to not be confused
bdObj <- res2 # from bdotsFit_test.R
bdObj <- res.b; diffGroup = "Group"; compareGroup <- "TrialType"   # bob's data
bdotsBoot <- function(ff, bdObj, N.iter = 1000, alpha = 0.05, p.adj = "oleson", cores = 0, ...) {

  if (cores < 1) cores <- detectCores()/2

  ## pretend ff is parsed. This gives us diffGroup, compareGroup
  diffGroup <- "LookType" # or NULL
  compareGroup <- "Group"

  #, ah, NOPE!
  ## Split by subject and diffgroup. Each resulting
  ## diffGroup can be null, and this splits correctly
  #test <- split(bdObj, by = c("Subject", NULL), drop = TRUE)
  diffList <- split(bdObj, by = c("Subject", diffGroup), drop = TRUE)
  diffList <- lapply(diffList, function(x) {class(x) <- class(bdObj); x})

  ##  First split by diffGroup (could be NULL). Then get corMat if necessary
  ## mm....should probably determine if this is paired now
  # which, actually would be true if each item if diffList is nrow2
  ## Determine if paired
  ## Should not need to error check this. If every split for diffList is length
  # 2, then paired. cor would only fail, then, if number of pars were different
  is.paired <- !any(vapply(diffList, nrow, numeric(1)) != 2)
  if (is.paired) {
    coefs <- lapply(split(bdObj, by = compareGroup, drop = TRUE), function(x) {
      class(x) <- class(bdObj)
      coef.bdots(x)
    })
    corMat <- do.call(cor, setNames(coefs, c("x", "y")))
  }


  ## Draw the bootstraps
  cl <- makePSOCKcluster(4)
  invisible(clusterEvalQ(cl, library(mvtnorm)))
  clusterExport(cl, c("getVarMat", "coef.bdots"), envir = parent.frame()) # <- used in bdotsBooter, but probably needs fixed
  system.time(res <- parLapply(cl, diffList, bdotsBooter, N.iter, corMat))
  stopCluster(cl)

  meanMat <- Reduce(`+`, res)/length(res)

  ## Get formula and turn into function
  ff <- attr(bdObj, "formula")
  f_bod <- deparse1(ff[[3]])
  f_args <- paste0(colnames(coef.bdots(bdObj)), collapse = ", ") # + colnames(dat)
  eval(parse(text = paste('curveFun <- function(', f_args, ', time', ') { return(' , f_bod , ')}', sep='')))

  ## Not the correct way to do this
  if (is.paired) {
    ## Here is where I need the `time` part
    # oh. I need each of these columns to go in seperately
    mm <- as.data.table(meanMat[, 1:4])
    time <- unique(currdata$Time)

    ## Basically, process looks like this
    # I need each column of meanMat to be a list, with additional elelemnt time
    # storing time in a list without it being a list is tricky
    mmList <- split(mm, by = colnames(mm))
    mmList <- lapply(mmList, function(x) {x <- as.list(x); x$time <- time; x}) #!# this solves issue in bdotsObj!
    ## Take away here is this - to make data.table/frame with a list/vector element,
    # start with a list for each row
    # append it to list AsIs
    # rbindlist

    ## This gives me a list of fitted curves
    ## Should really synthesize this with the split/apply/combine paper
    ## also, dude, the power of lists and lazy evaluation. This shit is so cool
    cl <- makePSOCKcluster(4)
    clusterExport(cl, "curveFun")
    res1 <- parLapply(cl, mmList, function(x) do.call(curveFun, x))
    stopCluster(cl)

    ## These didn't work, but illustrating example of why not, given above
    # mmListprep <- lapply(mmList, function(x) {x$time <- unlist(x$time); x})
    # res <- lapply(mmListprep, function(x) do.call(curveFun, as.list(x))) # <- issue here with time being a list

    ## Ok, let's do it again for 2, but obviously this is not how it will be
    ## its still embarassing to get the idea out. But why? I think this is ok as a
    # temporary solution. And it's a MUCH better way of keeping notes of ideas, really
    # This kind of part of my writing step

    ## Anyways, here's N.iter curves for the other group
    mm <- as.data.table(meanMat[, 5:8])
    mmList <- split(mm, by = colnames(mm))
    mmList <- lapply(mmList, function(x) {x <- as.list(x); x$time <- time; x})
    cl <- makePSOCKcluster(4)
    clusterExport(cl, "curveFun")
    res2 <- parLapply(cl, mmList, function(x) do.call(curveFun, x))
    stopCluster(cl)

    ### This is, in general, the idea of differencing curves
    ### so we could wrap this up into a single idea as well
    ## This is paired, so next steps
    # i. take mean of each element in the list
    # ii. subtract them
    #   a. these had to be matrices, not vectors first. False. Use `+` not sum (good to know)
    res11 <- Reduce(`+`, res1) / length(res1)
    res22 <- Reduce(`+`, res2) / length(res2)

    res111 <- Reduce(rbind, res1)
    res222 <- Reduce(rbind, res2)

  }

}































