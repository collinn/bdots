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
reslist <- split(res2, by = "LookType")
dat <- reslist[[1]][Subject == 1, ]
class(dat) <- c("bdotsObj", "data.table", "data.frame") # split doesn't retain class


## Extract coef from  bdotsObj
## uh, this doesn't adress null
# Ah, mother fucker, that's fitcode 6!
#### Can't replace fit[[1]] since it's unnamed length 1 list. Could name it, I guess
coef.bdots <- function(dat) {
  #if (!inherits(dat, "bdotsObj")) stop('need bdotsObj')
  nnfit_v <- which(vapply(dat$fit, function(x) !is.null(x$fit), logical(1))) #dat$fitCode != 6 (change here and somewhere else I remember)
  if (!length(nnfit_v)) stop("No models contain valid coefficients")
  mm <- matrix(NA, nrow = nrow(dat), ncol = length(cc <- coef(dat[nnfit_v[1], ]$fit[[1]])))
  colnames(mm) <- names(cc)
  for (i in seq_along(1:nrow(mm))) {
    if (dat[i, ]$fitCode != 6) mm[i, ] <- coef(dat[i, ]$fit[[1]])
  }
  mm
}

## Same business here
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
    # sig <- Matrix::nearPD(sig, keepDiag = TRUE, eig.tol = ll*1e-2, maxit = 1e7)$mat
    pars <- mvtnorm::rmvnorm(N.iter, mean = c(t(mm)), sigma = sig)
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
bdObj <- res.b; outerDiff = "Group"; innerDiff <- "TrialType"   # bob's data
bdotsBoot <- function(ff, bdObj, N.iter = 1000, alpha = 0.05, padj = "oleson", cores = 0, ...) {

  if (cores < 1) cores <- detectCores()/2
  ## Could maybe list what was removed
  if(any(bdObj$fitCode == 6)) {
    warning("Some observations had NULL gnls fits. These will be removed")
    bdObj <- bdObj[fitCode != 6, ]
  }
  time <- attr(bdObj, "time")

  ## Get formula and turn into function
  # be warry of name this time
  ff <- attr(bdObj, "formula")
  f_bod <- deparse1(ff[[3]]) # would be cool to get f_args from this instead
  f_args <- paste0(colnames(coef.bdots(bdObj)), collapse = ", ") # + colnames(dat)
  eval(parse(text = paste('curveFun <- function(', f_args, ', time', ') { return(' , f_bod , ')}', sep='')))

  ## Ok, we assume that above, we have parsed formula and have correctly subset our data
  outerDiff <- "Group"; innerDiff <- "TrialType" # bob data
  ## If anything would be null, it would be innerDiff

  ################################################################
  ################################################################
  ########## Everything above is assumed to be working ########### (it's not, yet)
  ################################################################
  ################################################################
  bdObj <- res.b # bob's data
  outerDiff <- "Group"; innerDiff <- NULL; bdObj2 <- bdObj[TrialType == "M", ]
  outerDiff <- "TrialType"; innerDiff <- NULL; bdObj2 <- bdObj[Group == "LI", ]
  outerDiff <- "Group"; innerDiff <- "TrialType"; bdObj2 <- bdObj
  innerDiff <- "Group"; outerDiff <- "TrialType"; bdObj2 <- bdObj

  ################################################################

  ## Here, we get bootstraps of the curves
  # diffs computed there now too
  ## diffList here can be put in recursively defined curveList
  # it's in here that I think we should also return the t statistic
  # curveList last element should indicate if paired (which can be decided in non-base case part of curveBooter)
  # if(is.null(innerDiff)) {
  #   curveList <- curveBooter(bdObj2, splitby = outerDiff, N.iter = 1000, curveFun = curveFun) # nice if this is a list with name of Group
  # } else {
  #   diffList <- split(bdObj, by = outerDiff, drop = TRUE)
  #   curveList <- lapply(diffList, curveBooter, splitby = innerDiff, N.iter = 1000, curveFun = curveFun)
  # }

  ## If diff of diff, the innerdiffs will be labeled groupname.diff
  # with the diff of diff just called "diff"
  # if not diff of diff, the single diff list is just called diff
  # so diff is always the key of our analaysis
  curveList <- curveBooter(bdObj2,
                           outerDiff = outerDiff,
                           innerDiff = innerDiff,
                           N.iter = N.iter,
                           curveFun = curveFun)
  length(curveList)
  names(curveList)
  str(curveList$diff)

  ## Compute t statistic
  tval <- curveList[['diff']][['fit']] / curveList[['diff']][['sd']]
  pval <- 2 * (1 - pt(abs(tval), df = curveList[['diff']][['n']]))

  ## pval adjustment
  # (here's where I need to modify p.adjust to make method oleson)
  system.time(rho <- ar1Solver(tval))
  if (TRUE) {
    ## This is what takes a minute to run
    # (it's also not in parallel yet)
    system.time(alphastar <- findModifiedAlpha(rho,
                                               n = length(tval),
                                               df = curveList[['diff']][['n']]))
  }
  system.time(findModifiedAlpha(rho, n = 100, df = 49))
}

Obj <- bdObj2

## But let's re-evaluate that, because it might be more
# sensible to return curve length t, sd length t

## Returns length 2 nested list (outdated)
# 1. curve1
#   i. curveMat
#   ii. parMat
# 2. curve2
#   i. curveMat
#   ii. parMat
## Make make sense to change inner/outer diff
# The 'n' for diff is associated with the t-statistic, not actual count
curveBooter <- function(Obj, outerDiff, innerDiff = NULL, N.iter, curveFun) {

  if(!is.null(innerDiff)) {
    obj <- split(Obj, by = outerDiff, drop = TRUE)
    res <- lapply(obj, curveBooter, outerDiff = innerDiff,
                  N.iter = N.iter, curveFun = curveFun)
    res <- unlist(res, recursive = FALSE)

    idx <- grep("diff", names(res))
    if (length(idx) != 2) stop("something weird in curveBooter. Contact author")

    ## diff of diff (length one list)
    diffList <- Map(function(x, y) {
      Map(function(a, b) {
        a - b
      }, x, y)
    }, res[idx[1]], res[idx[2]])

    diffList <- diffList[[1]]

    ## snap, we can
    if (ip <- isPaired(obj)) {
      diffList$sd <- apply(diffList$curveMat, 2, sd)
      diffList$n <- nrow(obj[[1]]) - 1
    } else {
      diffList$sd <- nopairSD(res[idx])
      diffList$n <- sum(vapply(obj, nrow, numeric(1))) - 2
    }
    diffList$paired <- ip

    ## Maybe first do something to combine these?
    ## specfically, we should only return a single diff

    return(setNames(c(res, list(diffList)), c(names(res), "diff")))
  }
  oP <- split(Obj, by = outerDiff, drop = TRUE)
  if (ip <- isPaired(oP)) {
    cm <- lapply(oP, coef.bdots)
    corMat <- do.call(cor, setNames(cm, c("x", "y")))
  } else {
    corMat <- NULL
  }
  ## Should investigate how these are different
  # oh, on the split
  if (!is.null(corMat)) {
    outDiffL <- split(Obj, by = "Subject", drop = TRUE)
    bootPars <- lapply(outDiffL, bdotsBooter, N.iter, corMat)
    meanMat <- parMatSplit(Reduce(`+`,  bootPars)/length(bootPars))
  } else {
    outDiffL <- lapply(oP, split, by = "Subject")
    meanMat <- lapply(outDiffL, function(x) {
      bootPars <- lapply(x,  bdotsBooter, N.iter, corMat)
      meanMat <- Reduce(`+`,  bootPars)/length(bootPars)
    })
  }
  curveList <- lapply(seq_along(meanMat), function(i) {
    mm <- meanMat[[i]]
    parNames <- colnames(mm)
    mmList <- lapply(split(mm, row(mm)), function(x) {
      x <- as.list(x)
      x$time <- time
      setNames(x, c(parNames, "time"))
    })
    ## Note, use this  to get mean and sd for the curves
    res <- lapply(mmList, function(x) {force(x); do.call(curveFun, x)})
    res <- matrix(unlist(res, use.names = FALSE), nrow = length(res), byrow = TRUE)
    curveFit <- colMeans(res)
    curveSD <- apply(res, 2, sd)
    list(fit = curveFit, sd = curveSD, curveMat = res, parMat = mm, n = nrow(oP[[i]]))
  })

  ## This is a bit sloppy
  nn <- names(curveList[[1]])
  diffList <- structure(vector("list", length = length(nn) + 1),
                        names = c(nn, "paired"))

  ## Could probably do the double map here like above
  # maybe see which is faster?
  tmp <- unlist(curveList, recursive = FALSE, use.names = FALSE)
  for (i in seq_along(nn)) {
    diffList[[i]] <- tmp[[i]] - tmp[[i + length(nn)]]
  }

  if (ip) {
    diffList$sd <- apply(diffList$curveMat, 2, sd)
    diffList$n <- nrow(oP[[1]]) - 1
  } else {
    diffList$sd <- nopairSD(curveList)
    diffList$n <- sum(vapply(oP, nrow, numeric(1))) - 2
  }
  diffList$paired <- ip


  ## Let's return all that above, and do the diff stuff here (wasted computation if not needed, I guess, but it's only matrix calc)
  setNames(c(curveList, list(diffList)), c(unique(Obj[[outerDiff]]), "diff"))
}



























