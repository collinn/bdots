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
bdotsBoot <- function(ff, bdObj, N.iter = 1000, alpha = 0.05, p.adj = "oleson", cores = 0, ...) {

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
  ########## Everything above is assumed to be working ###########
  ################################################################
  ################################################################


  ## Let's just write out as a bunch of conditionals and then collapse later
  # we can assume subjects are sorted from bdotsFit  (probably true if split is ordered by name)
  if(!is.null(outerDiff)) {

    ## This splits by diffGroup, isPaired defined in helper
    outerDiffList <- split(bdObj, by = outerDiff, drop = TRUE)

    #is.paired <- isPaired(diffList)

    ## Determine if parameter bootstrap should be correlated
    corMat <- lapply(outerDiffList, function(x) {
      tt <- split(x, by = innerDiff, drop = TRUE)
      if (isPaired(tt)) {
        cm <- lapply(tt, coef.bdots)
        do.call(cor, setNames(cm, c("x", "y")))
      } else {
        NULL
      }
    })

    ## We are now going to make outerDiffList a list of lists by subject
    outerDiffList <- lapply(outerDiffList, split, by = "Subject")

    ## Each diffGroup of diffList will now have its subjects fit (go for parLapply since it's two larger groups)
    curveList <- lapply(names(outerDiffList), function(nn) {
      ## Bootstrapped parameters for each subject
      bootPars <- lapply(outerDiffList[[nn]], bdotsBooter, N.iter, corMat[[nn]])
      meanMat <- Reduce(`+`, bootPars)/length(bootPars) # N.iter x numPars (mean across subjects for each iteration)

      ## Here, we are asking - did we fit 4 or 8  parameters
      if (!is.null(corMat[[nn]])) { # fit 8
        mm <- meanMat[, 1:(ncol(meanMat)/2)]
        mmList <- split(mm, row(mm))

        parNames <- colnames(mm)
        mmList <- lapply(mmList, function(x) {
          x <- as.list(x)
          x$time <- time
          setNames(x, c(parNames, "time"))
          })

        ## I should unlist this to get numeric matrix, and then rebuild it
        res1 <- lapply(mmList, function(x){force(x); do.call(curveFun, x)}) # this is 1000 curve fits
        res1 <- matrix(unlist(res1, use.names = FALSE), nrow = length(res1), byrow = TRUE)
        # curve1 <- colMeans(res1) # re bobs data, this is curve for TrialType = M
        # sd1 <- apply(res1, 2, sd)


        mm <- meanMat[, (ncol(meanMat)/2 + 1):ncol(meanMat)]
        mmList <- split(mm, row(mm))
        mmList <- lapply(mmList, function(x) {
          x <- as.list(x)
          x$time <- time
          setNames(x, c(parNames, "time"))
          })
        res2 <- lapply(mmList, function(x){force(x); do.call(curveFun, x)})
        res2 <- matrix(unlist(res2, use.names = FALSE), nrow = length(res2), byrow = TRUE)
        # curve2 <- colMeans(res2) # re bobs data, this is curve for TrialType = W
        # sd2 <- apply(res2, 2, sd)
        list(bootCurve1 = res1, bootCurve2 = res2)
        } else {
          mmList <- split(meanMat, row(meanMat))
          parNames <- colnames(meanMat)
          mmList <- lapply(mmList, function(x) {
            x <- as.list(x)
            x$time <- time
            setNames(x, c(parNames, "time"))
            })
          res <- lapply(mmList, function(x) {force(x); do.call(curveFun, x)})
          res <- matrix(unlist(res, use.names = FALSE), nrow = length(res), byrow = TRUE)
          #curve <- colMeans(res)
          #sdd <- apply(res, 2, sd)
          list(bootCurve1 = res)
          }
      })


}































