
## Additional notes saved at bottom

## Outside function used
## located in bootHelper.R
# bdotsBooter - gets random sample for subject
# getVarMat - return parameter covariance matrix for observation
# coef.bdots - returns coefficient for observations
# curveBooter - monolithic function that fits/generates the curves


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
# bdObj <- res2 # from bdotsFit_test.R
# bdObj <- res.b; outerDiff = "Group"; innerDiff <- "TrialType"   # bob's data
# formula <- diffs(y, TrialType(M,W)) ~ Group(LI, TD)
bdotsBoot <- function(formula, bdObj, N.iter = 1000, alpha = 0.05, padj = "oleson", cores = 0, ...) {

  if (cores < 1) cores <- detectCores()/2
  ## Could maybe list what was removed
  if(any(bdObj$fitCode == 6)) {
    warning("Some observations had NULL gnls fits. These will be removed")
    bdObj <- bdObj[fitCode != 6, ] # WARNING: DO NOT MODIFY THIS OBJECT EVER (create idx in attr to determine which are valid for use)
  }

  ## This can maybe be two functions
  prs <- bootParser(formula, bdObj)
  bdObj <- bootSubset(prs, bdObj)
  innerDiff <- prs[["innerDiff"]] # unnamed character vector
  outerDiff <- prs[["outerDiff"]] # named character vector? should be consistent

  #time <- attr(bdObj, "time")

  ## Get formula and turn into function
  # be warry of name this time
  ff <- attr(bdObj, "formula")
  f_bod <- deparse1(ff[[3]]) # would be cool to get f_args from this instead
  f_args <- paste0(colnames(coef(bdObj)), collapse = ", ") # + colnames(dat)
  eval(parse(text = paste('curveFun <- function(', f_args, ', time', ') { return(' , f_bod , ')}', sep='')))


  # ################################################################
  # ################################################################
  # ########## Everything above is assumed to be working ########### (it's not, yet)
  # ################################################################
  # ################################################################
  # bdObj <- res.b # bob's data
  # outerDiff <- "Group"; innerDiff <- NULL; bdObj2 <- bdObj[TrialType == "M", ]
  # outerDiff <- "TrialType"; innerDiff <- NULL; bdObj2 <- bdObj[Group == "LI", ]
  # outerDiff <- "Group"; innerDiff <- "TrialType"; bdObj2 <- bdObj
  # innerDiff <- "Group"; outerDiff <- "TrialType"; bdObj2 <- bdObj
  #
  # ################################################################


  ## If diff of diff, the innerdiffs will be labeled groupname.diff
  # with the diff of diff just called "diff"
  # if not diff of diff, the single diff list is just called diff
  # so diff is always the key of our analaysis
  curveList <- curveBooter(bdObj,
                           outerDiff = outerDiff,
                           innerDiff = innerDiff,
                           N.iter = N.iter,
                           curveFun = curveFun)
  ip <- curveList[['diff']][['paired']] # paired?
  # length(curveList)
  # names(curveList)
  # str(curveList$diff)

  ## Compute t statistic
  # diff already set to be the correct values, based on what was determined
  # in the curveList function
  ## Abstract below into own function so that it can be used in plots
  tval <- curveList[['diff']][['fit']] / curveList[['diff']][['sd']]
  pval <- 2 * (1 - pt(abs(tval), df = curveList[['diff']][['n']]))

  ## pval adjustment
  # (here's where I need to modify p.adjust to make method oleson)
  rho <- ar1Solver(tval)
  if (TRUE) {
    ## This is what takes a minute to run
    # (35s on Bob's data, not splendid)
    # (it's also not in (windows) parallel yet)
    alphastar <- findModifiedAlpha(rho,
                                   n = length(tval),
                                   df = curveList[['diff']][['n']],
                                   cores = 2)
    k <- alphastar/alpha
    adjpval <- pval/k
  }
  time <- attr(bdObj, "time")
  sigTime <- bucket(pval <= alphastar, time)

  #bdAttr <- attributes(bdObj)
  structure(class = "bdotsBootObj",
            .Data = list(curveList = curveList,
                         alpha = alpha,
                         adj.alpha = alphastar,
                         adj.pval = adjpval,
                         rho = rho,
                         paired = ip,
                         diffs = c(outerDiff, innerDiff)),
            call = match.call(),
            bdObjAttr = attributes(bdObj))

  #### Actually, both of these can be computed from the returned object
  # so let's just make these their own methods for the return value
  # of bdotsBoot

  ## Add time test later

  ## Add param test later

  ## Add plotting/CI functions later

}


### Notes

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











