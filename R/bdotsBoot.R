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

# should finish bdotsFitter so I can replicate this
# taken from current bdotsFit
sdat <- tt2[subject == 1, ]
diffGroup <- c("LookType", "Cohort", "Unrelated_Cohort")
fitGroup <- c("Group", "50", "65")
bdotsBooter <- function(sdat, N.iter, diffGroup = NULL, fitGroup = NULL) {

  ## Get names of Groupvars being kept
  vv <- vapply(list(diffGroup, fitGroup), function(x) x[[1]], character(1))
  sdat[, .N, by = vv]

  ## Perhaps a better way is to determine number of permutations in fitGroup/diffGroup

  ## Probably ought to do better checks than nrow == 1
  # check diffGroup, for instance
  if (nrow(sdat) == 1) {
    fit <- sdat$fit[[1]]$fit
    bootPars <- rmvnorm(N.iter, coef(fit), fit$varBeta)
  }
}
