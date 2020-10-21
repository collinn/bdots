#' Create bootstrapped curves from bdotsObj
#'
#' Creates bootstrapped curves and performs alpha adjustment. Can perform
#' "difference of difference" for nested comparisons
#'
#' @param formula See details.
#' @param bdObj An object of class 'bdotsObj'
#' @param N.iter Number of iterations of bootstrap to draw
#' @param alpha Significance level
#' @param p.adj Adjustment to make to pvalues for significance. Will be able to
#' use anything from \code{p.adjust} function, but for now, just "oleson"
#' @param cores Number of cores to use in parallel. Default is zero, which
#' uses half of what is available.
#'
#' @details The formula is the only tricky part of this. There will be a minor
#' update to how it works in the future.  The three parts we will examine here
#' are Groups, the LHS, and the RHS. This is incompelete below, but is detailed
#' in the original vignette sent out
#'
#' ## Groups
#'
#' To be confusing, we use "groups" in two senses here. The first use, which we
#' will denote with captial "G" Group, indicates the column name of the dataset
#' representing one of the variables by which we will distinguish categories of
#' observations. Lowercase "g" group will be the values in this column.
#'
#' ## LHS
#'
#' If we are comparing two groups, we will just use the name of the outcome, i.e.,
#' `y ~` .
#'
#' @return  Object of class 'bdotsBootObj'
#'
#' @import data.table
#' @export

bdotsBoot <- function(formula, bdObj, N.iter = 1000, alpha = 0.05, p.adj = "oleson", cores = 0, ...) {

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

  if (cores < 1) cores <- detectCores()/2
  ## Could maybe list what was removed
  ## Also need to remove relevant things from X
  ## looking back at refit function here, we  need a function
  # that takes an expression, returns subject IDS for which that
  # expression is true, and then removes all of that subject
  # (including observations that did not meet that expression)
  if (any(bdObj[['fitCode']] == 6)) {
    warning("Some observations had NULL gnls fits. These and their pairs will be removed")
    bdObj <- bdRemove(bdObj, fitCode = 6, removePairs = TRUE)
  }

  ## This can maybe be two functions
  prs <- bootParser(formula, bdObj)
  bdObj <- bootSubset(prs, bdObj)
  innerDiff <- prs[["innerDiff"]] # unnamed character vector
  outerDiff <- prs[["outerDiff"]] # named character vector? should be consistent

  curveGrps <- setNames(prs[['subargs']], prs[['subnames']])

  #time <- attr(bdObj, "time")

  ## Get formula and turn into function
  # be wary of name this time
  # ff <- attr(bdObj, "formula")
  # f_bod <- deparse1(ff[[3]]) # would be cool to get f_args from this instead
  # f_args <- paste0(colnames(coef(bdObj)), collapse = ", ") # + colnames(dat)
  # eval(parse(text = paste('curveFun <- function(', f_args, ', ',
  #                         attr(bdObj, "call")$time, ') { return(' , f_bod , ')}',
  #                         sep='')))

  ## This replaces the commented out code above, which threw a NOTE
  # on the R CMD CHECK for having undefined global variable curveFun
  ## From adv R, kinda
  makeCurveFun <- function(bdObj) {
    time <- attr(bdObj, "call")[['time']]
    f_bod <- attr(bdObj, "formula")[[3]]
    f_args <- c(colnames(coef(bdObj)), time)
    f_args <- setNames(as.pairlist(rep("", length(f_args))), f_args)
    # f_args <- setNames(as.pairlist(c(f_args, time)), c(f_args, time))
    eval(call("function", f_args, f_bod), parent.frame())
  }


  curveFun <- makeCurveFun(bdObj)


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

  ## make curveList more compact if diff of diff
  if (!is.null(innerDiff)) {
    vals <- curveGrps[[outerDiff]]
    if (length(vals) != 2) stop("Error 987, contact package author")
    idx1 <- grep(paste0("^", vals[1]), names(curveList))
    idx2 <- grep(paste0("^", vals[2]), names(curveList))
    idxDiff <- grep("^diff$", names(curveList))
    curveList <- list((curveList[idx1]),
                 (curveList[idx2]),
                 (curveList[[idxDiff]]))
    names(curveList) <- c(vals, "diff")
    }

  # length(curveList)
  # names(curveList)
  # str(curveList$diff)

  ### Below has been replaced with alphaAdjust in bootHelper.R
  ## Compute t statistic
  # diff already set to be the correct values, based on what was determined
  # in the curveList function
  ## Abstract below into own function so that it can be used in plots
  # tval <- curveList[['diff']][['fit']] / curveList[['diff']][['sd']]
  # pval <- 2 * (1 - pt(abs(tval), df = curveList[['diff']][['n']]))
  #
  # ## pval adjustment
  # # (here's where I need to modify p.adjust to make method oleson)
  # rho <- ar1Solver(tval)
  # if (TRUE) {
  #   ## This is what takes a minute to run
  #   # (35s on Bob's data, not splendid)
  #   # (it's also not in (windows) parallel yet)
  #   alphastar <- findModifiedAlpha(rho,
  #                                  n = length(tval),
  #                                  df = curveList[['diff']][['n']],
  #                                  cores = 2)
  #   k <- alphastar/alpha
  #   adjpval <- pval/k
  # }
  res <- alphaAdjust(curveList, p.adj, alpha, cores)
  pval <- res[['pval']]
  rho <- res[['rho']]
  alphastar <- res[['alphastar']]
  adjpval <- res[['adjpval']]
  time <- attr(bdObj, "time")
  sigTime <- bucket(pval <= alphastar, time)
  dod <- ifelse(is.null(innerDiff), FALSE, TRUE)
  #attr(bdObj, 'X')$X <- NULL

  #bdAttr <- attributes(bdObj)
  structure(class = "bdotsBootObj",
            .Data = list(curveList = curveList,
                         alpha = alpha,
                         adj.alpha = alphastar,
                         adj.pval = adjpval,
                         sigTime = sigTime,
                         rho = rho,
                         paired = ip,
                         diffs = c("outerDiff" = outerDiff,
                                   "innerDiff" = innerDiff),
                         curveGroups = curveGrps,
                         dod = dod,
                         curveFun = curveFun),
            call = match.call(),
            bdObjAttr = attributes(bdObj))

  #### Actually, both of these can be computed from the returned object
  # so let's just make these their own methods for the return value
  # of bdotsBoot

  ## Add time test later

  ## Add param test later

  ## Add plotting/CI functions later


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






}








