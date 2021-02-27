
### Functions used in bdotsBoot ###
# 1. bdotsBooter - bootstraps par values
# 2. bdotsGroupSubset - subsets bdotsFitObj based on group values to bdotsBoot

## bdotsBooter
# Takes subset dat with iter and corMat, returns
# More specifically, it takes a bdotsObj
# N.iter x npars matrix of random draws

# Notes
## bdotsBooter
# This will have at most two curves, for when things need
# to be bivariate normal. It will also have an argument for
# correlation coefficient. It does not need to know diffGroup or fitGroup
# it only needs to know correlated or not
## The bivariate normal matrix is going to always be problematic
# because var for base1, base2, ht are nearly 0 (widely different scales than mu, sig1, sig2)
# just look at kappa(sig11)
# for now, I will leave it. I  will try to come up with my own solution before
# giving this to jake
## very much same issue with logistic (cross is HUGE)
#' @import mvtnorm
bdotsBooter <- function(bdo, N.iter, corMat = NULL) {

  ## for now
  if (nrow(bdo) > 2) stop("something weird in bdotsBooter")

  mm <- coef(bdo)
  ## Can only be one or two (will this be always true?)
  if (!is.null(corMat)) {
    sig11 <- getVarMat(bdo[1, ])
    sig22 <- getVarMat(bdo[2, ])
    sig12 <- 0 * corMat %*% sqrt((diag(sig11)) %*% t(diag(sig22)))
    sig <- cbind(rbind(sig11, sig12), rbind(t(sig12), sig22))
    # ee <- eigen(sig)$values
    # ll <- min(abs(ee))/max(abs(ee))
    # sig <- Matrix::nearPD(sig, keepDiag = TRUE, eig.tol = ll*1e-2, maxit = 1e7)$mat
    pars <- mvtnorm::rmvnorm(N.iter, mean = c(t(mm)), sigma = sig)
  } else {
    sig <- getVarMat(bdo)
    mm <- coef(bdo)
    pars <- mvtnorm::rmvnorm(N.iter, mm, sigma = sig)
  }
  colnames(pars) <- rep(colnames(mm), ncol(pars)/ncol(mm))
  pars
}

## Subset bdObj based on groups being compared
bootGroupSubset <- function(l, bdObj) {
  subnames <- l[["subnames"]]
  subargs  <- l[["subargs"]]
  resNames <- l[['resNames']]

  bd <- bdObj[, c(resNames, subnames), with = FALSE]

  ## This will iteratively subset itself
  for(i in seq_along(subnames)) {
    ss_vec <- bd[[subnames[i]]] %in% subargs[[i]]
    bd <- bd[ss_vec, ]  #subsets multiple times
  }

  ## I'm still going to keep order
  bd[, c(resNames[1], subnames, resNames[-1]), with = FALSE]
}



## alphaAdjust
## curveList - returned from curveBooter
## group - either Group name or Group value, i.e., LI = LI.M/LI.W or W = LI.W/TD.W
## For right now, group can only be group name (since I would need to recalculate diff for group value)
alphaAdjust <- function(curveList, p.adj = "oleson", alpha = 0.05, cores, group = NULL) {
  if (is.null(group)) {
    curve <- curveList[['diff']]
  } else {
    idx <- grep(group, names(curveList))
    if (length(idx) == 0) stop("Invalid group name")
    d_idx <- grep(paste0(group, "\\.diff"), names(curveList[[idx]]))
    curve <- curveList[[idx]][[d_idx]]
  }
  tstat <- curve[["fit"]]/curve[['sd']]
  pval <- 2 * (1 - pt(abs(tstat), df = curve[['n']]))

  p.adj <- match.arg(p.adj, c("oleson", stats::p.adjust.methods))
  rho <- ar1Solver(tstat)

  adjpval <- p_adjust(pval, p.adj, length(tstat), alpha,
                      df = curve[['n']], rho = rho, cores = cores)
  alphastar <- attr(adjpval, "alphastar")

  list(pval = pval, adjpval = adjpval, alphastar = alphastar, rho = rho)
}






