library(bdots)
library(ggplot2)

ci <- as.data.table(ci)
ci <- ci[LookType == "Target", ]

## To get average of original data using data.table syntax
# for `by` variable, include time as well as all grouping columns
ci_mean <- ci[, mean(Fixations), by = c("Time", "protocol")]
head(ci_mean)
ggplot(ci_mean, aes(x = Time, y = V1, color = protocol)) +
  geom_line() + ylab("Mean Fixation") +
  ggtitle("Mean Fixation by Protocol")


## Write bdots output (from email)

res <- bdotsFit(data = ci,
                subject = "Subject",
                time = "Time",
                y = "Fixations",
                group = "protocol",
                curveType = logistic(),
                cor = TRUE,
                numRefits = 2)

boot <- bdotsBoot(formula = y ~ protocol(CI, NH),
                  bdObj = res,
                  Niter = 1000,
                  alpha = 0.05,
                  padj = "oleson",
                  cores = 4)


bdotsresults <- boot

writeBoot <- function(bootObj) {
  timeName <- attr(bootObj, "bdObjAttr")$call$time
  cc <- bootObj$curveList

  ## Create matrix with fit values and sd
  makeMat <- function(ll) {
    fitList <- lapply(names(ll), function(x) {
      mm <- matrix(nrow = length(ll[[x]]$fit), ncol = 2)
      mm[, 1] <- ll[[x]]$fit
      mm[, 2] <- ll[[x]]$sd
      colnames(mm) <- c(x, paste0(x, "_CI"))
      mm
    })
    mat <- Reduce(cbind, fitList)
    time <- matrix(attr(bootObj, "bdObjAttr")$time, dimnames = list(NULL, timeName))
    mat <- cbind(time, mat)
  }

  ## Add indicator column for significance
  addSig <- function(x) {
    st <- bootObj$sigTime
    time <- attr(bootObj, "bdObjAttr")$time
    idx <- lapply(split(st, row(st)), function(y) {
      between(time, y[1], y[2])
    })
    idx <- Reduce(`+`, idx)
    nn <- colnames(x)
    x <- cbind(x, idx)
    colnames(x) <- c(nn, "sig")
    x
  }

  if (inherits(cc, "innerGroupCurveList")) {
    mm <- makeMat(cc)
    mm <- addSig(mm)
  } else {
    diff <- makeMat(cc['diff'])
    cc <- cc[names(cc) != "diff"]
    ccl <- lapply(cc, makeMat)
    ccl <- lapply(names(ccl), function(x) {
      nn <- colnames(ccl[[x]])[-1]
      colnames(ccl[[x]]) <- c(timeName, paste0(x, "__", nn))
      ccl[[x]]
    })
    names(ccl) <- names(cc)
    dmats <- cbind(ccl[[1]][, 6:7], ccl[[2]][, 6:7])
    dmats <- cbind(diff[, 1], dmats, diff[, -1])
    colnames(dmats)[1] <- timeName
    dmats <- addSig(dmats)
    mm <- c("diffs" = list(dmats), ccl)
  }
  return(mm)
}

bootObj <- boot2

#### Try for double gauss
fit2 <- bdotsFit(data = cohort_unrelated,
                 subject = "Subject",
                 time = "Time",
                 y = "Fixations",
                 group = c("Group", "LookType"),
                 curveType = doubleGauss2(concave = TRUE),
                 cor = TRUE,
                 numRefits = 2,
                 cores = 8,
                 verbose = FALSE)

## Don't bother refitting, just to test
fit2 <- bdRemove(fit2, fitCode = 4)

## Difference of difference
boot2 <- bdotsBoot(formula = diffs(Fixations, LookType(Cohort, Unrelated_Cohort)) ~ Group(50, 65),
                        bdObj = fit2,
                        Niter = 1000,
                        alpha = 0.05,
                        padj = "oleson",
                        cores = 8)
