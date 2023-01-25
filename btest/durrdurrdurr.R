library(bdots)
library(eyetrackSim)
library(mvtnorm)

## Start by generating an empirical distribution
ci <- as.data.table(ci)
ci <- ci[LookType == "Target", ]
#ci <- ci[LookType == "Target" & protocol != "NH", ]
fit <- bdotsFit(data = ci,
                y = "Fixations",
                subject = "Subject",
                time = "Time",
                group = "protocol", curveType = logistic())

## Except this time I only need one empirical dist
## Here I am using BOTH TD/ND to create more variability
cc <- coef(fit[])
vv <- var(cc)
cc <- colMeans(cc)
cc[1] <- abs(cc[1])
cc[2] <- pmin(cc[2], 1)
pars <- list(mean = cc, sigma = vv)

## Will always create twice what N is (so two groups by default)
createData <- function(n = 25, trials = 100, pars, paired = FALSE, pairMag = 0.05) {
  time <- seq(0, 2000, by = 4)
  newpars <- do.call(rmvnorm, as.list(c(n, pars)))
  newpars[,1] <- abs(newpars[,1]) # need base > 0
  newpars[,2] <- pmin(newpars[,2], 1) # need peak < 1
  spars <- split(newpars, row(newpars))
  dts1 <- lapply(seq_len(n), function(x) {
    pp <- spars[[x]]
    dt <- data.table(id = x,
                     time = time,
                     group = "A",
                     true = eyetrackSim:::logistic_f(pp, time))
    dt[, fixations := mean(rbinom(trials, 1, true)), by = time]
  })

  dts1 <- rbindlist(dts1)

  ## Then we make our parameters for group 2
  if (!paired) {
    ## Basically just repeat above, exact same distribution
    newpars2 <- do.call(rmvnorm, as.list(c(n, pars)))
    newpars2[,1] <- abs(newpars2[,1]) # need base > 0
    newpars2[,2] <- pmin(newpars2[,2], 1) # need peak < 1
  } else {
    ## Keep the original pars from newpars
    orig_pars <- newpars
    ## Then make one with mean 0
    pars2 <- pars
    pars2$mean[] <- 0
    pars2$sigma <- pars2$sigma*pairMag
    ## This gets the variance
    varpars <- do.call(rmvnorm, as.list(c(n, pars2)))
    ## And then we make our paired parameters
    newpars2 <- orig_pars + varpars
    newpars2[,1] <- abs(newpars2[,1]) # need base > 0
    newpars2[,2] <- pmin(newpars2[,2], 1) # need peak < 1
  }
  spars2 <- split(newpars2, row(newpars2))
  ipn <- ifelse(paired, 0, n)
  # id not correct for paired here
  dts2 <- lapply(seq_len(n), function(x) {
    pp <- spars2[[x]]
    dt <- data.table(id = x + ipn, #ipn is 0 if paired
                     time = time,
                     group = "B",
                     true = eyetrackSim:::logistic_f(pp, time))
    dt[, fixations := mean(rbinom(trials, 1, true)), by = time]
  })
  dts2 <- data.table::rbindlist(dts2)
  dts <- data.table::rbindlist(list(dts1, dts2))

  return(list(dts = dts, parsA = newpars, parsB = newpars2))
}



tt <- createData(n = 10, pars = pars, paired = TRUE)$dts
tt2 <- createData(n = 10, pars = pars, paired = TRUE)$dts

tt$group2 <- "1"
tt2$group2 <- "2"

tt <- rbindlist(list(tt, tt2))

fit <- bdotsFit(data = tt,
                subject = "id",
                time = "time",
                y = "fixations",
                group = c("group","group2"),
                curveType = logistic(),
                cores = 7L)

bdObj <- fit
formula <- diffs(y, group2(1,2)) ~ group(A, B)
#formula <- y ~ group(A, B)
permutation <- FALSE
Niter <- 10
alpha <- 0.05
cores <- 7
#load_all("~/packages/bdots/")

boot <- bdotsBoot(formula = fixations ~ group(A, B),
                  bdObj = fit, permutation = FALSE)
# paired test 1
boot2 <- bdotsBoot(formula = diffs(y, group2(1,2)) ~ group(A, B),
                  bdObj = fit, permutation = FALSE)
# paired test 2
boot3 <- bdotsBoot(formula = diffs(y, group(A, B)) ~ group2(1, 2),
                  bdObj = fit, permutation = FALSE)

fit$id <- rep(1:20, times = 2)
# unpaired test 1
boot4 <- bdotsBoot(formula = diffs(y, group2(1,2)) ~ group(A, B),
                   bdObj = fit, permutation = FALSE)
# unpaired test 2
boot5 <- bdotsBoot(formula = diffs(y, group(A, B)) ~ group2(1, 2),
                   bdObj = fit, permutation = FALSE)
# unpaired test1

#
# x <- splitGroups
# prs <- prs
# alpha <- 0.05
# P <- 250
#
#
# x <- splitGroups
# sgroups <- splitGroups
# b <- Niter
# prs

#ok but what if dod???

#' Create bootstrapped distribution of groups
#'
#' Function to create bootstrapped distribution of groups depending
#' on paired status and difference of difference
#'
#' @param x list of bdObj
#' @param prs list from boot parser
#' @param b number of bootstraps
createGroupDists <- function(x, prs, b) {

  ## Now tend to parallel needs
  if (Sys.info()['sysname'] == "Darwin") {
    cl <- makePSOCKcluster(cores, setup_strategy = "sequential")
  } else {
    cl <- makePSOCKcluster(cores)
  }
  invisible(clusterEvalQ(cl, {library(bdots)}))

  ## Is this a difference of difference?
  dod <- !is.null(prs$innerDiff)

  if (!dod) {
    ## This means only length 2 so checking if paired is simple
    ip <- isPaired(x)
    if (!ip) {
      groupDists <- parLapply(cl, x, getBootDistUnpaired, b) # <- the old way, nice
      #stopCluster(cl)
      #return(groupDists) # done
    } else {
      # if it is paired
      groupDists <- getBootDistPaired(cl, x, b)
      #stopCluster(cl)
      #return(groupDists)
    }
  } else { # if difference of difference
    ## Find out if inner/outer paired
    ps <- isDODpaired(x, prs)

    ## Based on paired status
    if (sum(ps) == 2) {
      ## both paired
      groupDists <- getBootDistPaired(cl, x, b)
    } else if (sum(ps) == 0) {
      ## Neither paired
      groupDists <- parLapply(cl, x, getBootDistUnpaired, b)
    } else {
      ## This is the case if either inner or outer not but other is
      nx <- strsplit(names(x), split = "\\.")

      if (ps['ip']) {
        # paired inner group will share the same outer group in the split object
        ogn <- prs$subargs[[2]] # I hate this but this is values for outerDiff which also doesn't even make sense
        idx <- vapply(nx, function(y) ogn[1] %in% y, logical(1))
      } else if (ps['op']) {
        ign <- prs$subargs[[1]]
        idx <- vapply(nx, function(y) ign[1] %in% y, logical(1))
      }
      gd1 <- getBootDistPaired(cl, x[idx], b)
      gd2 <- getBootDistPaired(cl, x[!idx], b)
      groupDists <- c(gd1, gd2)
    }
  }
  stopCluster(cl)
  return(groupDists)
}































































