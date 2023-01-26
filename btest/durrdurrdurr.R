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
load_all("~/packages/bdots/")

# boot <- bdotsBoot(formula = fixations ~ group(A, B),
#                   bdObj = fit, permutation = TRUE)
# # paired test 1
# boot2 <- bdotsBoot(formula = diffs(y, group2(1,2)) ~ group(A, B),
#                   bdObj = fit, permutation = FALSE)
# # paired test 2
# boot3 <- bdotsBoot(formula = diffs(y, group(A, B)) ~ group2(1, 2),
#                   bdObj = fit, permutation = FALSE)
#
# fit$id <- rep(1:20, times = 2)
# # unpaired test 1
# boot4 <- bdotsBoot(formula = diffs(y, group2(1,2)) ~ group(A, B),
#                    bdObj = fit, permutation = FALSE)
# # unpaired test 2
# boot5 <- bdotsBoot(formula = diffs(y, group(A, B)) ~ group2(1, 2),
#                    bdObj = fit, permutation = FALSE)
# # unpaired test1



x <- splitGroups
prs


## Step 1. What kind of pair we got
ps <- isDODpaired(x, prs)
P <- 10

## Step 2. Combine into paired up inner groups
xx <- rbindlist(x)
xx <- split(xx, by = prs[["outerDiff"]])
group <- prs[['innerDiff']]

x <- xx[[1]]

getInnerPermuteMean <- function(x, idx, group) {
  newvec <-  x[[group]][idx]
  set(x, j = group, value = newvec)

  timeName <- attr(x, "call")$time
  TIME <- attributes(x)$time
  ff <- makeCurveFun(x)

  fit_s <- split(x, by = group)
  mvl <- lapply(fit_s, function(y) {
    cc <- coef(y)
    cl <- apply(cc, 1, function(z) {
      z <- as.list(z)
      z[[timeName]] <- TIME
      do.call(ff, z)
    })
    rowMeans(cl)
  })
  dm <- Reduce(`-`, mvl) # this won't work if they are different sizes so rowmeans first
  dm
}

## Step 3. Get indices for permutation, idxa and idxb
  # a. if paired in group find half permute then !flip (since boolean)
  # b. if paired in outer make idxa == idxb
  # c. if not paired at all generate four random idx based on n subjects
## I will do this by creating a function that gives inner diff index
# if (b) just use once, it will return conditional on (a) or (b)
innerIndex <- function(P, n, pair) {
  if (pair) {
    permmat <- replicate(P, sample(c(TRUE, FALSE), n/2, replace = TRUE))
    permmat <- rbind(permmat, !permmat)
    permmat <- apply(permmat, 2, bool2idx)
  } else {
    permmat <- replicate(P, sample(seq_len(n)))
  }
  permmat
}

x <- splitGroups
ps <- isDODpaired(x, prs)
xx <- rbindlist(x)
xx <- split(xx, by = prs[["outerDiff"]])
nv <- vapply(xx, nrow, numeric(1))
P <- 20

## if both paired
if (sum(ps) == 2) {
  permA <- innerIndex(P, nv[1], pair = TRUE)
  permB <- permA
} else if (sum(ps) == 0) {
  ## Neither paired
  permA <- innerIndex(P, nv[1], pair = FALSE)
  permB <- innerIndex(P, nv[2], pair = FALSE)
} else if (ps['ip']) {
  # both paired but outer different
  permA <- innerIndex(P, nv[1], pair = TRUE)
  permB <- innerIndex(P, nv[2], pair = TRUE)
} else if (ps['op']) {
  permA <- innerIndex(P, nv[1], pair = FALSE)
  permB <- permA
}

## Step 4. Compute the inner mean vectors based on this
# I guess do this in two steps at first
clusterExport(cl, c("getInnerPermuteMean", "group", "xx"))
clusterEvalQ(cl, library(bdots))
clusterEvalQ(cl, devtools::load_all("~/packages/bdots"))

## Could also do this with Map but does not lend itself to parallelization like this does
diffmatA <- parApply(cl, permA, 2, function(y) {
  getInnerPermuteMean(xx[[1]], y, group)
})
diffmatB <- parApply(cl, permB, 2, function(y) {
  getInnerPermuteMean(xx[[2]], y, group)
})
# I now need a t statisitic for each row which means i also need some sort of sd
diffmat <- diffmatA - diffmatB

















































