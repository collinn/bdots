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


x <- splitGroups
prs



library(eyetrackSim)


tt <- createData(ar1 = TRUE, trials = 10)$dts
tt2 <- createData(ar1 = FALSE, trials = 100)$dts
tt3 <- createData(ar1 = TRUE, trials = 100)$dts

plot(tt[id == 1, fixations], type = 'l')
plot(tt2[id == 1, fixations], type = 'l')
plot(tt3[id == 1, fixations], type = 'l')

tt <- createData(ar1 = TRUE, trials = 10, manymeans = FALSE, paired = TRUE)
tt2 <- createData(ar1 = FALSE, trials = 100, manymeans = FALSE, paired = TRUE)$dts
tt3 <- createData(ar1 = TRUE, trials = 100, manymeans = FALSE, paired = TRUE)$dts


# if we are not doing many means, i.e., only single means (gay)
if (FALSE) {

}


singleMeans <- function(n, trials, pars, paired, ar1, time) {
  pars <- pars[[1]]
  sigv <- 0.25 / sqrt(trials)

  group1 <- createSingleMeanSubs(n, ar1, pars = pars, sig = sigv, rho = 0.8,
                                 trials, time)

  if (paired) {
    # This just adds gaussian noise to first group
    group2 <- createSingleMeanSubs(n, ar1, pars, rho = 0, sig = sigv, gg = "B",
                                   trials, time)
  } else {
    group2 <- createSingleMeanSubs(n, ar1, pars, rho = 0.8, sig = sigv, gg = "B",
                                   trials, time)
  }
  dts <- rbindlist(list(group1, group2))
  parsA <- matrix(pars, ncol = 4, nrow = n, byrow = TRUE)
  return(list(dts = dts, parsA = parsA, parsB = parsA))
}





tt <- tt$dts

tt <- createData(ar1 = TRUE, trials = 100, manymeans = FALSE)$dts
tt2 <- createData(ar1 = FALSE, trials = 100, manymeans = FALSE)$dts

# tt <- tt[id %in% 1:4, ]
# tt2 <- tt2[id %in% 1:4, ]


colMeans(tt3$parsA)
colMeans(coef(fit1)) # cor true should be true
colMeans(coef(fit2)) # cor FALSE should be true
colMeans(coef(fit3)) # cor true should be false
colMeans(coef(fit4)) # cor FALSE # should be false


fit1 <- bdotsFit(data = tt,
                 subject = "id",
                 group = "group",
                 y = "fixations",
                 time = "time",
                 curveType = logistic())


fit2 <- bdotsFit(data = tt,
                 subject = "id",
                 group = "group",
                 y = "fixations",
                 time = "time",
                 cor = FALSE,
                 curveType = logistic())

fit3 <- bdotsFit(data = tt2,
                 subject = "id",
                 group = "group",
                 y = "fixations",
                 time = "time",
                 curveType = logistic())


fit4 <- bdotsFit(data = tt2,
                 subject = "id",
                 group = "group",
                 y = "fixations",
                 time = "time",
                 cor = FALSE,
                 curveType = logistic())























