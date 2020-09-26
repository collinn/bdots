
## Let's test bdotsFit funcitoin

rm(list = ls())
library(bdots)
library(data.table)
library(parallel)
library(nlme)
library(mvtnorm)
library(pryr)
load("~/packages/bdots/data/testRunData.RData")
# Not for parallel, I could do clusterEvalQ(cl, source(".."))
source("~/packages/bdots/R/bdotsFit.R")
source("~/packages/bdots/R/bdotsFitter.R")
source("~/packages/bdots/R/parser.R")
source("~/packages/bdots/R/bootParser.R")
source("~/packages/bdots/R/curveFitter.R")
source("~/packages/bdots/R/doubleGauss.R")
source("~/packages/bdots/R/logistic.R")
source("~/packages/bdots/R/helper.R")
source("~/packages/bdots/R/effectiveAlpha.R")
source("~/packages/bdots/R/findModifiedAlpha.R")
source("~/packages/bdots/R/ar1Solver.R")
source("~/packages/bdots/R/bucket.R")
source("~/packages/bdots/R/bootHelper.R")
source("~/packages/bdots/R/bdotsBoot.R")
source("~/packages/bdots/R/plotFunctions.R")
source("~/packages/bdots/R/summaryFunctions.R")
source("~/packages/bdots/R/bdotsRefit.R")


##############

load("~/bdots/demo.RData")
data <- cohort_unrelated
group = c("Group", "LookType")
time <- "Time"
subject <- "Subject"
y <- "Fixations"
concave <-  TRUE
conc <- TRUE
rho <-  0.9
numRefits <- 0
cores <-  1
verbose <-  FALSE
alpha <- 0.05
curveType <- c("doubleGauss")

#head(data)


## was missing concave smh
res <- bdotsFit(data = cohort_unrelated,
                subject = "Subject",
                time = "Time",
                y = "Fixations",
                group = c("Group", "LookType"),
                curveType = doubleGauss(concave = TRUE),
                cor = TRUE,
                numRefits = 2,
                cores = 0,
                verbose = FALSE)
res2 <- res[Subject %in% c(1, 2, 3, 5, 7:11, 14:21, 23:26)]

boot.test <- bdotsBoot(formula = diffs(Fixations, LookType(Cohort, Unrelated_Cohort)) ~ Group(50, 65),
                                  bdObj = res2,
                                  N.iter = 1000,
                                  alpha = 0.05,
                                  p.adj = "oleson",
                                  cores = 4)
#plotCompare(boot.test)

## Doing this, in boot.test2$curveList$Cohort, there are entries for both cohort and unrelated cohort?
# it's because I grep for "Cohort" which is in "Cohort" AND "Unrelated_Cohort".  Ah, subset of name
boot.test2 <- bdotsBoot(formula = diffs(y, Group(50, 65)) ~ LookType(Cohort, Unrelated_Cohort),
                                   bdObj = res2,
                                   N.iter = 1000,
                                   alpha = 0.05,
                                   p.adj = "oleson",
                                   cores = 4)
## Keep only valid pairs for now
#res2 <- res[, 1:7]

# this removes possible pairs
#res2 <- res2[fitCode != 6, ]




## logistic
library(bdots)
data(ci)
ci <- as.data.table(ci)
ci <- ci[LookType == "Target", ]
res.l <- bdotsFit(data = ci,
                              subject = "Subject",
                              time = "Time",
                              y = "Fixations",
                              group = "protocol",
                              curveType = logistic(),
                              cor = TRUE,
                              numRefits = 2)
boot.l <- bdotsBoot(formula = y ~ protocol(CI, NH),
                                  bdObj = res.l,
                                  N.iter = 1000,
                                  alpha = 0.05,
                                  p.adj = "oleson",
                                  cores = 4)
#plotCompare(boot)

#   ## For logistic
#   data(ci)
#   ci <- as.data.table(ci)
#   ci <- ci[LookType == "Target", ]
#   group <- "protocol"
#   dat <- data.table()
#   dat$subject <- ci[[subject]]
#   dat$time <- ci[[time]] %>% as.numeric()
#   dat$y <- ci[[y]] %>% as.numeric()
#   for(gg in seq_along(group)) {
#     set(dat, j = group[gg], value = ci[[group[gg]]])
#   }


## Bdots for bob's  data

bobdat <- read.delim("~/bdots/mcmurray_folder/bdots_stuff/example bdots code/data.txt")
bobdat <- as.data.table(bobdat)
currdata <- bobdat [bobdat$TrialType == "W" | bobdat$TrialType == "M" , ]
names(currdata)[names(currdata) == 'dx'] <- 'Group'
currdata2 <- as.data.table(currdata)
currdata2 <- currdata2[Subject != 405, ]
currdata2 <- currdata2[Subject != 1699, ]
currdata2 <- currdata2[Subject != 1526, ]
currdata2 <- currdata2[Subject != 1647, ]
currdata2 <- currdata2[Subject != 1688, ]
currdata2 <- currdata2[Subject != 1582, ]
N.iter <- 1000

## About 32 seconds
system.time(res.b <- bdotsFit(data = currdata2,
                              subject = "Subject",
                              time = "Time",
                              y = "Looks",
                              group = c("Group", "TrialType"),
                              curveType = logistic(),
                              cor = TRUE,
                              numRefits = 2))

#debugonce(bdotsFit)

# head(sort(res.b$R2))
#bdObj <- res.b
#tt <- currdata[Subject == 405, ]
#tt1 <- tt[TrialType == "M", ]
#tt2 <- tt[TrialType == "W", ]

###
# actually surprised here - about ~ 18 - 24 seconds
# Ok, but another time it took like k 47 seconds?
# third time 18?
# Another 31?
## Something about the diff function is very wrong
# should investigate nopairSD in helper.R
system.time(boot.res <- bdotsBoot(formula = diffs(Looks, TrialType(M,W)) ~ Group(LI, TD),
                      bdObj = res.b,
                      N.iter = 1000,
                      alpha = 0.05,
                      p.adj = "oleson",
                      cores = 4))

## Without diff, need to make sure we are subsetting correctly
# easy to forget to add TrialType(M) or TrialType(W)
# so be sure to check it
boot.res2 <- bdotsBoot(formula = y ~ Group(LI, TD) + TrialType(M),
                                  bdObj = res.b,
                                  N.iter = 1000,
                                  alpha = 0.05,
                                  p.adj = "oleson",
                                  cores = 4)


save(res2, boot.test, boot.test2, boot.res, res.b, res.l, boot.l,
     file = c("~/packages/bdots/data/testRunData.RData"))


#debugonce(bdotsBoot)
bdBootObj <- boot.res
rr <- boot.res$curveList
rr1 <- rr$diff

plot(boot.res)

plotCompare(boot.res, group = "LI")
plotCompare(boot.res2)
plotCompare(boot.res, diffs = FALSE)

# > boot.res$sigTime
# [,1] [,2]
# [1,]   92 1200
# [2,] 1472 2500


dt <- data.table(x = 9)

rr <- structure(.Data = dt, class = class(dt), someatt = "test")
y <- rr
f <- function(y) {
  attributes(y) <- c(attributes(y), list("newat" = "newtest"))
  y
}
f(rr)

boot.res2 <- bdotsBoot(formula = y ~ Group(LI, TD) + TrialType(M),
                       bdObj = bdo2,
                       N.iter = 1000,
                       alpha = 0.05,
                       p.adj = "oleson",
                       cores = 4)
