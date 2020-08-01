
## Let's test bdotsFit funcitoin

rm(list = ls())
library(nlme)
library(mvtnorm)
# Not for parallel, I could do clusterEvalQ(cl, source(".."))
source("~/packages/bdots/R/bdotsFit.R")
source("~/packages/bdots/R/bdotsFitter.R")
source("~/packages/bdots/R/parser.R")
source("~/packages/bdots/R/curveFitter.R")
source("~/packages/bdots/R/doubleGauss.R")
source("~/packages/bdots/R/logistic.R")
source("~/packages/bdots/R/helper.R")
source("~/packages/bdots/R/effectiveAlpha.R")
source("~/packages/bdots/R/findModifiedAlpha.R")
source("~/packages/bdots/R/ar1Solver.R")
library(nlme) # fuck me, not knowing I was missing this made debugging hard

load("~/bdots/demo.RData")
data <- cohort_unrelated
group = c("Group", "LookType")
time <- "Time"
subject <- "Subject"
y <- "Fixations"
concave <-  TRUE
conc <- TRUE
rho <-  0.9
refits <- 0
cores <-  1
verbose <-  FALSE
curveType <- c("doubleGauss")

head(data)

## was missing concave smh
system.time(res <- bdotsFit(cohort_unrelated,
                subject = "Subject",
                time = "Time",
                y = "Fixations",
                group = c("Group", "LookType"),
                curveType = doubleGauss(concave = TRUE),
                cor = TRUE,
                refits = 2,
                cores = 0,
                verbose = FALSE))
res2 <- res[Subject %in% c(1, 2, 3, 5, 7:12, 14:21, 23:26)]
## Keep only valid pairs for now
#res2 <- res[, 1:7]

# this removes possible pairs
#res2 <- res2[fitCode != 6, ]




## logistic
data(ci)
ci <- as.data.table(ci)
ci <- ci[LookType == "Target", ]
system.time(res.l <- bdotsFit(data = ci,
                              subject = "Subject",
                              time = "Time",
                              y = "Fixations",
                              group = "protocol",
                              curveType = logistic(),
                              cor = TRUE,
                              refits = 2))

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
#currdata$Curve <- ifelse(currdata$TrialType == "M", 1, 2)
currdata2 <- as.data.table(currdata)
currdata2 <- currdata2[Subject != 405, ]
currdata2 <- currdata2[Subject != 1699, ]
currdata2 <- currdata2[Subject != 1526, ]
N.iter <- 1000

## About 32 seconds
system.time(res.b <- bdotsFit(data = currdata2,
                              subject = "Subject",
                              time = "Time",
                              y = "Looks",
                              group = c("Group", "TrialType"),
                              curve = logistic(),
                              cor = TRUE,
                              refits = 2))

head(sort(res.b$R2))
#tt <- currdata[Subject == 405, ]
#tt1 <- tt[TrialType == "M", ]
#tt2 <- tt[TrialType == "W", ]



