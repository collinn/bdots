
## Let's test bdotsFit funcitoin

rm(list = ls())
source("~/packages/bdots/R/bdotsFit.R")
source("~/packages/bdots/R/bdotsFitter.R")
source("~/packages/bdots/R/parser.R")
source("~/packages/bdots/R/curveFitter.R")
source("~/packages/bdots/R/doubleGauss.R")
source("~/packages/bdots/R/logistic.R")

load("~/bdots/demo.RData")
data <- cohort_unrelated
group = c("Group", "LookType")
time <- "Time"
subject <- "Subject"
y <- "Fixations"
concave <-  TRUE
conc <- TRUE
rho <-  0.9
cor <-  TRUE
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






