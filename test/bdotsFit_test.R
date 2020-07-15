
## Let's test bdotsFit funcitoin

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
                curveType = "doubleGauss",
                cor = TRUE,
                refits = 2,
                cores = 0,
                concave = TRUE,
                verbose = FALSE))
