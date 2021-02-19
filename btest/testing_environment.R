## All of the variables that are typically used  during testing that can be loaded

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

dat <- data.table()
dat$subject <- data[[subject]] %>% as.character()
dat$time <- data[[time]] %>% as.numeric()
dat$y <- data[[y]] %>% as.numeric()

## Set group variables in data.table
for(gg in seq_along(group)) {
  set(dat, j = group[gg], value = data[[group[gg]]])
}

fakefunctionwithlongname <- function(dat, ct, cc, rho, eee, eee3) {
  q1 <- list(res = list(a = dat, b = class(dat), c = cc, curveType = rho))
  return(q1)
}

## Set key
setkeyv(dat, c("subject", group))

## HOLY SHIT THIS WORKS
## DOUBLE HOLY SHIT, USING SD PASSES ENTIRE DATA TABLE SUBSET !!
## just replace f with bdotsFitter, once finished, goddamn, doing my job for me
## .SD going in only contains exactly what is needed to fit function
rr <- dat[, list(fakefunctionwithlongname(.SD, curveType, concave, cor, rh0, verbose)), by = group,
          .SDcols = c("subject", "time", "y")]

## test dataset - only contains subject, time, and y
tdat <- rr$V1[1][[1]][[1]]

ff <- list.files("~/packages/bdots/R", full.names = TRUE)
for(f in ff) source(f)
rm(f)
rm(ff)

save.image("~/packages/bdots/data/test_env.RData")






