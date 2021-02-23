data <- cohort_unrelated
dat <- data.table()
dat$subject <- data[[subject]]
dat$time <- data[[time]] %>% as.numeric()
dat$y <- data[[y]] %>% as.numeric()
#group <- c("Group", "LookType")

## Set group variables in data.tableb+
for(gg in seq_along(group)) {
    set(dat, j = group[gg], value = data[[group[gg]]])
}
#
dd <- dat[subject == 1 & LookType == "Cohort", ]
dd <- droplevels(dd)

## Getting parameters once
pars1 <- dgaussPars(dat[Group == 50, ], TRUE)
names(pars1) <- paste0(names(pars1), 1)
pars2 <- dgaussPars(dat[Group == 65, ], TRUE)
names(pars2) <- paste0(names(pars2), 2)

tt <- copy(dd)
tt$LookType <- NULL

tt2 <- dcast(tt, y + time ~ .)
           
           
ff <- quote(y ~ (Group == 50) * ((time < mu1) * (exp(-1 * (time - mu1) ^ 2 / (2 * sig11 ^ 2)) * (ht1 - base11) + base11) + (mu1 <= time) * (exp(-1 * (time - mu1) ^ 2 / (2 * sig21 ^ 2)) * (ht1 - base21) + base21)) + 
                (Group == 65) * ((time < mu2) * (exp(-1 * (time - mu2) ^ 2 / (2 * sig12 ^ 2)) * (ht2 - base12) + base12) + (mu2 <= time) * (exp(-1 * (time - mu2) ^ 2 / (2 * sig22 ^ 2)) * (ht2 - base22) + base22)))

test <- curveFitter(dd, ff, c(pars1, pars2), rho = 0.9, refits = 2)

test <- gnls(eval(ff), start = c(pars1, pars2), data = tt, correlation = corAR1(0.9), 
             control = gnlsControl(nlsTol = 0.01))
