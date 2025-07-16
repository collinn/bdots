

library(bdots)

res <- bfit(data = cohort_unrelated,
            subject = "Subject",
            time = "Time",
            y = "Fixations",
            group = c("Group", "LookType"),
            curveFun = doubleGauss2(concave = TRUE),
            numRefits = 2,
            cores = 8)

res <- brefit(res, quickRefit = TRUE)
res <- brefit(res)

plot(res[1:4, ])

#############################
## Single difference first ##
#############################

boot1 <- bboot(formula = Fixations ~ Group(50, 65) + LookType(Cohort),
               bdObj = res,
               Niter = 500,
               alpha = 0.05,
               permutation = FALSE,
               padj = "oleson",
               cores = 8)

boot2 <- bboot(formula = Fixations ~ Group(50, 65) + LookType(Cohort),
               bdObj = res,
               Niter = 500,
               alpha = 0.05,
               padj = "oleson",
               cores = 8)

plot(boot1)
plot(boot2)


## Diff of diff (eww)

bdiff <- bboot(formula = diffs(Fixations, LookType(Cohort, Unrelated_Cohort)) ~ Group(50, 65),
               bdObj = res,
               Niter = 200,
               alpha = 0.05,
               permutation = FALSE,
               padj = "oleson",
               cores = 8)

bdiff2 <- bboot(formula = diffs(Fixations, LookType(Cohort, Unrelated_Cohort)) ~ Group(50, 65),
               bdObj = res,
               Niter = 200,
               alpha = 0.05,
               cores = 8)
