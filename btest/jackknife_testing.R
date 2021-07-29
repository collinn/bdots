
setwd("~/projects/noisy_sentence/")

jk <- fread("~/projects/noisy_sentence/data/c-uJK.txt")
njk <- fread("~/projects/noisy_sentence/data/c-u.txt")

jk <- jk[rTime >= 0 & rTime <= 1500, ]
njk <- njk[rTime >= 0 & rTime <= 1500, ]


## Only expect noise
jk <- jk[Trialcodecondition == "ExpectNoise", ]
njk <- njk[Trialcodecondition == "ExpectNoise", ]

jk <- jk[, CUR := mean(CUR), by = c("subjectID", "Trialcodecondition", "rTime")]
njk <- njk[, CUR := mean(CUR), by = c("subjectID", "Trialcodecondition", "rTime")]
njk <- unique(njk)

sub3_jk <- jk[subjectID == 3, ]
sub3_njk <- njk[subjectID != 3, ]
sub3_njk <- sub3_njk[subjectID != 3, ]
sub3_njk$subjectID <- 3
sub3_njk <- unique(sub3_njk)

## Confirms that jk happens as expected
plot(sub3_jk$CUR, sub3_njk$CUR)

tt <- jk[subjectID %in% c(3, 4), ]

tt3 <- tt[subjectID == 3, ]
tt3s <- copy(tt)
tt3s$subjectID <- 3

res <- bdotsFit(data = sub3_jk,
                subject = "subjectID",
                time = "rTime",
                y = "CUR",
                group = c("Trialcodecondition"),
                curveType = doubleGauss(concave = TRUE),
                cor = TRUE,
                numRefits = 1,
                cores = 8,
                verbose = FALSE)

res2 <- bdotsFit(data = sub3_njk,
                subject = "subjectID",
                time = "rTime",
                y = "CUR",
                group = c("Trialcodecondition"),
                curveType = doubleGauss(concave = TRUE),
                cor = TRUE,
                numRefits = 1,
                cores = 8,
                verbose = FALSE)

fit1 <- res$fit[[1]]
fit2 <- res2$fit[[1]]


## Try with 3 subjects and see how vars look
## nope here need to jackknife with only whats included - 4,6,7 are built with more than each other
ttjk <- jk[subjectID %in% c(4, 6, 7), ]
ttnjk <- njk[subjectID %in% c(4, 6, 7), ]

res <- bdotsFit(data = ttjk,
                subject = "subjectID",
                time = "rTime",
                y = "CUR",
                group = c("Trialcodecondition"),
                curveType = doubleGauss(concave = TRUE),
                cor = TRUE,
                numRefits = 1,
                cores = 8,
                verbose = FALSE)

res2 <- bdotsFit(data = ttnjk,
                 subject = "subjectID",
                 time = "rTime",
                 y = "CUR",
                 group = c("Trialcodecondition"),
                 curveType = doubleGauss(concave = TRUE),
                 cor = TRUE,
                 numRefits = 2,
                 cores = 8,
                 verbose = FALSE)


## ---
theta <- colMeans(coef(res))
theta_m <- matrix(theta, nrow = 3, ncol = 6, byrow = TRUE)

n <- 3
P <- n*theta_m - (n-1)*coef(res)
Pbar <- colMeans(P)

V <- t(P-Pbar) %*% (P - Pbar)
Vp <- (1 / (n * (n-1))) * V
Vinv <- MASS::ginv(Vp)

fit1 <- res[1, ]$fit[[1]]$varBeta
fit2 <- res2[1, ]$fit[[1]]$varBeta

cfit1 <- chol(fit1)
t(cfit1) %*% (cfit1)

t(cfit1) %*% Vp %*% cfit1
fit1
fit2



makeCurveFun <- function(bdObj) {
  time <- attr(bdObj, "call")[['time']]
  f_bod <- attr(bdObj, "formula")[[3]]
  f_args <- c(colnames(coef(bdObj)), time)
  f_args <- setNames(as.pairlist(rep("", length(f_args))), f_args)
  eval(call("function", f_args, f_bod), parent.frame())
}

ff <- makeCurveFun(res)
gr <- numDeriv::grad(ff, coef(res)[1, ])



### Let's just brute try it

jk <- fread("~/projects/noisy_sentence/data/c-uJK.txt")
njk <- fread("~/projects/noisy_sentence/data/c-u.txt")

jk <- jk[rTime >= 0 & rTime <= 1500, ]
njk <- njk[rTime >= 0 & rTime <= 1500, ]


## Only expect noise
jk <- jk[Trialcodecondition == "ExpectNoise", ]
njk <- njk[Trialcodecondition == "ExpectNoise", ]

jk <- jk[, CUR := mean(CUR), by = c("subjectID", "Trialcodecondition", "rTime")]
njk <- njk[, CUR := mean(CUR), by = c("subjectID", "Trialcodecondition", "rTime")]
njk <- unique(njk)

njk_en <- njk[Trialcodecondition == "ExpectNoise", ]
nsub <- unique(njk_en$subjectID)
df_list <- vector("list", length = length(nsub))
for (i in seq_along(nsub)) {
  df_list[[i]] <- njk_en[subjectID != nsub[i], ]
  df_list[[i]]$subjectID <- nsub[i]
}

new_df <- data.table::rbindlist(df_list)

## "unjackknifed"
res <- bdotsFit(data = new_df,
                subject = "subjectID",
                time = "rTime",
                y = "CUR",
                group = c("Trialcodecondition"),
                curveType = doubleGauss(concave = TRUE),
                cor = TRUE,
                numRefits = 2,
                cores = 8,
                verbose = FALSE)

#  jackknifed
res2 <- bdotsFit(data = jk,
                 subject = "subjectID",
                 time = "rTime",
                 y = "CUR",
                 group = c("Trialcodecondition"),
                 curveType = doubleGauss(concave = TRUE),
                 cor = TRUE,
                 numRefits = 2,
                 cores = 8,
                 verbose = FALSE)

res3 <- bdotsFit(data = njk,
                 subject = "subjectID",
                 time = "rTime",
                 y = "CUR",
                 group = c("Trialcodecondition"),
                 curveType = doubleGauss(concave = TRUE),
                 cor = TRUE,
                 numRefits = 2,
                 cores = 8,
                 verbose = FALSE)

fit1 <- res[2, ]$fit[[1]]
fit2 <- res2[2, ]$fit[[1]]
fit3 <- res3[2, ]$fit[[1]]

y3 <- jk[subjectID == 3, CUR]
yy3 <- fitted.values(res[1, ]$fit[[1]])
yy3 <- matrix(yy3, ncol = 376, byrow = TRUE)
yy3 <- colMeans(yy3)
rr <- matrix(resid(res[1, ]$fit[[1]]), ncol = 376, byrow = TRUE)
rr <- sum(colMeans(rr)^2)
ssy <- sum((y3-mean(y3))^2)
1 - rr/ssy

y3 <- new_df[subjectID == 3, CUR]
rr <- resid(res[1, ]$fit[[1]])
y3m <- matrix(y3, ncol = length(unique(jk$rTime)))
ssy <- sum((y3 - mean(y3))^2)
1 - rr/ssy


## "unjackknifed"
res4 <- bdotsFit(data = new_df[subjectID == 3, ],
                subject = "subjectID",
                time = "rTime",
                y = "CUR",
                group = c("Trialcodecondition"),
                curveType = doubleGauss(concave = TRUE),
                cor = TRUE,
                numRefits = 2,
                cores = 8,
                verbose = FALSE)


test <- copy(res)

test[subjectID %in% 10:20, Trialcodecondition := "Fake"]

boottest <- bdotsBoot(CUR ~ Trialcodecondition(ExpectNoise, Fake),
                      bdObj = test,
                      Niter = 1000,
                      alpha = 0.05,
                      cores = 8)



jk <- fread("~/projects/noisy_sentence/data/c-uJK.txt")
njk <- fread("~/projects/noisy_sentence/data/c-u.txt")

jk <- jk[rTime >= 0 & rTime <= 1500, ]
njk <- njk[rTime >= 0 & rTime <= 1500, ]

njk3 <- njk[subjectID %in% 4:10, ]
jk3 <- jk[subjectID %in% 4:10, ]

res_njk3 <- bdotsFit(data = njk3,
                 subject = "subjectID",
                 time = "rTime",
                 y = "CUR",
                 group = c("Trialcodecondition"),
                 curveType = doubleGauss(concave = TRUE),
                 cor = TRUE,
                 numRefits = 2,
                 cores = 8,
                 verbose = FALSE,
                 jackknife = TRUE)

res_jk3 <- bdotsFit(data = jk3,
                    subject = "subjectID",
                    time = "rTime",
                    y = "CUR",
                    group = c("Trialcodecondition"),
                    curveType = doubleGauss(concave = TRUE),
                    cor = TRUE,
                    numRefits = 2,
                    cores = 8,
                    verbose = FALSE,
                    jackknife = FALSE)

## does this work?
# no, because not using the original data matrix (damn)
refit <- bdotsRefit(res_njk3, quickRefit = TRUE)
