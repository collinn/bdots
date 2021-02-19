


tt <- boot.l$curveList
X1 <- tt$CI$curveMat # 1000 x 501
X2 <- tt$NH$curveMat # 1000 x 501



fastT <- function(otu, indg1, indg2){
  l1 <- length(indg1)
  l2 <- length(indg2)
  S1 <- (rowSums((otu[,indg1, drop = FALSE] - rowMeans(otu[,indg1, drop = FALSE]))^2) / (l1 - 1))
  S2 <- (rowSums((otu[,indg2, drop = FALSE] - rowMeans(otu[,indg2, drop = FALSE]))^2) / (l2 - 1))
  sigmas <- sqrt ( S1/l1 + S2/l2 )
  obsT <- (rowMeans(otu[,indg1, drop = FALSE]) - rowMeans(otu[,indg2, drop = FALSE])) / sigmas
  return(obsT)
}

t1 <- vector("list", ncol(X1))
alpha <- 0.05
for (i in seq_len(ncol(X1))) {
  t1[[i]] <- t.test(X1[, i], X2[, i], conf.level = 1 - alpha)
}

x1 <- X1[, 1]
x2 <- X2[, 2]

t1 <- x1 - mean(x1)
t2 <- x2 - mean(x2)

# s1 <- sd(x1)
# ss1 <- sqrt(mean((x1 - mean(x1))^2))
s1 <- sum((x1 - mean(x1))^2)/(length(x1) - 1)
s2 <- sum((x2 - mean(x2))^2)/(length(x2) - 1)

(mean(x1) - mean(x2)) / sqrt((s1/length(x1) + s2/length(x2))) -> myt
t.test(x1, x2)

s1
ss1
sss1

library(microbenchmark)
microbenchmark(
  s1 <- sd(x1),
  ss1 <- sqrt(mean((x1 - mean(x1))^2)),
  sss1 <- sqrt(sum((x1 - mean(x1))^2)/(length(x1) - 1)),
  times = 5000
)
