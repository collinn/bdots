
x <- matrix(1:8, ncol = 2)

stripMat <- function(x) {
  #browser()
  if(nrow(x) > 1) {
    i <- 1 # indicates it happened
    rr <- vector("list", 1)
    for (j in seq_len(nrow(x))) {
      rr[j] <- stripMat(matrix(x[j, ], nrow = 1))
    }
    res <- c(list(rr), list(i))
  } else {
    if (!exists("i", parent.frame())) { # troublesome if length 1 and i exists. hmm
      i <- 0
      res <- c(list(x), list(i))
    } else {
      res <- list(x)
    }
  }
  res
}

stripMat(x)

i <- 4
y <- matrix(1:2, nrow = 1)
stripMat(y)
