

g <- function(...) {
  c(as.list(environment()), list(...))
}

f <- function(a = 1, b, ...) {
  c <- 3
  x <- 1:10
  res <- g(...)
}

r <- f(a = 2, b = 3, cat = 4)

h <- function(a = 1, b = 2) {
  a + b
}

k <- function(...) {
  arggs <- c(as.list(environment()), list(...))
  do.call(h, arggs)
}
k(a = 4, b = 9, cat = "dog")


## This returns class 'call'
f <- function(a) {
  substitute(a)
}
rr <- f(doubleGauss(concave = FALSE, balls = TRUE))
names(rr)


f <- function(a, b, ...) {
  match.call()
}

f(a = 2, b = 4, cat = "dog")



## Save function as attrbute

g <- function(dog = TRUE) {
  if (dog) print("dog") else print("not dog")
}

f <- function(a = 1, func) {
  dt <- as.data.table(matrix(rnorm(9), nrow = 3))
  b

  structure(dt, class = c("myclass", "data.table", "data.frame"),
            retFun = g)
}

b <- f(a = 2, g)

library(microbenchmark)

microbenchmark(
  x <- dat$time,
  x <- dat[["time"]]
)














