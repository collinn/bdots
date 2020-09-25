

## res.b is bdotsFit object
bdo <- copy(res.b)
bdCall <- attr(bdo, "call")
nn <- nn <- c(eval(bdCall[['subject']]), eval(bdCall[['group']]))



# split2 <- function(bdo, by, ...) {
#   oldAttr <- attributes(bdo)
#   res <- lapply(data.table:::split.data.table(bdo, by = by, ...), function(x) {
#     attributes(x) <- oldAttr
#     x
#   })
#   structure(.Data = res, class = c("bdObjList"))
# }
#
# rbindlist2 <- function(bdo, by, ...) {
#   oldAttr <- attributes(bdo[[1]])
#   class(bdo) <- "list"
#   bdo <- rbindlist(bdo)
#   attributes(bdo) <- oldAttr
#   bdo
# }
#
#
#
#
# setkeyv(bdo, nn)
mem_used()
## Let's verify this works
pryr::mem_change(r1 <- split(bdo, by = nn))
r1 <- split(bdo, by = nn)
class(r1)
rr1 <- r1[[1]]
x <- attr(rr1, "X")$X
head(x)
attr(rr1, "X")$X$dogmeat <- "dogmeat"
r2 <- rbindlist(r1)
attr(r2, "X")$X
### End verification

bdo2 <- copy(bdo)
class(bdo2) <- c("data.table", "data.frame")

library(microbenchmark)

## Apparently, split.data.table does preserve attributes, albiet reorderd
microbenchmark(
  r1 <- split2(bdo, by = nn),
  r2 <- split(bdo2, by = nn),
  times = 5
)

class(r1[[1]])
class(r2[[1]])

rr1 <- r1[[2]]
rr2 <- r2[[2]]

x1 <- attr(rr1, "X")$X
x2 <- attr(rr2, "X")$X

head(x1)
head(x2)

attr(rr1, "X")$X$dogmeat <- "dogmeat"
attr(rr2, "X")$X$dogmeat <- "dogmeat"

## Update was successful
mm1 <- r1[[2]]
mm2 <- r2[[2]]

x1 <- attr(mm1, "X")$X
x2 <- attr(mm2, "X")$X

head(x1)
head(x2)

test1 <- rbindlist2(r1)
test2 <- rbindlist(r2)

head(test1)
head(test2)

## And here, xx2 is null.
xx1 <- attr(test1, "X")$X
xx2 <- attr(test2, "X")$X

head(xx1)
head(xx2)


## Investigating size of objects
tt <- res.b[1, ]$fit[[1]]
class(tt)

rr <- tt$modelStruct
object_size(rr)

rr2 <- rr$corStruct
class(rr2)
ff <- attr(rr2, "factor")

class(tt)

object_size(tt)
attr(tt$modelStruct$corStruct, "factor") <- NULL
object_size(tt)

cc <- coef(tt)
fitt <- fitted.values(tt)
vb <- tt$varBeta


