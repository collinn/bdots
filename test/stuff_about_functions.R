# Meant as a playground/resource for parser, not indended to be exported or  used later

# as.list(formula) - length 2 or 3, depending on existence of LHS (which we will always have)
# each element of this list will be a call, except `~` which is a name
# in ff1 group(DLD, DT) is a call. In subject | group + group2,
# `subject` and `|` are names, `group + group2` is a call
# as.character(call) will turn it into whatever length (and I'm not quite sure
# yet what determines that number), HOWEVER, deparse will turn it directly into a
# character vector of length 1, which might be ideal.

ff1 <- diffs(y, condition(M,W)) ~ group(DLD, TD) + g2(n1) + g3(m1, m2, m3, m4) + g4(r1)
ff2 <- y ~ group(DLD, TD)
ff <- ff1
ff <- ff2


(lhs <- ff[[2]])
(rhs <- ff[[3]])

as.list(lhs)
lapply(lhs, class)

## Rhs can be of arbitrary length, seperated by `+`
as.list(rhs)
lapply(rhs, class)

## This is ideal for RHS
# here's whats to note. If we have  arbitrary g1 + g2 + g3 + g4
# the elements are [1] `+`, [2] g1 + g2 + g3, [3] g4
# so it goes from rhs to lhs of the first function on the right
x <- rhs
recurToSubset <- function(x) {
  allnames <- identical(rep("name", length(x)),
                        vapply(x, class, character(1)))
  if (length(x) == 1 | allnames) {
    outer <- bdCall2Subset(x)
  } else {
    if (!identical(as.symbol("+"), x[[1]])) stop("invalid formula syntax on rhs")
    outer <- lapply(x[-1], recurToSubset)
  }
  unzipList(outer)
}


(tt <- recurToSubset(rhs))




## If length(lhs) == 1, this is just the response variable
## If length(lhs) == 3, it has to be diffs (should still check)
if(length(lhs) == 1) {
  resp <- deparse1(lhs)
} else if (length(lhs) == 3) {
  diffs <- deparse1(lhs[[1]])
  if(!identical(diffs, "diffs")) stop("invalid formula at diffs")
  resp <- deparse(lhs[[2]]) # don't actually need this

  ## This part below needs to be generalized (item a1)
  inner <- bdCall2Subset(lhs[[3]])
} else {
  stop("invalid formula sytax")
}

# item a1
(lhs <- ff[[2]])
outer <- lhs[[3]]

## used for subsetting dt based on bdotsBoot formula
# takes g(n1, n2, ...) and returns c("g", "n1", "n2", ...)
bdCall2Subset <- function(x) {
  if(!is.call(x)) stop(paste0("invalid syntax:", x))
  x <- vapply(x, deparse1, character(1L))
  vv <- paste0("val", 1:(length(x) - 1))
  setNames(x, c("col",))
}


## Stolen from purrrrrrrr
vec_depth <- function(x) {
  if (is.null(x)) {
    0L
  } else if (is.atomic(x)) {
    1L
  } else if (is.list(x)) {
    depths <- vapply(x, vec_depth, numeric(1))
    1L + max(depths, 0L)
  } else {
    stop("'x' must be a vector")
  }
}





### BELOW THIS IS OBSOLETE ###

# function to turn f(a, b, ...) into  c(f, a, b, ...)
# but this could also be diffs(y, g1(n1, n2)),
# length 2 if it's normal call, i.e., [1] f  [2] a, b ...)
# and length 3 if its diffs, [1] diffs [2] y, group [3] n1, n2) (can't be diffs and length 2)
# ugh. What if group only had two names. Then we wouldn't want to specify g1(n1, n2), we would just want g1
# anything else would be invalid syntax
cc <- lhs
turnCallintoPars <- function(cc) {
  ccc <- deparse(cc)
  ccsp <- strsplit(ccc, "\\(")
}
