# Meant as a playground/resource for parser, not indended to be exported or  used later

# as.list(formula) - length 2 or 3, depending on existence of LHS (which we will always have)
# each element of this list will be a call, except `~` which is a name
# in ff1 group(DLD, DT) is a call. In subject | group + group2, 
# `subject` and `|` are names, `group + group2` is a call
# as.character(call) will turn it into whatever length (and I'm not quite sure
# yet what determines that number), HOWEVER, deparse will turn it directly into a 
# character vector of length 1, which might be ideal.

ff1 <- diffs(y, condition(M,W)) ~ group(DLD, TD) + g2(n1, n2)
ff2 <- y ~ group(DLD, TD)
ff3 <- c(y, time) ~ subject | group + group2


ff <- ff3

lhs <- ff[[2]]
rhs <- ff[[3]]


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