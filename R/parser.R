
# ff <- diffs(y, condition(M,W)) ~ group(DLD, TD) + g2(n1, n2)
# ff <- y ~ group(DLD, TD)

## Get group and subset from something of form grp(n1, n2)
# needs to be length 1 or 2
breakupExpression <- function(ee) {
  if(!is.call(ee)) stop("Argument must be call")
  ee <- vapply(ee, as.character, character(1))
  if(length(ee) > 3) stop("Can only subset by one or two values in group: ", ee[[1]], .call = FALSE)
  ee
}

### Still needs some work
# behaves differently when there is a + on the RHS
# probably crazy ass errors that can occur too, idk


## Parse formula passed to bdotsBoot
## Can probably make something similar for bdotsFit?
# y ~ Subject * Time | groups ?
# need to figure out the subject/time arrangement
# fuuuh, or do it like survival packave
# (y, time) ~ subject | g1 + ...  <- yup

## Also possible that there are no RHS groups (shit)
## Actually, no, could be no groups to LHS
# or no, we wouldn't have diffs(y, g1(n1, n2)) ~ 1, it would be y ~ g1(n1, n2))



## important to keep in  mind, this function should return
# same object, regardless of input.
## Also. For final, export to CRAN bdots package,
# should rewrite everything from this there, so that I can
# keep all of the comments I've made detailing my thoughts on certain aspects of the code
# while retaining dignity of original authors, I do think it would be a helpful
# exercise to detail what was wrong, and how it change. How can I do that in
# a purely academic sense, with no discredit to Michael (Brad, eh. Nobody knows he did anything
# and what he added made it marginally better. It was an improvement to poorly written code,
# while retaining the logic and paradigm/style of the original. Ex. See, he saw this problem in v.1,
# he made v.2, see how v.2 > v.1, but then see why v3 is better suited to this problem)
bdotsParser <- function(ff) {
  print(class(ff))
  if(length(ff) != 3) stop("need y ~ x")

  ## Since formula always list, with first object `~`
  lhs <- ff[[2]]
  rhs <- ff[[3]]

  ## RHS First
  # Remove `+` from expression
  rhs <- Filter(function(x) !identical(as.character(x), "+"), rhs)

  rhs <- lapply(rhs, identity)


  fitGroups <- lapply(rhs, breakupExpression)
  grplen <- vapply(fitGroups, length, numeric(1))

  ## Let's be careful here
  if(sum(grplen == 3) > 1) stop("Only one grouping variable on LHS can have 2 values, i.e., y ~ g1(n1,n2) + g2(m1, m2) is not allowed")
  if(any(grplen > 3) | any(grplen < 2)) stop("Formula must have one or two values, i.e., y ~ g1(n1, n2, n3) is not allowed")

  ## LHS next
  ## Determine if diffs called
  diff <- grep("diffs", as.character(lhs[[1]]))

  if(diff) {

    ## Make sure diffs of correct formn
    if(length(lhs) != 3) {
      if(length(lhs) < 3) stop("Must supply grouping variable with diffs(), i.e., diffs(y, g1(n1, n2))")
      if(length(lhs) > 3) stop("Must only supply one grouping variable in diffs(), i.e., diffs(y, g1(n1, n2))")
    }
    y <- as.character(lhs[[2]])
    diffGroup <- breakupExpression(lhs[[3]])

    ## Must be 3, for group name and 2 values
    if(length(diffGroup) != 3) stop("Must supply 2 group values with diffs(), i.e., diffs(y, g1(n1, n2))")
    return(list(y = y, diffGroup = diffGroup, fitGroups = fitGroups, diffs = TRUE))
  } else {
    y <- as.character(lhs)
    return(list(y = y, diffGroup = NULL, fitGroups = fitGroups, diffs = FALSE))
  }

}

#bdotsParser(ff <- diffs(y, condition(M,W)) ~ group(DLD, TD) + g2(n1, n2))

## parses things like doubleGauss(concave = TRUE), logistic(), poly(degree = n), etc.

## This version of f/g works
f <- function(a) { # <- curve parse function
  #class(a) %>% print()
  a
}
g <- function(a) { # <- bdots, with a = doubleGauss(curve = TRUE)
  f(substitute(a))
}

ee <- g(doubleGauss(concave = TRUE, lickmyballs = TRUE, n = 8))


# bdotsFit(...) { curveParser(substitute(expr))}

#strsplit(ee, "[\\(\\)]")
## expr comes in as a call
expr <- ee

## This takes expression for curve type
# i.e., doubleGauss(concave = TRUE)
# poly(n = w/e)
# maybe someday get sophisticated with this
# return value is list named with curve, elements are arguments
## Keep
curveParser <- function(expr) {
  if(!is.call(expr)) stop("Invalid curve expression")
  expr <- unlist(strsplit(deparse(expr), "[\\(\\)]"))

  ## First of these is function name
  curve <- expr[1]
  if (!(curve %in% c("doubleGauss", "logistic", "poly"))) stop("Invalid curve type") # can probably remove this (probably should)

  arggs <- strsplit(expr[-1], ",") # get all arguments
  arggs <- lapply(arggs, function(x) gsub("^[ \t]+|[ \t]+$", "", x))

  # sanity check
  for(v in unlist(arggs)) {
    if(!grep("[=]", v)) stop(paste0(v, " is an invalid argument"))
    if(length(unlist(strsplit(v, "[=]"))) != 2) stop(paste0(v, " is an invalid assignment"))
  }

  setNames(arggs, curve)
}

### Lets test that output
# Shit yes, this evaluates everything and returns environment
## But instead assume it has input from curveParser
## Keep
# takes list, returns environment
makeCurveEnv <- function(val) {
  curveType <- names(val)
  val <- as.list(unlist(val, use.names = FALSE))
  myenv <- new.env()
  for(i in seq_along(val)) {
    #eval(parse(text = val[i]), envir = sys.frame(sys.nframe()))
    eval(parse(text = val[i]), envir = myenv)
  }
  myenv
}



## It's important that we call with substitute in function
ww <- function(a) {
  curveParser(substitute(a))
}

(ww(doubleGauss(concave = FALSE, eatmyshitter = TRUE)))
### Lets test that output
# Shit yes, this evaluates everything and returns environment
ww <- function(a) {
  val <- curveParser(substitute(a))
  curveType <- names(val)
  val <- as.list(unlist(val, use.names = FALSE))
  myenv <- new.env()
  for(i in seq_along(val)) {
    #eval(parse(text = val[i]), envir = sys.frame(sys.nframe()))
    eval(parse(text = val[i]), envir = myenv)
  }
  myenv
}

## This returns an environment
aa <- (ww(doubleGauss(concave = FALSE, eatmyshitter = TRUE)))

## Maybe I can combien this all in one place
curveParser2 <- function(expr) {
  if(!is.call(expr)) stop("Invalid curve expression")
  expr <- unlist(strsplit(deparse(expr), "[\\(\\)]"))

  ## First of these is function name
  curve <- expr[1]
  if (!(curve %in% c("doubleGauss", "logistic", "poly"))) stop("Invalid curve type") # can probably remove this (probably should)

  arggs <- strsplit(expr[-1], ",") # get all arguments
  arggs <- lapply(arggs, function(x) gsub("^[ \t]+|[ \t]+$", "", x))

  # sanity check
  for(v in unlist(arggs)) {
    if(!grep("[=]", v)) stop(paste0(v, " is an invalid argument"))
    if(length(unlist(strsplit(v, "[=]"))) != 2) stop(paste0(v, " is an invalid assignment"))
  }

  myenv <- new.env()
  arggs <- as.list(unlist(arggs, use.names = FALSE))
  myenv <- new.env()
  myenv$curveType <- curve
  for(i in seq_along(arggs)) {
    #eval(parse(text = val[i]), envir = sys.frame(sys.nframe()))
    eval(parse(text = arggs[i]), envir = myenv)
  }
  myenv
  #setNames(arggs, curve)
}

ww2 <- function(a) {
  curveParser2(substitute(a))
}

aa <- (ww2(doubleGauss(concave = FALSE, eatmyshitter = TRUE)))

## What we find is that frame_n gives frame of current env
# sys.frame(frame_n) allows us to list2env in function
# this shit works with mclapply
# also works for  parLapply
# f <- function(ll = 1) {
#   a <- 2
#   b <- 3
#   frame_n <- sys.nframe() # what frame am I in?
#   list2env(ll, envir = sys.frame(frame_n))
#   tt <- ls()
#   mm
# }
#
# h <- function() {
#   ll <- list(list(m = 1, mm = 2), list(m = 1, mm = 9))
#   #aa <- mclapply(ll, f)
#
#   cl <- makePSOCKcluster(2)
#   clusterExport(cl, c("ll", "f"))
#   aa <- parLapply(cl, ll, f)
#   stopCluster(cl)
#   aa
# }
#
# (h())

### Can we assign to env and export in parallel? YES
f <- function(x, ...) {
  #cat(doggy)
}
ww <- function(a) {
  val <- curveParser(substitute(a))
  curveType <- names(val)
  val <- as.list(unlist(val, use.names = FALSE))
  myenv <- new.env()
  for(i in seq_along(val)) {
    #eval(parse(text = val[i]), envir = sys.frame(sys.nframe()))
    eval(parse(text = val[i]), envir = myenv)
  }

  cl <- makePSOCKcluster(2) #, outfile = "")
  clusterExport(cl, c("f"))
  clusterExport(cl, ls(envir = myenv), envir = myenv) # this will export named  variables in myenv
  aa <- parLapply(cl, 1:10, f)
  stopCluster(cl)

}
(ww(doubleGauss(doggy = "doggy")))



