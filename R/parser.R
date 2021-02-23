
## This is for parsing formula for bdotsBoot
## Syntax as follows:
## LHS can be only response, or it can be diffs, with a group of two elements
## RHS  MUST contain at least one group with 2 elements and arbitrary number of single element groups
# y ~ g(n1, n2)
# diffs(y, g1(n1, n2)) ~ g2(m1, m2)
# Note: from above 2 lines, length(lhs) can only be 1 or 3
# ... ~ g1(n1, n2) + g2(m1) + g3(m2) + ...

# ff1 <- diffs(y, condition(M,W)) ~ group(DLD, TD) + g2(n1) + g3(m1, m2, m3, m4) + g4(r1)
# ff1 <- diffs(y, TrialType(M,W)) ~ Group(LI, TD)
# ff1 <- y ~ Group(LI, TD) + TrialType(M)
#
# ff <- ff1
bootParser <- function(ff, bdObj) {
  if (!inherits(ff, "formula")) stop("Must supply a formula to bdotsBoot")

  ## is this dangerous
  lhs <- ff[[2]]
  rhs <- ff[[3]]

  ## Process LHS
  if(length(lhs) == 1) {
    resp <- deparse1(lhs)
    diffs <- NULL
    inner <- NULL
  } else if (length(lhs) == 3) {
    diffs <- deparse1(lhs[[1]])
    if (!identical(diffs, "diffs")) stop("invalid formula at diffs")
    resp <- deparse(lhs[[2]]) # don't actually need this
    inner <- bdCall2Subset(lhs[[3]])
    if (length(inner) != 3) stop("Must subset by exactly two group values for diff")
  } else {
    stop("invalid formula sytax on LHS")
  }

  ## RHS nicely wrapped in functions
  ss <- recurToSubset(rhs)

  ## ensure valid syntax for rhs
  vv <- vapply(ss, length, numeric(1))
  if ((sum(vv == 3) != 1) | any(vv > 3)) stop("Exactly one group on RHS of formula must have two values")

  ## Get outer group
  outerDiff <- ss[vv == 3][[1]]["col"]
  names(outerDiff) <- NULL

  ## Prep for subset
  # clear case where having inner <-bdCall2Subset OK to be length 1 list
  # returns a list of length (numargs). Each list element is a list containing
  # the column subset name and the values (perhaps better as a vector?)
  # As in, what's wrong with  just c(list(inner), ss) ? Should verify later
  ww <- lapply(c(list(inner), ss), function(x) {
    c(x[1], list(x[-1]))
  })
  if (is.null(inner)) ww <- ww[-1]

  ## Reserved names from bdObj
  # this feels precarious as best
  resNames <- c(attr(bdObj, "names")[1],
                colnames(bdObj)[(ncol(bdObj) - 3):ncol(bdObj)])

  ## get names
  subnames <- vapply(ww, `[[`, character(1), 1)
  subargs <- lapply(ww, `[[`, 2)
  nn <- intersect(colnames(bdObj), subnames)
  if (!identical(intersect(subnames, nn), subnames)) stop("Provided group names not in bdObj")

  list(subnames = subnames,
       subargs = subargs,
       resNames = resNames,
       outerDiff = outerDiff,
       innerDiff = inner[['col']])
}




# ff <- diffs(y, condition(M,W)) ~ group(DLD, TD) + g2(n1, n2)
# ff <- y ~ group(DLD, TD)

## Get group and subset from something of form grp(n1, n2)
# needs to be length 1 or 2
# I think this function is obsolete. Not sure yet
breakupExpression <- function(ee) {
  if(!is.call(ee)) stop("Argument must be call")
  ee <- vapply(ee, as.character, character(1))
  if(length(ee) > 3) stop("Can only subset by one or two values in group: ", ee[[1]], .call = FALSE)
  ee
}


## important to keep in  mind, this function should return
# same object, regardless of input.
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



## Maybe I can combien this all in one place
curveParser2 <- function(expr) {
  if(!is.call(expr)) stop("Invalid curve expression")
  expr <- unlist(strsplit(deparse(expr), "[\\(\\)]"))

  ## First of these is function name
  curve <- expr[1]

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
    eval(parse(text = arggs[i]), envir = myenv)
  }
  myenv
  #setNames(arggs, curve)
}


