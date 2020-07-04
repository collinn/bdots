
ff <- diffs(y, condition(M,W)) ~ group(DLD, TD) + g2(n1, n2)
ff <- y ~ group(DLD, TD)

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

## Also possible that there are no RHS groups (shit)
# or no, we wouldn't have diffs(y, g1(n1, n2)) ~ 1, it would be y ~ g1(n1, n2))
bdotsParser <- function(ff) {
  print(class(ff))
  if(length(ff) != 3) stop("need y ~ x")
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

bdotsParser(ff <- diffs(y, condition(M,W)) ~ group(DLD, TD) + g2(n1, n2))



