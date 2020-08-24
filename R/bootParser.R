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

  ## Prep for subset
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



bootSubset <- function(l, bdObj) {
  subnames <- l[["subnames"]]
  subargs  <- l[["subargs"]]
  resNames <- l[['resNames']]

  # ouch (copy expensive)
  # remove columns we don't want
  # also, this shit is internal, order of columns doesn't matter
  # since it's not being returned
  bd <- subset(bdObj, select = c(resNames, subnames))
  # ss_vec <- vector("numeric", length = nrow(bd))
  for(i in seq_along(subnames)) {
    ss_vec <- bd[[subnames[i]]] %in% subargs[[i]]
    bd <- bd[ss_vec, ]
  }

  ## I'm still going to keep order
  bd[, c(resNames[1], subnames, resNames[-1]), with = FALSE]
}


###
## bdCall2Subset and recurToSubset should be
# combined. The only difference is that recurToSubset
# will return a length 1 list (as is the case for lhs[[3]])

## used for subsetting dt based on bdotsBoot formula
# takes g(n1, n2, ...) and returns c("g", "n1", "n2", ...)
bdCall2Subset <- function(x) {
  if(!is.call(x)) stop(paste0("invalid syntax:", x))
  x <- vapply(x, deparse1, character(1L))
  vv <- paste0("val", 1:(length(x) - 1))
  setNames(x, c("col", vv))
}

## Should document what this does
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
