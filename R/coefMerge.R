
#
# dat <- fread(file = "btest/AUDIO_7_Cohort_out.txt")
#
# bdObj <- copy(res)
# infile <- coefWriteout(bdObj)
#
# infile <- infile[1:20, ]

mergecoef <- function(bdObj, infile) {

  infile <- as.data.table(infile)

  idcols <- getIdentifierCols(bdObj)

  setkeyv(bdObj, idcols)
  setkeyv(infile, idcols)

  bdo_obs <- do.call(paste, bdObj[, idcols, with = FALSE])
  infile_obs <- do.call(paste, infile[, idcols, with = FALSE])

  ## Perfect new idx (double assignment SUCH a bad idea)
  idx <- which(bdo_obs %in% infile_obs)

  bdObj2 <- split(bdObj[idx, ], by = idcols)
  newpars <- split(infile, by = idcols)

  if (!identical(order(names(bdObj2)), order(names(newpars)))) {
    stop("bad news bears in mergecoef")
  } else {
    newpars <- newpars[names(bdObj2)]
  }

  newpars <- lapply(newpars, function(x) {
    x <- subset(x, select = setdiff(colnames(x), idcols))
    x <- unlist(x)
    if (any(is.na(x))) x <- NULL
    x
  })

  cores <- attr(bdObj, "call")$cores
  if (Sys.info()['sysname'] == "Darwin") {
    cl <- makePSOCKcluster(cores, setup_strategy = "sequential")
  } else {
    cl <- makePSOCKcluster(cores)
  }
  invisible(clusterEvalQ(cl, {library(bdots)}))

  new_bd <- clusterMap(cl, bdRefitter, bdo = bdObj2, params = newpars, getCovOnly = TRUE)
  stopCluster(cl)

  return(list(new_bd = new_bd, idx = idx))
}





