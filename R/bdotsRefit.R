
# bdObj - object from bdotsFit
# fitCode - min code for refit, i.e.,  fitCode = 3 will do 3:6
bdObj <- res.b
fitcode <- 1

## Probably want to create diagnostic information
# that can be output or saved along with this
# also probably worth recording which were updated?
bdotsRefit <- function(bdObj, fitcode = 0, ...) {
  if (is.null(fitcode)) {
    fitcode <- readline(prompt = "fitcode: ")
  }
  X_orig <- attr(bdObj, "X")
  attr(bdObj, "X") <- NULL
  if (is.null(X)) {
    stop("Dataset must be provided")
  }

  ## Because of factors, which I should really think about removing
  idx <- which(as.numeric(bdObj$fitCode) >= fitcode + 1)
  #bdObj <- bdObj[as.numeric(fitCode) >= fitcode + 1, ]

  ## Subset X by selected bdObj
  bdCall <- attr(bdObj, "call")
  nn <- c(eval(bdCall[['subject']]), eval(bdCall[['group']]))
  bdNames <- do.call(paste, c(bdObj[idx, nn, with = FALSE], sep = "-"))
  XNames <- do.call(paste, c(X_orig[, nn, with = FALSE], sep = "-"))
  x_idx <- which(XNames %in% bdNames)
  X <- X_orig[x_idx, ]

  X <- split(X, by = nn)
  bdObj2 <- split(bdObj[idx, ], by = nn)

  ## They should be identical, but who knows
  if (!identical(names(X), names(bdObj2))) {
    if (intersect(names(X), names(bdObj2)) != names(X)) stop("error 45367")
    bdObj2 <- bdObj2[sort(names(bdObj2))]
    X <- X[sort(names(bdObj2))]
  }

  new_bd <- lapply(names(bdObj2), bdUpdate)
  new_bd <- rbindlist(new_bd)

  for (i in seq_len(nrow(test))) {
    bdObj[idx[i], ] <- new_bd[i, ]
  }
  # ## possibly do by names so as not to have to map?
  # test <- Map(function(bdo, x) {
  #   # print current status of observation
  #   # give plot
  #   # prompt to give set of parameters
  #   # display new fit/plot
  #   # keep new pars or reprompt
  #
  # }, bdObj, X)

  # then, the results of test will be
  # used to update the original bdObj

}


# bdo <- bdObj[[1]]
# x <- X[[1]]
# nn <- names(bdObj)[1]

## A lot of bits of these can be abstracted to functions
# that way we can just jitter as often as we want (or w/e else)
# note, this will also be found in bdotsFit too
# Perhaps throw in  option to "Trash" an observation
bdUpdate <- function(nn) {
  ## get whatever
  bdo <- bdObj2[[nn]]
  x <- X[[nn]]
  attr(bdo, "X") <- x

  ## Getting curve function
  bdCall <- attr(bdo, "call")
  time <- attr(bdo, "time")
  crvFun <- curve2Fun(bdCall[['curveType']])

  ## This is a terrible work around
  #  to get correct names in x for curveFitter
  set(x, j = c("y", "time"),
      value = x[,c(bdCall[['y']], bdCall[['time']]), with = FALSE])

  rho <- attr(bdo, "rho")

  accept <- FALSE
  while (!accept) {
    plot(bdo, gridSize = 1)
    r2 <- round(bdo[['R2']], 3)
    ar1 <- bdo[['AR1']]
    fc <- bdo[['fitCode']]
    fit <- bdo[['fit']][[1]]
    msg <- paste0("Subject: ", nn, "\nR2: ", r2, "\nAR1: ",
                  as.logical(ar1), "\nrho: ", rho,
                  "\nfitCode: ", fc, "\n\n")
    msg <- c(msg, "Model Parameters:\n")
    cat(msg)
    print(oldPars <- coef(fit))

    ## Maybe add in future ability to change row
    rf_msg <- paste0("\nActions:\n",
                     "1) Keep as is\n",
                     "2) Jitter parameters\n",
                     "3) Adjust parameters manually\n")
    cat(rf_msg)
    resp <- readline("Choose (1-3): ")
    if (resp == 1) {
      accept <- TRUE
      break
    } else if (resp == 2) {
      newPars <- jitter(coef(fit))
    } else if (resp == 3) {
      newPars <- oldPars
      for (pname in names(oldPars)) {
        cat("Current value:\n")
        print(oldPars[pname])
        newPars[pname] <- readline(paste0("New value for ", pname, ": "))
      }
      class(newPars) <- "numeric"
    }

    ## repeat of inside of bdotsFit (make this a function? - yes)
    # then need to update to bdFit as well
    result <- bdotsFitter(dat = x, curveType = crvFun, rho = rho, params = newPars)
    # Or just copy to preserve relevant bits

    new_bdo <- copy(bdo)
    new_bdo$fit <- I(list(result['fit']))
    new_bdo$R2 <- result[['R2']]
    new_bdo$AR1 <- (result[['fitCode']] < 3)
    new_bdo$fitCode <- result[['fitCode']]

    both_bdo <- rbindlist(list(bdo, new_bdo))
    attributes(both_bdo) <- attributes(bdo)
    plot(both_bdo, gridSize = "refit")

    keep <- readline("Keep new fit? (y/n): ")
    if (length(grep("y", keep, ignore.case = TRUE)) == 1) {
      bdo <- new_bdo
      accept <- TRUE
    }
  }
  bdo
}











