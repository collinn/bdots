#' Refit Observations Returned from bdotsFit
#'
#'
#'
#' @param bdObj An object of class 'bdotsObj' returned from \code{bdotsFit}
#' @param fitCode A length one numeric indicating observations to refit. See Details
#' @param ... not used
#'
#' @return Returns bdObj with updated fit
#'
#' @details fitCode indicates lower bound on observations to refit. For example,
#' if \code{fitCode = 4}, \code{bdotsRefit} will prompt user to refit all
#' observations with fitCode = 4, 5, 6.
#' @import data.table
#' @export
bdotsRefit <- function(bdObj, fitCode = 1, ...) {
  if (is.null(fitCode)) {
    fitCode <- readline(prompt = "fitCode: ")
  }
  # X <- attr(bdObj, "X")
  # attr(bdObj, "X") <- NULL
  if (is.null(attr(bdObj, "X")$X)) {
    stop("Dataset must be provided")
  }

  ## These uniquely identify each fit
  bdCall <- attr(bdObj, "call")
  nn <- c(eval(bdCall[['subject']]), eval(bdCall[['group']]))

  ## Because of factors, which I should really think about removing
  fitcode <- fitCode
  idx <- which(bdObj$fitCode >= fitcode)
  if (length(idx) == 0) {
    message("All observations have fitCode = 0. Nothing to refit :)")
    return(bdObj)
  }
  bdObj2 <- split(bdObj[idx, ], by = nn)
  new_bd <- lapply(bdObj2, bdUpdate)

  ## Remove deleted observations
  # null_idx - which fits were removed
  # new_bd - compacted new_bd list
  # rmv_names - names from bd list being removed
  # rmv_sub_id - Subject identifier of removed
  # -- this is to provide option to remove all of subject
  # idx - idx representing position for new_bd list
  null_idx <- which(vapply(new_bd, is.null, logical(1)))
  if (length(null_idx) != 0) {
    new_bd <- compact(new_bd)
    rmv_names <- names(bdObj2)[null_idx]
    rmv_sub_id <- vapply(strsplit(rmv_names, "\\."), `[`, character(1), 1)
    idx <- idx[-null_idx]
  }

  ## Update original bdObj with changes
  new_bd <- rbindlist(new_bd)
  for (i in seq_len(nrow(new_bd))) {
    bdObj[idx[i], ] <- new_bd[i, ]
  }

  ## This is gross
  # no, its super gross
  if (length(null_idx) != 0) {

    ## First, get all bdNames, remove those in rmv_names
    bdNames <- do.call(paste, c(bdObj[, nn, with = FALSE], sep = "."))
    rmv_idx <- which(bdNames %in% rmv_names)
    bdObj <- bdObj[-rmv_idx, ]

    ## Now determine which pairs might be left of reduced bdObj
    bdNames2 <- do.call(paste, c(bdObj[, nn[1], with = FALSE], sep = "."))
    rmv_pairs <- which(bdNames2 %in% rmv_sub_id)

    if (length(rmv_pairs)) {
      ll <- length(rmv_sub_id)
      if (ll < 2) {
        msg <- paste0("\n1 observation was deleted during the update process.\n",
                      "This subject has other paired entries in the bdObj dataset.\n",
                      "Would you like to remove their remaining observations?\n",
                      "(may be necessary for paired t-test in bdotsBoot)\n")
      } else {
        msg <- paste0("\n", ll, " observations were deleted during the update process.\n",
                     "These subjects have other paired entires in the bdObj dataset.\n",
                     "Would you like to remove their remaining observations?\n",
                     "(may be necessary for paired t-test in bdotsBoot)\n")
      }

      cat(msg)
      corr_resp <- FALSE
      while (!corr_resp) {
        resp <- readline("Remove all associated observations? (Y/n): ")
        if (resp  %in% c("Y", "n")) {
          corr_resp <- TRUE
        } else {
          cat("Please enter 'Y' or 'n'\n")
        }
      }
      if (resp == "Y") {
        bdObj <- bdObj[-rmv_pairs, ]
      }
    }
  }
  bdObj
}




## A lot of bits of these can be abstracted to functions
# that way we can just jitter as often as we want (or w/e else)
# note, this will also be found in bdotsFit too
# Perhaps throw in  option to "Trash" an observation
## Have to use nn, since it's used
bdUpdate <- function(bdo) {

  # class(bdo) <- c("bdotsObj", class(bdo))
  # x <- X[[nn]]
  # attr(bdo, "X") <- x

  ## Getting curve function
  bdCall <- attr(bdo, "call")
  nn <- c(eval(bdCall[['subject']]), eval(bdCall[['group']])) # this is split vars!
  time <- attr(bdo, "time")
  rho <- attr(bdo, "rho")
  crvFun <- curve2Fun(bdCall[['curveType']])
  fit <- bdo[['fit']][[1]]

  ## This is a terrible work around
  #  to get correct names in x for curveFitter
  #x <- setDT(attr(bdo, "X")$X)
  x <- getSubX(bdo)
  set(x, j = c("y", "time"),
      value = x[,c(bdCall[['y']], bdCall[['time']]), with = FALSE]) # this needs to be subset


  plot(bdo, gridSize = 1)
  oldPars <- printRefitUpdateInfo(bdo)

  accept <- FALSE
  while (!accept) {

    ## Maybe add in future ability to change row
    rf_msg <- paste0("\nActions:\n",
                     "1) Keep original fit\n",
                     "2) Jitter parameters\n",
                     "3) Adjust starting parameters manually\n",
                     "4) Remove AR1 assumption\n",
                     "5) See original fit metrics\n",
                     "6) Delete subject")
    cat(rf_msg)
    resp <- NA
    while (!(resp %in% 1:6)) {
      resp <- readline("Choose (1-6): ")
    }
    if (resp == 1) {
      accept <- TRUE
      break
    } else if (resp == 2) {
      newPars <- jitter(coef(fit))
    } else if (resp == 3) {
      newPars <- oldPars
      cat("Press Return to keep original value\n")
      for (pname in names(oldPars)) {
        cat("Current value:\n")
        print(oldPars[pname])
        tmpval <- readline(paste0("New value for ", pname, ": "))
        if (!is.na(as.numeric(tmpval))) {
          newPars[pname] <- tmpval
        } else if (tmpval != "") {
          warning("Invalid entry, keeping old value")
        }
      }
      class(newPars) <- "numeric"
    } else if (resp == 4) {
      rho <- 0
      newPars <- oldPars
    } else if (resp == 5) {
      printRefitUpdateInfo(bdo)
      next
    } else if (resp == 6) {
        corr_resp <- FALSE
        while (!corr_resp) {
          dd <- readline("Delete observation? (Y/n): ")
          if (dd  %in% c("Y", "n")) {
            corr_resp <- TRUE
          } else {
            cat("Please enter 'Y' or 'n'\n")
          }
        }
        if (dd == "Y") {
          bdo <- NULL
          break
        }
        next
     }

    ## repeat of inside of bdotsFit (make this a function? - yes)
    # then need to update to bdFit as well
    ## Needs splitVars and datVarNames
    result <- bdotsFitter(dat = x, curveType = crvFun, rho = rho, params = newPars,
                          splitVars = nn, datVarNames = bdCall, numRefits = 5)

    new_bdo <- copy(bdo)
    new_bdo$fit <- result$fit
    new_bdo$R2 <- result[['R2']]
    new_bdo$AR1 <- (result[['fitCode']] < 3)
    new_bdo$fitCode <- result[['fitCode']]

    both_bdo <- rbindlist(list(bdo, new_bdo))
    attributes(both_bdo) <- attributes(bdo)
    plot(both_bdo, gridSize = "refit")

    cat("Refit Info:\n")
    printRefitUpdateInfo(new_bdo)
    keep <- readline("Keep new fit? (y/n): ")
    if (length(grep("y", keep, ignore.case = TRUE)) == 1) {
      bdo <- new_bdo
      accept <- TRUE
    } else {
      plot(bdo, gridSize = 1)
    }
  }
  bdo
}


printRefitUpdateInfo <- function(bdo) {
  bdCall <- attr(bdo, "call")
  rho <- attr(bdo, "rho")
  r2 <- round(bdo[['R2']], 3)
  ar1 <- bdo[['AR1']]
  fc <- bdo[['fitCode']]
  fit <- bdo[['fit']][[1]]
  subname <- bdo[, eval(bdCall[['subject']]), with = FALSE]
  msg <- paste0("Subject: ", subname, "\nR2: ", r2, "\nAR1: ",
                as.logical(ar1), "\nrho: ", rho,
                "\nfitCode: ", fc, "\n\n")
  msg <- c(msg, "Model Parameters:\n")
  cat(msg)
  print(oldPars <- coef(fit))
}


## Delete subjects with specified fitCode value
#' bdots Remove Function
#'
#' Remove observations with a specified fitCode and optionally all pairs
#'
#' @param bdObj bdots object
#' @param fitCode min fitCode to remove. Default is 6, which removes all subjects with NULL fits (fitCode = 5 would remove 5 and 6)
#' @param removePairs Boolean. Remove subject pairs is one of pair is removed.
#' Default is TRUE to retain paired t-test
#'
#' @details This function is used to remove all bdots observations with a fit code
#' equal to or larger than the argument passed to \code{fitCode} without refitting.
#' If \code{removePairs == TRUE}, all entries for a subject will be removed if their
#' fit failed in any of the groups in which they were a member
#'
#' @export
bdRemove <- function(bdObj, fitCode = 6L, removePairs = TRUE) {

  ## Checking for decency
  fitCode <- min(6L, max(fitCode, 1L))
  idx <- bdObj[['fitCode']] >= fitCode

  if (removePairs) {
    sub <- attr(bdObj, "call")[['subject']]
    subRmv <- bdObj[[sub]][idx]
    idx <- bdObj[[sub]] %in% subRmv
  }

  bdObj[!idx, ]
}







