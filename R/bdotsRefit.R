#' Refit Observations Returned from bdotsFit
#'
#'
#'
#' @param bdObj An object of class 'bdotsObj' returned from \code{bdotsFit}
#' @param fitCode A length one integer indicating observations to refit. See Details
#' @param quickRefit Boolean indicating if a quick refit should be used. If TRUE,
#' rather than prompting the user for adjustments for each observation, \code{bdotsReft}
#' will jitter the parameters of all observations indicated by \code{fitCode} and attempt
#' to refit. Between the original and the refitted curve, this will place priority on
#' the higher \code{fitCode}. If these are equal, R2 will take precedence. Otherwise,
#' the original fit will be kept.
#' @param numRefits Integer indicating the number of refit attempts after jittering
#' parameters, either with quickRefit or when done individually
#' @param ... not used
#'
#' @return Returns bdObj with updated fits
#'
#' @details fitCode indicates lower bound on observations to refit. For example,
#' if \code{fitCode = 4}, \code{bdotsRefit} will prompt user to refit all
#' observations with fitCode = 4, 5, 6. The \code{quickRit} option will attempt
#' to jitter and refit all observations selected by \code{fitCode}. Otherwise, the
#' user will be prompted through a menu to individually refit observations
#' @import data.table
#' @export
bdotsRefit <- function(bdObj, fitCode = 1L, quickRefit = FALSE, numRefits = 2L, ...) {

  if (is.null(fitCode)) fitCode <- 1L


  if (is.null(attr(bdObj, "X")$X)) {
    stop("Dataset must be provided")
  }

  ## These uniquely identify each fit
  bdCall <- attr(bdObj, "call")
  nn <- c(eval(bdCall[['subject']]), eval(bdCall[['group']]))

  ## Because of factors, which I should really think about removing
  fitcode <- fitCode
  idx <- which(bdObj$fitCode >= fitcode)
  if (length(idx) == 0L) {
    message(paste0("All observations fitCode greater than ",
                   fitcode, ". Nothing to refit :)"))
    return(bdObj)
  }
  bdObj2 <- split(bdObj[idx, ], by = nn)

  if (quickRefit) {
    ## Oddly, I get errors running parLapply not lapply. Will investigate
    # cores <- attr(bdObj, "call")$cores
    # if (cores == 0) cores <- detectCores()/2
    # cl <- makePSOCKcluster(cores)
    # invisible(clusterEvalQ(cl, library(bdots)))
    # new_bd <- parLapply(cl, bdObj2, bdQuickRefit, numRefits)
    # stopCluster(cl)
    new_bd <- lapply(bdObj2, bdQuickRefit, numRefits)
  } else {
    new_bd <- lapply(bdObj2, bdUpdate, numRefits)
  }

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
        if (resp  %in% c("Y", "n"))
          corr_resp <- TRUE
        else
          cat("Please enter 'Y' or 'n'\n")
      }
      if (resp == "Y")
        bdObj <- bdObj[-rmv_pairs, ]
    }
  }
  bdObj
}


bdQuickRefit <- function(bdo, numRefits) {



  if (bdo$fitCode != 6L) {
    njitter <- 5L
    newPars <- coef(bdo[['fit']][[1]])
  } else {
    newPars <- NULL
  }

  while(numRefits != 0L) {

    ## jitterbug
    if (bdo$fitCode != 6L) {
      for (i in seq_along(newPars)) {
        newPars[i] <- jitter(newPars[i], factor = njitter)
      }
      njitter <- njitter + 1L
    }

    new_bdo <- bdRefitter(bdo, params = newPars)
    if (new_bdo$fitCode < bdo$fitCode) {
      return(new_bdo)
    } else if (new_bdo$fitCode == bdo$fitCode) {
      v <- ifelse(is.na(new_bdo$R2), 0, new_bdo$R2) - ifelse(is.na(bdo$R2), 0, bdo$R2)
      if (v > 0.05) return(new_bdo)
    }
    numRefits <- numRefits - 1L
  }
  return(bdo)
}


## Given a single bdotsFit observation, computes a refit
bdRefitter <- function(bdo, numRefits = 0L, rho = NULL, params = NULL, ...) {
  if (nrow(bdo) != 1L) stop("bdRefitter can only take a single observation")
  bdCall <- attr(bdo, "call")
  nn <- c(eval(bdCall[['subject']]), eval(bdCall[['group']])) # this is split vars!
  if (is.null(rho)) rho <- attr(bdo, "rho")
  crvFun <- curve2Fun(bdCall[['curveType']])

  x <- getSubX(bdo)
  set(x, j = c("y", "time"),
      value = x[,c(bdCall[['y']], bdCall[['time']]), with = FALSE])

  new_bdo <- bdotsFitter(dat = x, curveType = crvFun, rho = rho,
                         params = params, splitVars = nn, datVarNames = bdCall,
                         numRefits = numRefits)
  prob <- tryCatch(attributes(new_bdo) <- attributes(bdo), error = function(e) 2)
  if (is.numeric(prob)) browser()
  new_bdo
}



## A lot of bits of these can be abstracted to functions
## Need to be able to handle fitCode == 6
bdUpdate <- function(bdo, numRefits) {

  ## This will attempt to fit without. If a fit is made, then fitCode
  # will no longer be 6. In that case, bdUpdate_NULL will call bdUpdate
  # which will then skip this step
  if (bdo$fitCode == 6) {
    bdo <- bdUpdate_NULL(bdo, numRefits)
    return(bdo)
  }

  rho <- attr(bdo, "rho")
  plot(bdo, gridSize = 1)

  ## Adding newPars here allows jitter to bounce around more
  oldPars <- newPars <- printRefitUpdateInfo(bdo)

  ## Increase jitter each  time it occurs
  njitter <- 5L
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
      for (i in seq_along(newPars)) {
        newPars[i] <- jitter(newPars[i], factor = njitter)
      }
      njitter <- njitter + 1L
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
      next # reset while loop
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
        next # reset while loop
    }

    new_bdo <- bdRefitter(bdo, rho, numRefits, params = newPars)
    both_bdo <- structure(.Data = list(bdo, new_bdo),
                          class = "bdObjList")

    both_bdo <- rbindlist(both_bdo)
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

## Prints info for update
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


## This should probably move to a different file.
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
#' If \code{removePairs = TRUE}, all entries for a subject will be removed if their
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
  ## Should I also subset the data matrix?
  bdObj[!idx, ]
}




bdUpdate_NULL <- function(bdo, numRefits) {
  plot(bdo, gridSize = 1)

  rounds <- 1L
  while (TRUE) {

    if (rounds == 1L) {
      rf_msg <- paste0("\nActions:\n",
                       "1) Adjust starting parameters manually\n",
                       "2) Delete subject")
      cat(rf_msg)
      resp <- NA
      while (!(resp %in% 1:2)) {
        resp <- readline("Choose (1-2): ")
      }
      rounds <- rounds + 1L
    } else {
      rf_msg <- paste0("\nActions:\n",
                       "1) Adjust starting parameters manually\n",
                       "2) Print previous parameter attempt\n",
                       "3) Update previous parameter attempt\n",
                       "4) Delete subject")
      cat(rf_msg)
      resp <- NA
      while (!(resp %in% 1:4)) {
        resp <- readline("Choose (1-4): ")
      }
    }

    if (resp == 1) {
      parnames <- attributes(attr(bdo, "formula"))[['parnames']]
      newPars <- setNames(rep(NA, length(parnames)), parnames)
      for (pname in names(newPars)) {
        tmpval <- NA
        while (is.na(as.numeric(tmpval))) {
          tmpval <- readline(paste0("New value for ", pname, ": "))
          if (!is.na(as.numeric(tmpval)))
            newPars[pname] <- as.numeric(tmpval)
          else
            warning("Invalid entry, please enter numeric value")
        }
      }

      new_bdo <- bdRefitter(bdo, rho = 0.9, numRefits, params = newPars)

      if (new_bdo$fitCode == 6) {
        cat("Fit unsuccessful. Plotting curve from your input parameters in red.\n",
            "Use this to adjust parameter estimates accordingly\n")

        curveFun <- makeCurveFun(bdo)
        Time <- attr(bdo, "time")
        TimeName <- attr(bdo, "call")[['time']]
        parList <- as.list(newPars)
        parList[[TimeName]] <- Time
        suggestFit <- do.call(curveFun, parList)
        plot(bdo, gridSize = 1)
        lines(x = Time, y = suggestFit, lwd = 2, col = 'tomato')
        next
      } else {
        cat("Fit success! Moving to standard refit options for current observation\n")
        readline("Press Enter to continue")
        bdo <- bdUpdate(new_bdo, numRefits)
        break
      }
    } else if (resp == 2) {
      print(newPars)
    } else if (resp == 3) {
      oldPars <- newPars
      cat("Press Enter to keep current value\n")
      for (pname in names(oldPars)) {
        cat("Current value:\n")
        print(oldPars[pname])
        tmpval <- NA
        while (is.na(as.numeric(tmpval))) {
          tmpval <- readline(paste0("New value for ", pname, ": "))
          if (!is.na(as.numeric(tmpval)))
            newPars[pname] <- as.numeric(tmpval)
          else if (tmpval == "") {
            newPars[pname] <- oldPars[pname]
            break
          } else
            warning("Invalid entry, please enter a numeric value")
        }
      }
      new_bdo <- bdRefitter(bdo, rho = 0.9, numRefits, params = newPars)
      if (new_bdo$fitCode == 6) {
        cat("Fit unsuccessful. Plotting curve from your input parameters in red.\n",
            "Use this to adjust parameter estimates accordingly\n")

        curveFun <- makeCurveFun(bdo)
        Time <- attr(bdo, "time")
        TimeName <- attr(bdo, "call")[['time']]
        parList <- as.list(newPars)
        parList[[TimeName]] <- Time
        suggestFit <- do.call(curveFun, parList)
        plot(bdo, gridSize = 1)
        lines(x = Time, y = suggestFit, lwd = 2, col = 'tomato')
        next
      } else {
        cat("Fit success! Moving to standard refit options for current observation\n")
        readline("Press Enter to continue")
        bdo <- bdUpdate(new_bdo, numRefits)
        break
      }
    } else if (resp == 4) {
      corr_resp <- FALSE
      while (!corr_resp) {
        dd <- readline("Delete observation? (Y/n): ")
        if (dd  %in% c("Y", "n"))
          corr_resp <- TRUE
        else
          cat("Please enter 'Y' or 'n'\n")
      }
      if (dd == "Y") {
        bdo <- NULL
        break
      }
      next # reset while loop
    }
  }
  return(bdo) # this will return to orig. bdUpdate call and then return again there
}




