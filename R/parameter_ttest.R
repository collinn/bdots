
# groups are ExpectNoise, NoiseBoth, Quiet

## Group would be trial, the column name
# vals would be the values in group. For vals,
# if group has two values, it can be empty
# if group has more than two, it will do all
# otherwise, it has to be a length two or more of values located
# within the group column
# bdObj <- fit_cohort
# group <- "trial"
# vals <- c("ExpectNoise", "NoiseBoth", "Quiet")

#' Parameter t-test
#'
#' Perform t-test on curve parameters of bdotsFit object
#'
#' @param bdObj Object of class \code{bdObj}
#' @param group Length one character of grouping column in which to perform t-test
#' @param vals Character vector of values within grouping column in which to perform the
#' test. If \code{NULL}, it will do all pairwise tests
#'
#' @return List of t-test results of class \code{bdotsPars_ttest}
#'
#' @details Performs pairwise t-test. Currently only tests at alpha = 0.95. Also
#' currently only allows t-test within single grouping column. Ability to test
#' across grouping columns to come later
#' @examples
#' \dontrun{
#' res <- bdotsFit(data = cohort_unrelated,
#'                 subject = "Subject",
#'                 time = "Time",
#'                 y = "Fixations",
#'                 group = c("Group", "LookType"),
#'                 curveType = doubleGauss(concave = TRUE),
#'                 cor = TRUE,
#'                 numRefits = 2,
#'                 cores = 0,
#'                 verbose = FALSE)
#' tstats <- parTest(res, group = "LookType", vals = c("Cohort", "Unrelated_Cohort"))
#' }
#' @importFrom utils combn
#' @export
parTest <- function(bdObj, group, vals = NULL) {

  #### This whole section needs to be logic validated
  if (is.null(vals)) {
    vals <- names(table(bdObj[[group]]))
  }
  valid_vals <- vals %in% names(table(bdObj[[group]]))
  nvals <- length(vals)
  if (nvals < 2) {
    stop(paste("Must have two or more values to subset", group))
  }
  #### End section that needs to be logic validated

  idx <- bdObj[[group]] %in% vals
  bdgroups <- split(bdObj[idx, ], by = group)
  coef_mats <- lapply(bdgroups, coef)

  ## need all unique pairwise elements in vector 1:n
  tidx <- utils::combn(seq_len(nvals), 2)
  ttests <- apply(tidx, 2, function(x) {
    sub <- lapply(bdgroups[x], function(y) {
      ss <- attr(y, "call")$subject
      y[[ss]]
    })
    isPaired <- Reduce(identical, sub)

    pm <- lapply(bdgroups[x], coef)
    xx <- pm[[1]]
    yy <- pm[[2]]
    pars <- vector("list", length = ncol(xx))
    pars <- setNames(pars, colnames(xx))
    for (i in seq_along(pars)) {
      pars[[i]] <- t.test(xx[, i], yy[, i], paired = isPaired) # don't know if paired
    }
    attr(pars, "isPaired") <- isPaired
    pars
  })

  names(ttests) <- apply(tidx, 2, function(x) {
    paste(vals[x], collapse = " vs. ")
    })

  make_tabs <- function(xx) {
    res <- lapply(names(xx), function(nn) {
      tt <- xx[[nn]]
      est <- ifelse(length(tt$estimate) == 2, diff(tt$estimate), tt$estimate)
      vv <- signif(c(est, tt$conf.int, tt$statistic, tt$p.value), digits = 2)
      vv <- matrix(vv, nrow = 1)
      colnames(vv) <- c("diff", "95L", "95U", "tstat", "pval")
      if (vv[, "pval"] < 0.05) {
        nn <- paste0(nn, "*")
      }
      rownames(vv) <- nn
      vv
    })
    res <- Reduce(rbind, res)
    attr(res, "isPaired") <- attr(xx, "isPaired")
    res
  }

  if (length(ttests) == 1) {
    tabs <- list(make_tabs(ttests[[1]]))
    names(tabs) <- names(ttests)
  } else {
    tabs <- lapply(ttests, make_tabs)
  }

  structure(.Data = tabs, class = "bdotsPars_ttest")
}

#' Print Parameter Test Summary
#'
#' Print Parameter Test Summary
#'
#' @param x object to be printed
#' @param ... not used
#'
#' @details That's pretty much it. This is a print method, so there is likely
#' not much need to call it directly
#' @export
print.bdotsPars_ttest <- function(x, ...) {

  ## Figure out header, i.e., degrees of freedom + isPaired
  for (i in seq_along(x)) {
    isPaired <- attr(x[[i]], "isPaired")
    if (isPaired) {
      cat(paste0(names(x)[i], " (paired) :\n"))
    } else {
      cat(paste0(names(x)[i], " (unpaired):\n"))
    }
    attr(x[[i]], "isPaired") <- NULL
    print(x[[i]])
    attr(x[[i]], "isPaired") <- isPaired
    cat("\n\n")
  }
  return(invisible(x))
}

# rr <- parTest(fit_cohort, "trial")
# rr
# qq <- parTest(res.b, "Group")
# qq
