## Actual bdots fitting function, not exported to user

## We could almost assume that this takes subjectDat <- dat[subject == w/e, ]

## test environment
#load(file = "~/packages/bdots/data/test_env.RData")
#dat <- tdat[subject == "2", ]



## bdotsFitter takes a dat of a single individual and their curve type
# It will be simple and robust, and any changes at the subject level
# can be done here, in one place

## Actually, this function CAN be exported. Once the person knows which thing
# failed, they can inspect it OR they can take names from that list and just pass
# their own subset portion of the dataset to do it themselves.
# The modularity of this is coming along nicely, I think
for(i in seq_along(newdat)) {
  bdotsFitter(newdat[[i]])
}



bdotsFitter <- function(dat, curveType, concave, cor, rho, refits = FALSE, verbose, params = NULL) {



  ## This will need to be handled differently, as dat may not contain
  # y, time, etc
  if(any(dat[, .N, by = .(time)]$N > 1)) {
    warning("Some subjects have multiple observations for unique time. These will be averaged")
    dat[, y := mean(y), by = .(time)]
    dat <- unique(dat, by = c("time"))
  }

  ## Well, determine what curveType, but for now, just  gauss
  # which means I need to make estLogistic
  # also the order of these variables (and in curveFitter and elsewhere) is fucked
  #fit <- estDgaussCurve(dat, rho = 0.9, get.cov.only = FALSE, cor = TRUE, jitter = 0, conc = TRUE)
  fit <- estLogisticCurve(dat, rho = 0.9, get.cov.only = FALSE, cor = TRUE, refits = 4)
  #str(fit)

  #fit <- estLogisticCurve(dat, rho, cor = FALSE, jitter = 5)
  ## Assign fail codes next <- this can be done in bdotsFit, as we may not typically need these (or would we?)
  # I guess we need to compute other shit too, like r2 or w/e
  # if fit[1] null and cor = false, then w/e
  # fit [1] not null, cor = t , w/e,
  # etc


  ## Now determine what exactly I need for bootstrapping step or for reporting
  ## also need to revist  verbose
  ## This also potentially could be done in bdotsFit?
  ## It could, but do it here, for verbose option. There is no need to compute anything else
  # we can have helper functions that do that for us, which benefits us in two ways
  # 1 - summaries are divorced from these functions
  # 2 - keeps output tidy for use with other things
  # pretty much, that's it. Anyone using this can grab a particular subject
  # instance and just do their own assessment OR, even better, it makes it easier
  # to 'swap out' entries in refitting OR manually, i.e, for one specific subject
  # they can recall bdotsFitter with different jitter or rho, if necessary
  SSE <- sum(resid(fit[[1]])^2)
  SSY <- sum((dat$y - mean(dat$y))^2)
  R2 <- 1 - SSE/SSY

  ## possibly move this to bdotsFit, as the simpler this function return value is,
  # the easier it is for user to update themselves.
  ## Fail Code ==
  # + 1 - 0.8 < R2 < 0.95
  # + 2 - R2 < 0.8
  # + 3 - AR1 == FALSE
  ## This gives
  # 0 - AR1 TRUE and R2 > .95
  # 1 - AR1 TRUE and 0.8 < R2 < 0.95
  # 2 - AR1 TRUE and R2 < 0.8
  # 3 - AR1 FALSE and R2 > .95
  # 4 (3 + 1) - AR1 FALSE and 0.8 < R2 < 0.85
  # 5 (3 + 2) - AR1 FALSE and R2 < 0.8
  fitCode <- 3*(!fit[["cor"]]) + 1*(R2 < 0.95)*(R2 > 0.8) + 2*(R2 < 0.8)

  ## Practice verbose
  if(TRUE)

  ## Don't actually need fit[[2]], since we can check if failCode > 3
  ## Returns list(model_fit, AR1cor (boolean), R2, failCode)
    ## Consider also returning parameters start position?
  list(fit = fit[["fit"]], R2 = R2, fitCode = fitCode)
}




