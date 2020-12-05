

# dat <- data.table(cohort_unrelated)
# dat <- dat[Subject == 1 & LookType == "Cohort" & Group == 50, ]
# xx <- dat
# time <- "Time"
# y <- "Fixations"
#
# yy <- dat[[y]]
# tt <- dat[[time]]
# degree <- 5
# raw <- TRUE
####

polynomial <- function(dat, y, time, degree, raw = TRUE, params = NULL, ...) {

  polyPars <- function(dat, y, time, degree, raw, ...) {
    yy <- dat[[y]]
    tt <- dat[[time]]
    pp <- lm(yy ~ poly(tt, degree = degree, raw = raw))
    rr_names <- c(paste0("beta", seq(degree + 1L)))
    rr <- setNames(coef(pp), rr_names)
  }

  if (is.null(params)) {
    params <- polyPars(dat, y, time, degree, raw)
  } else {
    if (length(params) != degree + 1) stop("poly requires that number of params matches degree + 1 (for intercept)")
    if (!all(names(params) %in% paste0("beta", seq(degree)))) {
      stop("polynomial parameters must be labeled beta1, ..., beta[degree + 1]")
    }
  }

  ## This needs to be appended to dat. It changes dat. Is this dangerous? Probably
  # mm <- model.matrix(~poly(dat[[time]], degree = degree, raw = raw) - 1, data = dat)
  # mm_names <- paste0("bdpoly_", letters[seq(degree)])
  # colnames(mm) <- mm_names
  # dat[, (mm_names) := lapply(mm_names, function(x) mm[, x])]

  time_names <- paste0("I(Time^", seq(degree + 1L) - 1L , ")")

  ff <- paste(names(params), time_names, sep = "*", collapse = "+")
  ff <- str2lang(ff)
  y <- str2lang(y)
  ff <- bquote(.(y) ~ .(ff))
  return(list(formula = ff, params = params))
}

# res <- polynomial(dat, y, time, degree, raw)
# ff <- res[['formula']]
# params <- res[['params']]
# fit <- curveFitter(dat, ff, params, rho, numRefits)
#
#
#
#
# fit <- gnls(eval(ff), data = dat, start = params, correlation = corAR1(rho))
# plot(x, fitted.values(fit), type = 'l')
# lines(x,dat$Fixations,col='red')
#
#
#
#
#
#
# dat2 <- copy(dat)
#
# ww <- function() {
#   dat2[, (mm_names) := lapply(mm_names, function(x) mm[, x])]
# }
# rr_names <- names(rr)
#
# time_names <- paste0("I(Time^", seq(degree), ")")
# ff <- paste(rr_names, time_names, sep = "*", collapse = "+")
# ff <- paste0("Fixations~", ff)
# ff <- as.formula(ff)
# fit <- gnls(ff, data = dat, start = rr,
#             correlation = corAR1(0.9))
#
#
# timett <- "Time"
# f_bod <- ff[[3]]
# f_args <- c(rr_names, timett)
# f_args <- setNames(as.pairlist(rep("", length(f_args))), f_args)
#
# new_f <- eval(call("function", f_args, f_bod))
#
# myargs <- as.list(rr)
# myargs$Time <- dat$Time
#
# test1 <- do.call(new_f, myargs)
