
bdObj <- res2

plotFits2 <- function(bdObj, fitCode, gridSize = NULL, ...) {

  bdCall <- attr(bdObj, "call")
  splitVars <- c(bdCall$subject, eval(bdCall$group))
  y <- bdCall$y # what observed values are in X
  if (!is.character(y)) stop("Error 123")
  time <- attr(bdObj, "time")



  X <- attr(bdObj, "X")$X
  dfname <- deparse1(bdCall$data)
  if (is.null(X) & exists(dfname)) {
    X <- get(dfname)
  } else if (is.null(X) & !exists(dfname)) {
    stop("Cannot find dataset used to construct bdObj, please pass as argument")
  }

  Xs <- split(X, by = splitVars)

  old_gg <- theme_set(theme_bw())
  on.exit(theme_set(old_gg))

  for (i in seq_len(nrow(bdObj))) {
    code <- bdObj[i, ]$fitCode
    r2 <- round(as.numeric(bdObj[i, ]$R2), 3)
    if (code == 6) next
    obs <- unlist(bdObj[i, splitVars, with = FALSE])
    obs2 <- paste(obs, collapse = ".")

    dtt <- data.table(time = time)
    dtt$Observed <- Xs[[obs2]][[y]]
    dtt$Fit <- fitted.values(bdObj[i, ]$fit[[1]])
    dtt <- melt(dtt, id.vars = 'time')

    ## Janky fix for update
    if (gridSize == "refit") {
      if (i == 1) {
        title <- paste0("Original Fit\nfitCode = ", code, ", R2 = ", r2)
      } else {
        title <- paste0("Updated Fit\nfitCode = ", code, ", R2 = ", r2)
      }
    } else {
      title <- paste0(paste0(obs, collapse = " "), "\nfitCode = ", code, ", R2 = ", r2)
    }

    test <- melt(dtt, id.vars = 'time')

    ggplot(dtt, aes(x=time, y=value, fill = variable, color = variable)) +
      geom_line(size = 1) +
      #geom_line(linetype = "dashed", aes(time, Observed), size = 1, color = 'cornflowerblue') +
      ylab(y) + xlab(attr(bdObj, "call")$time) +
      ggtitle(title)

  }
}

library(microbenchmark)

microbenchmark(
  seq(1e5),
  seq_len(1e5),
  seq.int(1e5)
)
