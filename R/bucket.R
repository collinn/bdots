# values of sig are TRUE/FALSE
#realsig <- sig
# sig <- realsig
# sig <- c(rep(TRUE, 50), rep(FALSE, 59), sig[110:length(sig)])

## Determines time regions where difference significant
# sig are boolean values
bucket <- function(sig, time) {
  rr <- rle(sig)
  nr <- sum(rr$values) # number of sig regions
  if (nr == 0) return(NULL)
  rr$idx <- cumsum(rr$lengths)
  mm <- rr$values*(rr$idx) + (1 - rr$values)*(rr$idx + 1)

  ## Alt condition if first region true
  if (rr$values[1] == 1) {
    mm <- mm[c(length(mm), 1:(length(mm)-1))]
    mm[1] <- 1
  }
  matrix(time[mm], ncol = 2, byrow = TRUE)
}

