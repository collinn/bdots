## I do not know where these derivations come from - this needs to be followed up on.
## Otherwise, code modifed from bdots_v_0.1.30
## Fairly confident this takes a vector of t-values, generated from t distributio for test
## Internal function, and is used to get estimate of phi0 which is subsequently
## passed to a modify.alpha function in bdots used to do oleson adjustment for pvalues
## My goodness. This could be internal to a modified p.adjust
## At any rate, it makes sense that these functions all live together here

ar1Solver <- function(y) {
  ## Not sure what motivates these values
  n <- length(y)
  v1 <- sum(y[2:n]*y[1:(n-1)])
  v2 <- n - 2 * sum(y^2) + y[n]^2

  phiEst <- polyroot(c(v1, v2, v1, -n))
  phiEst <- Re(phiEst[abs(Im(phiEst)) < 0.0001]) # Discard complex roots
  phiEst <- phiEst[phiEst >= 0]
  phiEst <- phiEst[ar1ScoreDerivative(phiEst, y) < 0] # Discard roots where not concave
  phiEst <- phiEst[phiEst < 1]
  if(length(phiEst) > 1) {
    phiEst <- phiEst[which.max(ar1Loglikelihood(phiEst, y))]
  } else if (length(phiEst) == 0) {
    warning("AR1 solver could not find a proper root. Using Yule-Walker solver.")
    phiEst <- ar(y, order.max = 1, aic = FALSE)$ar
  }
  phiEst
}

## Is this the score? The derivative of the score?
## Need to ask Jake about these ar1 functions
ar1ScoreDerivative <- function(phi, y) {
  phi2 <- phi ^ 2
  n <- length(y)

  ## Not sure where these values come from either
  v1 <- (1 + phi2) / (1 - phi2) ^ 2 * n
  v2 <- (1 + 3 * phi2) / (1 - phi2)^3 * (2 * sum(y ^ 2) - y[n] ^ 2)
  v3 <- 2 * phi * (3 + phi2) / (1 - phi2) ^ 3 * sum(y[2:n] * y[1:(n - 1)])
  v1 - v2 + v3
}

## Ditto above. Clearly it's normal, though
ar1Loglikelihood <- function(phi, y) {
  n <- length(y)
  v1 <- - n / 2 * log(2 * pi * (1 - phi ^ 2))
  v2 <- y[1] ^ 2 + sum((y[2:n] - phi * y[1:(n - 1)]) ^ 2)
  v3 <- 2 * (1 - phi ^ 2)

  v1 - v2 / v3
}
