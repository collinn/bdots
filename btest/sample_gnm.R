

Logistic <- function(x) {
  list(predictors = list(mini = 1, peak = 1, slope = 1, cross = 1),
       variables = list(substitute(x)),
       term = function(predLabels, varLabels) {
         paste(predLabels[1], " + (", predLabels[2],"-", predLabels[1], ")/(1 + exp(4*",
               predLabels[3], "*(", predLabels[4], "-", varLabels[3], ")/(",
               predLabels[2], "-", predLabels[1], ")))")
       },
       start = function(theta) {
         theta[1] <- 0
         theta[2] <- 0.81355932
         theta[3] <- 0.00079449
         theta[4] <- 776
       })
}
class(Logistic) <- "nonlin"

fit <- gnm(Fixations ~ Logistic(dat$Time), data = dat, family = "binomial")


Logistic <- function(x){
  list(predictors = list(Asym = 1, xmid = 1, scal = 1),
       variables = list(substitute(x)),
       term = function(predLabels, varLabels) {
         paste(predLabels[1], "/(1 + exp((", predLabels[2], "-",
               varLabels[3], ")/", predLabels[3], "))")
       },
       start = function(theta){
         theta[3] <- 1
         theta
       })
}
class(Logistic) <- "nonlin"

debugonce(gnm)

fit <- gnm(Fixations ~ Logistic(Time), data = dat, family = "binomial")


##### with gnls

fit <- gnls(eval(ff), start = params, data = dat,
            correlation = corAR1(rho),
            control = gnlsControl(maxIter = 0, nlsMaxIter = 0, msMaxIter = 0, returnObject = TRUE,
                                  optimMethod="L-BFGS-B", lower = c(0,0,0,0), upper = c(1, 1, 1, 2000)))

fit2 <- gnls(eval(ff), start = params, data = dat,
            correlation = corAR1(rho),
            control = gnlsControl(maxIter = 0, nlsMaxIter = 0, msMaxIter = 0, returnObject = TRUE))


## Try with perfect logistic
library(eyetrackSim)

newfix <- eyetrackSim:::logistic_f(params, dat$Time)
dat2 <- copy(dat)
dat2$Fixations <- newfix


newpars <- params
for (i in seq_along(newpars)) newpars[i] <- jitter(params[i])

fit <- gnls(eval(ff), start = newpars, data = dat2,
             correlation = corAR1(rho),
             control = gnlsControl(maxIter = 0, nlsMaxIter = 0, msMaxIter = 0, returnObject = TRUE))

plot(x = dat2$Time, y = fitted.values(fit), type = 'l')
lines(x = dat2$Time, y = dat2$Fixations, col = 'blue')
