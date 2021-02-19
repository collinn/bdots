
fit <- gnls(Fixations ~ poly(Time, 6), data = xx,
            correlation = corAR1(0.9))

mm <- model.matrix(~poly(Time, 5), data = xx)
colnames(mm) <- c("int", letters[1:5])

fit <- gnls(xx$Fixations ~ int + a + b + c + d + e,
            data = mm, correlation = corAR1(0.9))



x <- xx$Time
y <- xx$Fixations

pp <- lm(y ~ poly(x, 5))
mm <- model.matrix(~poly(Time, 5), data = xx)
colnames(mm) <- c("int", letters[1:5])
rr <- setNames(coef(pp), colnames(mm))[-1]
mm <- cbind(mm, y)
fit <- gnls(y ~ a*x + b*I(x^2) + c*I(x^3) + d*I(x^4) + e*I(x^5), data = data.frame(mm), start = rr,
            correlation = corAR1(0.9),
            control = gnlsControl(nlsTol = .25))

plot(x, fitted.values(fit), type = 'l')
lines(x,y,col='red')

pp1 <- lm(y ~ poly(x, 5, raw = TRUE))
mm1 <- model.matrix(~poly(Time, 5, raw = TRUE), data = xx)
colnames(mm1) <- c("int", letters[1:5])
rr1 <- setNames(coef(pp1)[-1], LETTERS[1:5])
mm1 <- cbind(mm1, y)
fit1 <- gnls(y ~ a*A + b*B + c*C + d*D + e*E, data = data.frame(mm1), start = rr1,
            correlation = corAR1(0.9),
            control = gnlsControl(nlsTol = .001))

plot(x, fitted.values(fit1), type = 'l')
lines(x,y,col='red')
