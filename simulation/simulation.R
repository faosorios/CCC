## reading R codes (files must be at the working directory)
source("ccc.R")
source("ccc.influence.R")
source("ccc.simul.R")

## load 'heavy' package required to simulate multivariate normal deviates
library(heavy)

## setting parameters
mu    <- c(0,0)
Sigma <- matrix(c(1,.95,.95,1), ncol = 2)

## creating Table 2..
# seeds: 5, 23, 53, 167, 239, 347, 433 (OMG! all are prime numbers...)

## 1st column: delta = 0.5
set.seed(5)
out <- simul.ccc(mean = mu, Sigma = Sigma, delta = 0.5)
out[c("percentage","cutoffs","speed")] # just 'percentages' are reported at Table 2

## 2nd column: delta = 1.0
set.seed(23)
out <- simul.ccc(mean = mu, Sigma = Sigma, delta = 1.0)
out[c("percentage","cutoffs","speed")] # just 'percentages' are reported at Table 2

## 3rd column: delta = 1.5
set.seed(53)
out <- simul.ccc(mean = mu, Sigma = Sigma, delta = 1.5)
out[c("percentage","cutoffs","speed")] # just 'percentages' are reported at Table 2

## 4th column: delta = 2.0
set.seed(167)
out <- simul.ccc(mean = mu, Sigma = Sigma, delta = 2.0)
out[c("percentage","cutoffs","speed")] # just 'percentages' are reported at Table 2
x <- out$x # saving a 'typical' dataset for posterior analyses

## 5th column: delta = 2.5
set.seed(239)
out <- simul.ccc(mean = mu, Sigma = Sigma, delta = 2.5)
out[c("percentage","cutoffs","speed")] # just 'percentages' are reported at Table 2

## 6th column: delta = 3.0
set.seed(347)
out <- simul.ccc(mean = mu, Sigma = Sigma, delta = 3.0)
out[c("percentage","cutoffs","speed")] # just 'percentages' are reported at Table 2

## 7th column: delta = 3.5
set.seed(347)
out <- simul.ccc(mean = mu, Sigma = Sigma, delta = 3.0)
out[c("percentage","cutoffs","speed")] # just 'percentages' are reported at Table 2

## Figure 4
par(pty = "s")
plot(x, xlim = c(-3.5,3.5), ylim = c(-3.5,3.5), xlab = "X1", ylab = "X2")
abline(c(0,1), col = "red", lty = 2, lwd = 2)
obs <- c(1,15,27,46,60)
text(x[obs,1], x[obs,2], as.character(obs), pos = 3)

## fitting postulated model
fm <- fit.ccc(x)

# Fig 5.a
z  <- influence.ccc(fm, method = "normal")
plot(z, idn = 4)

# Fig 5.b
z  <- influence.ccc(fm, method = "conformal")
plot(z, idn = 4)

# Fig 5.c
z  <- influence.ccc(fm, method = "FI")
plot(z)

# Fig 5.d
z  <- influence.ccc(fm, method = "SI")
plot(z)

## Fig 6
nobs <- 100 # nrow(x)
rhoc <- rep(0, nobs)
accu <- rep(0, nobs)
r    <- rep(0, nobs)
for (i in 1:nobs) {
  f0 <- fit.ccc(x, subset = -i)
  rhoc[i] <- f0$ccc
  accu[i] <- f0$accuracy
  r[i]    <- f0$precision
}
par(mfrow = c(3,1))
plot(rhoc, type = "b", ylim = c(.935,.965), ylab = "CCC estimate")
abline(h = fm$ccc, col = "red", lty = 2, lwd = 2)
text(1, rhoc[1], as.character(1), pos = 3)
plot(accu, type = "b", ylim = c(.9991,1), ylab = "accuracy estimate")
abline(h = fm$accuracy, col = "red", lty = 2, lwd = 2)
obs <- c(15,46)
text(obs, accu[obs], as.character(obs), pos = 1)
obs <- c(27,60)
text(obs, accu[obs], as.character(obs), pos = 3)
plot(r, type = "b", ylim = c(.935,.965), ylab = "precision estimate")
abline(h = fm$precision, col = "red", lty = 2, lwd = 2)
text(1, r[1], as.character(1), pos = 3)
