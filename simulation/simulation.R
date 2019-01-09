## reading R sources
source("../code/ccc.R")
source("../code/ccc.influence.R")
source("../code/poa.influence.R")
source("ccc.simul.R")
source("poa.simul.R")

## load 'heavy' package required to simulate multivariate normal deviates
library(heavy)

## setting parameters
mu    <- c(0,0)
Sigma <- matrix(c(1,.95,.95,1), ncol = 2)

## Table 3 (seeds: 5, 23, 53, 167, 239, 347, 433)

summary.ccc(Nsize = 500, nobs = 25, mean = mu, Sigma = Sigma)
#      0.5  1.0  1.5  2.0  2.5   3.0   3.5
#C    11.4 33.8 66.4 77.8 86.4  93.4  95.4
#B    11.4 33.8 66.4 77.8 86.4  93.4  95.4
#FI   28.4 74.4 96.0 99.0 99.8 100.0 100.0
#SI   11.4 33.8 66.4 77.8 86.4  93.4  95.4
#both  6.4 27.4 63.2 76.8 86.2  93.4  95.4

summary.ccc(Nsize = 500, nobs = 50, mean = mu, Sigma = Sigma)
#      0.5  1.0  1.5   2.0 2.5   3.0 3.5
#C     9.6 34.2 59.4  77.4  90  96.4  97
#B     9.6 34.2 59.4  77.4  90  96.4  97
#FI   26.0 76.0 96.4 100.0 100 100.0 100
#SI    9.6 34.2 59.4  77.4  90  96.4  97
#both  5.6 28.4 58.0  77.4  90  96.4  97

summary.ccc(Nsize = 500, nobs = 100, mean = mu, Sigma = Sigma)
#      0.5  1.0  1.5   2.0   2.5   3.0 3.5
#C    11.8 30.8 53.8  74.8  92.4  98.4  98
#B    11.8 30.8 53.8  74.8  92.4  98.4  98
#FI   27.0 77.2 97.2 100.0 100.0 100.0 100
#SI   11.8 30.8 53.8  74.8  92.4  98.4  98
#both  8.8 26.8 52.8  74.8  92.4  98.4  98

summary.ccc(Nsize = 500, nobs = 200, mean = mu, Sigma = Sigma)
#      0.5  1.0  1.5 2.0   2.5   3.0 3.5
#C    16.0 42.8 58.2  68  89.6  98.2  99
#B    16.0 42.8 58.2  68  89.6  98.2  99
#FI   26.4 79.2 97.2 100 100.0 100.0 100
#SI   16.0 42.8 58.2  68  89.6  98.2  99
#both  9.6 36.6 57.0  68  89.6  98.2  99

## Table 4 (seeds: 5, 23, 53, 167, 239, 347, 433)

summary.poa(Nsize = 500, nobs = 25, mean = mu, Sigma = Sigma)
#      0.5  1.0  1.5 2.0 2.5 3.0 3.5
#C    24.6 77.6 98.2 100 100 100 100
#B    24.6 77.6 98.2 100 100 100 100
#FI   26.4 81.4 98.6 100 100 100 100
#SI   24.6 77.6 98.2 100 100 100 100
#both 24.6 77.4 98.2 100 100 100 100

summary.poa(Nsize = 500, nobs = 50, mean = mu, Sigma = Sigma)
#      0.5  1.0  1.5 2.0 2.5 3.0 3.5
#C    23.4 74.0 97.8 100 100 100 100
#B    23.4 74.0 97.8 100 100 100 100
#FI   27.0 81.4 98.2 100 100 100 100
#SI   23.4 74.0 97.8 100 100 100 100
#both 22.4 74.0 97.6 100 100 100 100

summary.poa(Nsize = 500, nobs = 100, mean = mu, Sigma = Sigma)
#      0.5  1.0  1.5   2.0 2.5 3.0 3.5
#C    17.2 58.8 94.6  99.8 100 100 100
#B    17.2 58.8 94.6  99.8 100 100 100
#FI   26.8 82.2 98.4 100.0 100 100 100
#SI   17.2 58.8 94.6  99.8 100 100 100
#both 16.8 58.2 94.6  99.8 100 100 100

summary.poa(Nsize = 500, nobs = 200, mean = mu, Sigma = Sigma)
#      0.5  1.0  1.5   2.0 2.5 3.0 3.5
#C    16.8 59.4 89.2  98.8 100 100 100
#B    16.8 59.4 89.2  98.8 100 100 100
#FI   27.4 83.6 98.8 100.0 100 100 100
#SI   16.8 59.4 89.2  98.8 100 100 100
#both 15.8 57.4 89.0  98.8 100 100 100

## 'typical' dataset used in the simulation experiment
set.seed(167)
out <- simul.ccc(Nsize = 500, nobs = 100, mean = mu, Sigma = Sigma, delta = 2.0)
x <- out$x

## Figure 4
par(pty = "s")
plot(x, xlim = c(-4,4), ylim = c(-4,4), xlab = "X1", ylab = "X2")
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

# Fig 6.a
z  <- influence.poa(fm, method = "normal")
plot(z)

# Fig 6.b
z  <- influence.poa(fm, method = "conformal")
plot(z)

# Fig 6.c
z  <- influence.poa(fm, method = "FI")
plot(z)

# Fig 6.d
z  <- influence.poa(fm, method = "SI")
plot(z)
