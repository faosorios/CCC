## Id: plot.ccc.R
## Author: Felipe Osorio
## Last update: 2018-04-20

load("x.rda")

delta <- seq(0, 3.5, length = 200)
rhoc  <- double(200)
rhoi  <- double(200)
x.orig <- x[1,2]

for(i in seq_along(delta)) {
  x[1,2] <- x[1,2] + delta[i]
  rhoc[i] <- fit.ccc(x)$ccc
  rhoi[i] <- fit.ccc(x, subset = -1)$ccc
  x[1,2] <- x.orig
}

# Figure
par(pty = "s")
plot(delta, rhoc, type = "l", lwd = 2, lty = 2, xlab = "delta", ylab = "CCC")
lines(delta, rhoi, col = "red", lwd = 2)
