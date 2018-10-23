## loading dataset and reading R sources
load("../data/PSG.rda")
source("../code/ccc.R")
source("../code/poa.R")
source("../code/ccc.influence.R")
source("../code/poa.influence.R")

## Fig 1.a
par(pty = "s")
plot(automated ~ manual, data = PSG, xlim = c(0,7), ylim = c(0,7))
abline(c(0,1), col = "red", lty = 2, lwd = 2)
obs <- c(1,30,79)
text(PSG$manual[obs], PSG$automated[obs], as.character(obs), pos = 3)

## Fig 1.b
average <- apply(PSG[,1:2], 1, mean)
diff <- PSG$manual - PSG$automated
qz <- qnorm(c(0.025, 0.975))
limits <- mean(diff) + sd(diff) %o% qz
par(pty = "s")
plot(average, diff, xlab = "average between methods", ylab = "difference between methods", ylim = c(-4,4))
abline(h = mean(diff), col = "red", lwd = 2)
abline(h = c(limits), col = "red", lty = 2, lwd = 2)
obs <- c(1,30,35,79)
text(average[obs], diff[obs], as.character(obs), pos = 3)

## Fig 1.c
par(pty = "s")
diff <- PSG$manual - PSG$automated
o <- boxplot(PSG$manual, plot = FALSE) # saving summaries
o <- c(min(PSG$manual), o$stats, max(PSG$manual))
o <- o[-3] # removing quantile 25%
boxplot(diff, col = "gray", xlab = "Difference between methods", ylim = c(-4,4), boxwex = 0.5)
text(rep(1.3,6), round(o,3), round(o,3), cex = 0.7)


## compute CCC estimate and confidence interval
fm <- fit.ccc(~ manual + automated, data = PSG)
fm # print results
confint(fm, method = "asymp") # print 95% CI

# Output: (1st row of Table 1)
#Call:
#fit.ccc(x = ~manual + automated, data = PSG)
#
#Coefficients:
#  estimate  variance  accuracy precision
#   0.6744    0.0032    0.9430    0.7152
#confint(fm, method = "asymp")
#    parameter         SE    lower     upper
#CCC 0.6744407 0.05626925 0.564155 0.7847264

## removing obs 1,30 and 79
f0 <- fit.ccc(~ manual + automated, data = PSG, subset = -c(1,30,79))
f0 # print results
confint(f0, method = "asymp") # print 95% CI

# Output: (2nd row of Table 1)
#Call:
#fit.ccc(x = ~manual + automated, data = PSG, subset = -c(1, 30, 79))
#
#Coefficients:
#  estimate  variance  accuracy precision
#   0.8603    0.0008    0.9669    0.8897
#confint(f0, method = "asymp")
#    parameter         SE     lower    upper
#CCC 0.8602722 0.02795145 0.8054883 0.915056

## Fig 2
nobs <- 82 # nrow(PSG)
rhoc <- rep(0, nobs)
for (i in 1:nobs)
  rhoc[i] <- fit.ccc(x = ~manual + automated, data = PSG, subset = -i)$ccc
obs <- c(1,30,79)
#pdf(width = 11, height = 5)
plot(rhoc, type = "b", ylim = c(.64,.76), ylab = "CCC estimate")
abline(h = fm$ccc, col = "red", lty = 2, lwd = 2)
text(obs, rhoc[obs], as.character(obs), pos = 3)

## Fig 3
poac <- rep(0, nobs)
for (i in 1:nobs) {
  f1 <- fit.ccc(x = ~manual + automated, data = PSG, subset = -i)$ccc
  poac[i] <- confint.poa(f1, level = 0.95, cad = 2)[1]
}
obs <- c(1,30,35,79)
#pdf(width = 11, height = 5)
plot(rhoc, type = "b", ylab = "Estimated probability of agreement")
abline(h = confint.poa(fm)[1], col = "red", lty = 2, lwd = 2)
text(obs, poac[obs], as.character(obs), pos = 3)

## plotting influence measures..

## Fig 7.a
z <- influence.ccc(fm, method = "normal")
plot(z, idn = 4)

## Fig 7.b
hmax <- z$hmax
h2nd <- z$vectors[,2]
obs  <- c(1,30,35,79)
cut1 <- mean(hmax) + 2 * sd(hmax)
cut2 <- mean(h2nd) + 2 * sd(h2nd)
plot(h2nd, hmax, xlim = c(0,1), ylim = c(0,1), xlab = "h2nd", ylab = "hmax")
abline(h = cut1, col = "red", lty = 2, lwd = 2)
abline(v = cut2, col = "red", lty = 2, lwd = 2)
text(h2nd[obs], hmax[obs], as.character(obs), pos = 3)

## Fig 7.c
z <- influence.ccc(fm, method = "FI")
plot(z)

## Fig 7.d
z <- influence.ccc(fm, method = "SI")
plot(z)

## Fig 8.a
z <- influence.poa(fm, method = "normal")
plot(z)

## Fig 8.b
z <- influence.poa(fm, method = "conformal")
plot(z)

## Fig 8.c
z <- influence.poa(fm, method = "FI")
plot(z)

## Fig 8.d
z <- influence.poa(fm, method = "SI")
plot(z)
