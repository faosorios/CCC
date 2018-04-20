## loading dataset and reading R code (files must be at the working directory)
load("PSG.rda")
source("ccc.R")
source("ccc.influence.R")

## Fig 1
#pdf(width = 10)
par(mfrow = c(1,2))
o <- boxplot(PSG$manual, plot = FALSE) # saving summaries
o <- c(min(PSG$manual), o$stats, max(PSG$manual))
boxplot(PSG$manual, col = "gray", xlab = "Manual", ylim = c(0,6.2), boxwex = 0.5)
text(rep(1.3,7), round(o,3), round(o,3), cex = 0.7)
o <- boxplot(PSG$automated, plot = FALSE) # saving summaries
o <- c(min(PSG$automated), o$stats, max(PSG$automated))
boxplot(PSG$automated, col = "gray", xlab = "Automated", ylim = c(0,6.2), boxwex = 0.5)
text(rep(1.3,7), round(o,3), round(o,3), cex = 0.7)

## Fig 2.a
par(pty = "s")
plot(automated ~ manual, data = PSG, xlim = c(0,7), ylim = c(0,7))
abline(c(0,1), col = "red", lty = 2, lwd = 2)
obs <- c(1,30,79)
text(PSG$manual[obs], PSG$automated[obs], as.character(obs), pos = 3)

## Fig 2.b
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

# Output: (2nd row of Table 2)
#Call:
#fit.ccc(x = ~manual + automated, data = PSG, subset = -c(1, 30, 79))
#
#Coefficients:
#  estimate  variance  accuracy precision
#   0.8603    0.0008    0.9669    0.8897
#confint(f0, method = "asymp")
#    parameter         SE     lower    upper
#CCC 0.8602722 0.02795145 0.8054883 0.915056

## Fig 3
nobs <- 82 # nrow(PSG)
rhoc <- rep(0, nobs)
for (i in 1:nobs)
  rhoc[i] <- fit.ccc(x = ~manual + automated, data = PSG, subset = -i)$ccc
obs <- c(1,30,79)
#pdf(width = 11, height = 5)
plot(rhoc, type = "b", ylim = c(.64,.76), ylab = "CCC estimate")
abline(h = fm$ccc, col = "red", lty = 2, lwd = 2)
text(obs, rhoc[obs], as.character(obs), pos = 3)

## plotting influence measures..

## Fig 6.a
z <- influence.ccc(fm, method = "normal")
plot(z, idn = 4)

## Fig 6.b
hmax <- z$hmax
h2nd <- z$vectors[,2]
obs  <- c(1,30,35,79)
cut1 <- mean(hmax) + 2 * sd(hmax)
cut2 <- mean(h2nd) + 2 * sd(h2nd)
plot(h2nd, hmax, xlim = c(0,1), ylim = c(0,1), xlab = "h2nd", ylab = "hmax")
abline(h = cut1, col = "red", lty = 2, lwd = 2)
abline(v = cut2, col = "red", lty = 2, lwd = 2)
text(h2nd[obs], hmax[obs], as.character(obs), pos = 3)

## Fig 6.c
z <- influence.ccc(fm, method = "FI")
plot(z)

## Fig 6.d
z <- influence.ccc(fm, method = "SI")
plot(z)
