## loading dataset and reading R sources
load("../data/PSG.rda")
source("../code/ccc.R")
source("../code/poa.R")
source("../code/ccc.influence.R")
source("../code/poa.influence.R")

## Summary statistics at the beginning of Section 2.1
fm <- fit.ccc(~ manual + automated, data = PSG)
fm$xbar # print mean
#   manual automated
# 2.553891  2.308982
cov(fm$x) # print sample covariance matrix
#             manual automated
#manual    0.7711002 0.7027743
#automated 0.7027743 1.2522002
det(cov(fm$x)) # print det(S)
#[1] 0.47168

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
o <- boxplot(diff, plot = FALSE) # saving summaries
o <- c(min(diff), o$stats, max(diff))
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

## confidence interval for the probability of agreement (PoA)
## (fixme: asymptotic confidence interval must be contained in [0,1])
confint.poa(fm, cad = 2)
# Output: (1st row of Table 2)
#    parameter         SE    lower    upper
#PoA 0.9753448 0.03410157 0.908507 1.042183

confint.poa(f0, cad = 2)
# Output: (2nd row of Table 2)
#    parameter          SE     lower    upper
#PoA 0.9999223 0.000821773 0.9983117 1.001533

## Fig 2
nobs <- 82 # nrow(PSG)
rhoc <- rep(0, nobs)
for (i in 1:nobs)
  rhoc[i] <- fit.ccc(x = ~ manual + automated, data = PSG, subset = -i)$ccc
obs <- c(1,30,79)
cutoff <- fm$ccc
#pdf(width = 11, height = 5)
plot(rhoc, type = "b", ylim = c(.64,.76), ylab = "CCC estimate")
abline(h = cutoff, col = "red", lty = 2, lwd = 2)
text(obs, rhoc[obs], as.character(obs), pos = 3)

## Fig 3
poac <- rep(0, nobs)
for (i in 1:nobs) {
  f1 <- fit.ccc(x = ~manual + automated, data = PSG, subset = -i)
  poac[i] <- confint.poa(f1, level = 0.95, cad = 2)[1]
}
obs <- c(1,30,35,79)
cutoff <- confint.poa(fm)[1]
#pdf(width = 11, height = 5)
plot(poac, type = "b", ylim = c(.974,.995), ylab = "Estimated probability of agreement")
abline(h = cutoff, col = "red", lty = 2, lwd = 2)
text(obs, poac[obs], as.character(obs), pos = 3)

## Table 5
fm <- fit.ccc(~ manual + automated, data = PSG)
f0 <- fit.ccc(~ manual + automated, data = PSG, subset = -c(1,30,79))
sm <- c(fm$xbar, det(fm$S))
s0 <- c(f0$xbar, det(f0$S))
rel <- 100 * (s0 - sm) / sm
out <- cbind(sm, s0, rel)
rownames(out) <- c("manual","automated","det")
colnames(out) <- c("all","removed","change (%)")
out
# Output:
#                all   removed  change (%)
#manual    2.5538907 2.5260423  -1.0904269
#automated 2.3089822 2.3133702   0.1900416
#det       0.4602458 0.1556957 -66.1711855

## plotting influence measures..

## Fig 7.a
z <- influence.ccc(fm, method = "normal")
plot(z, idn = 4)

## Fig 7.b
z <- influence.ccc(fm, method = "normal")
hmax <- abs(z$hmax)
h2nd <- abs(z$vectors[,2])
obs  <- c(1,30,35,79)
cut1 <- mean(hmax) + 2 * sd(hmax)
cut2 <- mean(h2nd) + 2 * sd(h2nd)
par(pty = "s")
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

## Table 6
out  <- matrix(0, nrow = 8, ncol = 6)
cols <- c(1,3,5,6)
fm <- fit.ccc(~ manual + automated, data = PSG)
out[1,cols] <- c(fm$ccc, confint.poa(fm)[1], fm$logLik, fm$AIC)
f0 <- fit.ccc(~ manual + automated, data = PSG, subset = -1)
out[2,cols] <- c(f0$ccc, confint.poa(f0)[1], f0$logLik, f0$AIC)
f0 <- fit.ccc(~ manual + automated, data = PSG, subset = -30)
out[3,cols] <- c(f0$ccc, confint.poa(f0)[1], f0$logLik, f0$AIC)
f0 <- fit.ccc(~ manual + automated, data = PSG, subset = -79)
out[4,cols] <- c(f0$ccc, confint.poa(f0)[1], f0$logLik, f0$AIC)
f0 <- fit.ccc(~ manual + automated, data = PSG, subset = -c(1,30))
out[5,cols] <- c(f0$ccc, confint.poa(f0)[1], f0$logLik, f0$AIC)
f0 <- fit.ccc(~ manual + automated, data = PSG, subset = -c(1,79))
out[6,cols] <- c(f0$ccc, confint.poa(f0)[1], f0$logLik, f0$AIC)
f0 <- fit.ccc(~ manual + automated, data = PSG, subset = -c(30,79))
out[7,cols] <- c(f0$ccc, confint.poa(f0)[1], f0$logLik, f0$AIC)
f0 <- fit.ccc(~ manual + automated, data = PSG, subset = -c(1,30,79))
out[8,cols] <- c(f0$ccc, confint.poa(f0)[1], f0$logLik, f0$AIC)
rel <- 100 * (out[,1] - out[1,1]) / out[1,1]
out[,2] <- rel
rel <- 100 * (out[,3] - out[1,3]) / out[1,3]
out[,4] <- rel
rownames(out) <- c("none","1","30","79","1,30","1,79","30,79","1,30,79")
colnames(out) <- c("CCC","change (%)","PoA","change (%)","logLik","AIC")
out
# Output:
#              CCC change (%)       PoA change (%)    logLik      AIC
#none    0.6744407   0.000000 0.9753448   0.000000 -200.8901 411.7803
#1       0.7152318   6.048128 0.9860342   1.095963 -192.9896 395.9791
#30      0.7487737  11.021425 0.9922783   1.736157 -186.4483 382.8966
#79      0.7240161   7.350593 0.9895527   1.456703 -185.9074 381.8147
#1,30    0.7946708  17.826635 0.9974509   2.266486 -175.7411 361.4822
#1,79    0.7701287  14.187752 0.9960953   2.127501 -175.9815 361.9630
#30,79   0.8078412  19.779421 0.9988321   2.408099 -166.6881 343.3761
#1,30,79 0.8602722  27.553415 0.9999223   2.519882 -150.7281 311.4563

## Table 7
fm <- fit.ccc(~ manual + automated, data = PSG)
confint(fm, method = "asymp")
#    parameter         SE    lower     upper
#CCC 0.6744407 0.05626925 0.564155 0.7847264

confint(fm, method = "z-transform")
#  parameter        SE     lower     upper
#z 0.6744407 0.1028567 0.5487102 0.7703369

f0 <- fit.ccc(~ manual + automated, data = PSG, subset = -c(1,30,79))
confint(f0, method = "asymp")
#    parameter         SE     lower    upper
#CCC 0.8602722 0.02795145 0.8054883 0.915056

confint(f0, method = "z-transform")
#  parameter        SE     lower     upper
#z 0.8602722 0.1071212 0.7945408 0.9060751
