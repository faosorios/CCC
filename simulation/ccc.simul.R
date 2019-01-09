## Id: ccc.simul.R, last updated 2019/01/07
## Author: originally coded by Carla Leal, with contributions of Felipe Osorio

simul.ccc <- function(Nsize = 500, nobs = 100, mean, Sigma, delta = 3)
{ ## function to perform the simulation experiment (Section 4 of the manuscript)
  ok <- matrix(FALSE, nrow = Nsize, ncol = 5) # results container
  cutoff <- rep(0, 4)
  now <- proc.time()

  # Monte Carlo iterations
  for (i in 1:Nsize) {
    x <- rmnorm(n = nobs, mean = mean, Sigma = Sigma)
    # introducing an outlier to 2nd instrument
    x[1,2] <- x[1,2] + delta
    # fitting postulated model and computing influence measures
    fm <- fit.ccc(x)
    # normal curvature
    z <- influence.ccc(fm, method = "normal")
    hmax <- abs(z$hmax)
    cutoff[1] <- mean(hmax) + 2 * sd(hmax)
    ok[i,1] <- hmax[1] > cutoff[1]
    # conformal curvature
    z <- influence.ccc(fm, method = "conformal")
    hmax <- abs(z$hmax)
    cutoff[2] <- mean(hmax) + 2 * sd(hmax)
    ok[i,2] <- hmax[1] > cutoff[2]
    # FI measure
    z <- influence.ccc(fm, method = "FI")
    hmax1 <- abs(z$hmax)
    cutoff[3] <- mean(hmax1) + 2 * sd(hmax1)
    ok[i,3] <- hmax1[1] > cutoff[3]
    # SI measure
    z <- influence.ccc(fm, method = "conformal")
    hmax2 <- abs(z$hmax)
    cutoff[4] <- mean(hmax2) + 2 * sd(hmax2)
    ok[i,4] <- hmax2[1] > cutoff[4]
    # FI & SI measures
    ok[i,5] <- (hmax1[1] > cutoff[3]) && (hmax2[1] > cutoff[4])
  }

  mnames <- c("C","B","FI","SI","both")
  percentage <- apply(ok, 2, sum) / Nsize
  names(percentage) <- mnames
  names(cutoff) <- mnames[1:4]
  speed <- proc.time() - now
  out <- list(percentage = 100 * percentage, cutoffs = cutoff, speed = speed, x = x, hmax = hmax)
  out
}

summary.ccc <- function(Nsize = 500, nobs = 100, mean, Sigma) {
  out <- matrix(0, nrow = 5, ncol = 7)

  set.seed(5)
  out[,1] <- simul.ccc(Nsize, nobs, mean, Sigma, delta = 0.5)$percentage
  set.seed(23)
  out[,2] <- simul.ccc(Nsize, nobs, mean, Sigma, delta = 1.0)$percentage
  set.seed(53)
  out[,3] <- simul.ccc(Nsize, nobs, mean, Sigma, delta = 1.5)$percentage
  set.seed(167)
  out[,4] <- simul.ccc(Nsize, nobs, mean, Sigma, delta = 2.0)$percentage
  set.seed(239)
  out[,5] <- simul.ccc(Nsize, nobs, mean, Sigma, delta = 2.5)$percentage
  set.seed(347)
  out[,6] <- simul.ccc(Nsize, nobs, mean, Sigma, delta = 3.0)$percentage
  set.seed(433)
  out[,7] <- simul.ccc(Nsize, nobs, mean, Sigma, delta = 3.0)$percentage

  rownames(out) <- c("C","B","FI","SI","both")
  colnames(out) <- c("0.5","1.0","1.5","2.0","2.5","3.0","3.5")
  out
}
