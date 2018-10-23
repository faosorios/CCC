## $ID: poa.R, last updated 2018/04/20
## Author: Felipe Osorio

confint.poa <- function(object, level = 0.95, cad = 2)
{
  ## extracting elements from object
  n <- nrow(object$x)
  mu1 <- object$xbar[1]
  mu2 <- object$xbar[2]
  sigma11 <- object$phi[1]
  sigma12 <- object$phi[2]
  sigma22 <- object$phi[3]

  ## computing PoA
  mu <- mu1 - mu2
  sigma2 <- sigma11 + sigma22 - 2. * sigma12
  z    <- (cad - mu) / sqrt(sigma2)
  dn   <- dnorm(z)
  prob <- 2. * pnorm(z) - 1. # pnorm(z) - pnorm(-z)

  ## asymptotic normal approximation
  var.psi <- sigma11^2 + sigma22^2 + 2. * sigma12^2 + 8. * (sigma11 - sigma12) * (sigma22 - sigma12)
  var.psi <- 2. * (sigma2 + z^2 * var.psi / sigma2)
  var.psi <- (dn^2 / sigma2) * var.psi / n

  ## computing CI
  a  <- (1 - level) / 2
  a  <- c(a, 1 - a)
  qz <- qnorm(a)
  SE <- sqrt(var.psi)
  ci <- prob + SE %o% qz
  ci <- cbind(prob, SE, ci)

  dimnames(ci) <- list("PoA", c("parameter", "SE", "lower", "upper"))
  ci
}
