## $ID: ccc.R, last updated 2018/04/20
## Author: Felipe Osorio

fit.ccc <-
function(x, data, subset, na.action)
{ # estimation of the concordance correlation coefficient
  Call <- match.call()
  if (missing(x))
    stop("'x' is not supplied")
  if (inherits(x, "formula")) {
    mt <- terms(x, data = data)
    if (attr(mt, "response") > 0)
      stop("response not allowed in formula")
    attr(mt, "intercept") <- 0
    mf <- match.call(expand.dots = FALSE)
    names(mf)[names(mf) == "x"] <- "formula"
    mf$qpar <- mf$maxiter <- mf$tol <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    na.act <- attr(mf, "na.action")
    z <- model.matrix(mt, mf)
  }
  else {
    z <- as.matrix(x)
    if (!missing(subset))
      z <- z[subset, , drop = FALSE]
    if (!missing(na.action))
      z <- na.omit(z)
    else
      z <- na.fail(z)
  }
  if (!is.numeric(z))
    stop("fit.ccc applies only to numerical variables")
  znames <- dimnames(z)[[2]]
  dz <- dim(z)
  n <- dz[1]
  p <- dz[2]
  if (p > 2)
    stop("fit.ccc is not implemented for p > 2")

  ## estimating mean vector and covariance matrix
  wts <- rep(1, n)
  o   <- cov.wt(z, wt = wts, center = TRUE, method = "ML")
  xbar <- o$center
  xcov <- o$cov
  mu   <- as.vector(xbar)
  phi  <- xcov[lower.tri(xcov, diag = TRUE)]

  ## computing CCC
  diff  <- mu[1] - mu[2]
  ratio <- phi[1] / phi[3]
  a <- diff / ((phi[1] * phi[3])^.25) # location shift
  b <- sqrt(ratio) # scale shift

  rhoc <- 2. * phi[2] / (phi[1] + phi[3] + diff^2)
  accu <- 2. / (b + 1 / b + a^2)
  r <- rhoc / accu

  ## asymptotic normal approximation (correction by Lin, Biometrics 56, pp.325, 2000)
  var.rhoc <- (1 - r^2) * rhoc^2 * (1 - rhoc^2) / r^2 + 2 * rhoc^3 * (1 - rhoc) * a^2 / r - 0.5 * (rhoc * a)^4 / r^2
  var.rhoc <- var.rhoc / (n - 2)
  ## z-transformation
  var.z <- var.rhoc / (1 - rhoc^2)^2

  ## creating the output object
  out <- list(call = Call, x = z, dims = dz, xbar = xbar, S = xcov, phi = phi,
              ccc = rhoc, var.ccc = var.rhoc, accuracy = accu, precision = r,
              z = atanh(rhoc), var.z = var.z)
  class(out) <- "ccc"
  out
}

print.ccc <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  dput(x$call, control = NULL)
  cf <- c(x$ccc, x$var.ccc, x$accuracy, x$precision)
  names(cf) <- c("estimate","variance","accuracy","precision")
  cat("\nCoefficients:\n ")
  print(format(round(cf, digits = digits)), quote = F, ...)
  invisible(x)
}

confint.ccc <- function(object, method = "z-transform", level = 0.95)
{ # confidence intervals for the concordance correlation coefficient
  z2rho <- function(x) (exp(2 * x) - 1) / (exp(2 * x) + 1)
  a  <- (1 - level) / 2
  a  <- c(a, 1 - a)
  qz <- qnorm(a)
  switch(method,
         "asymp" = { # asymptotic interval
           cname <- "CCC"
           cf <- object$ccc
           SE <- sqrt(object$var.ccc)
           ci <- cf + SE %o% qz
         },
         "z-transform" = { # z-transform (default)
           cname <- "z"
           cf <- object$z
           SE <- sqrt(object$var.z)
           ci <- cf + SE %o% qz
           cf <- z2rho(cf)
           SE <- z2rho(SE)
           ci <- z2rho(ci)
         },
         stop("method = ", method, " is not implemented."))
  ci <- cbind(cf, SE, ci)
  dimnames(ci) <- list(cname, c("parameter", "SE", "lower", "upper"))
  ci
}
