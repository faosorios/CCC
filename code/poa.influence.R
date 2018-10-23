## $ID: poa.influence.R, last updated 2018/09/18
## Author: Felipe Osorio

influence.poa <- function(object, cad = 2, method = "FI")
{ ## influence measures for the probability of agreement
  Call <- match.call()
  is.first <- FALSE

  switch(method,
         "FI" = { # first-order influence measure
           grad <- derivatives.poa(object, cad, which = 1)
           hmax <- grad / sqrt(sum(grad^2))
           is.first <- TRUE
         },
         "SI" = { # second-order influence measure
           derv <- derivatives.poa(object, which = 1:2)
           curv <- derv$hessian + diag(derv$gradient)
         },
         "normal" = { # normal curvature
           derv <- derivatives.poa(object, which = 1:2)
           grad <- derv$gradient
           n <- length(grad)
           scaling <- (1 + sum(grad^2)) * (diag(n) + outer(grad, grad))
           curv <- solve(scaling, derv$hessian)
         },
         "conformal" = { # normal conformal curvature
           derv <- derivatives.poa(object, which = 1:2)
           grad <- derv$gradient
           hess <- derv$hessian
           n <- length(grad)
           scaling <- sqrt(sum(diag(crossprod(hess)))) * (diag(n) + outer(grad, grad))
           curv <- solve(scaling, hess)
         },
         stop("method =", method, " is not implemented.")
        )

  ## compute largest eigenvectors and create the output object
  if (!is.first) {
    rs <- svd(curv, nv = 0)
    which <- abs(rs$d) < .Machine$double.eps
    hmax  <- rs$u[,1]
    out <- list(hmax = hmax, vectors = rs$u[,!which])
  } else
    out <- list(hmax = hmax, gradient = grad)

  out$call <- Call
  class(out) <- "influence.poa"
  out
}


print.influence.poa <- function(x, idn = 3, ...)
{ # print influential observations
  cat("Call:\n")
  dput(x$call, control = NULL)
  hmax <- abs(x$hmax)
  cutoff <- mean(hmax) + 2 * sd(hmax)
  n <- length(hmax)
  cat("\nCutoff:", format(cutoff))

  if (is.null(idn))
    idn <- 0
  else {
    idn <- as.integer(idn)
    if(idn < 0 || idn > n)
    stop(paste("`idn' must be in { 1,..,",n,"}"))
  }
  if(idn > 0) {
    idx   <- 1:idn
    show  <- order(hmax, decreasing = TRUE)[idx]
    show  <- intersect(show, (1:n)[hmax > cutoff])
    cat("\nInfluential observations:\n")
    print(sort(show))
  }
  invisible(x)
}

plot.influence.poa <- function(x, idn = 3, main = "")
{ # plotting direction of largest curvature
  hmax <- abs(x$hmax)
  cutoff <- mean(hmax) + 2 * sd(hmax)
  n <- length(hmax)

  if (is.null(idn))
    idn <- 0
  else {
    idn <- as.integer(idn)
    if(idn < 0 || idn > n)
    stop(paste("`idn' must be in { 1,..,",n,"}"))
  }
  if(idn > 0) {
    idx  <- 1:idn
    show <- order(-hmax)[idx]
    show <- intersect(show, (1:n)[hmax > cutoff])
  }

  plot(hmax, ylim = c(0,1), ylab = "hmax", main = main, font.main = 1)
  abline(h = cutoff, col = "red", lty = 2, lwd = 2)
  if (idn > 0) {
    which <- 1:n
    text(which[show], hmax[show], as.character(show), pos = 3)
  }
  invisible(x)
}

derivatives.poa <- function(object, cad = 2, which =  1:2)
{ ## compute 1st- and 2nd-order derivatives of probability of agreement (PoA)
  ## objective (influence) function
  o <- object
  if (!inherits(o, "ccc"))
    stop("use only with \"ccc\" objects")
  if (!is.numeric(which) || any(which < 1) || any(which > 2))
    stop("'which' must be in 1:2")
  show <- rep(FALSE, 2)
  show[which] <- TRUE
  both <- ifelse(all(show), TRUE, FALSE)

  ## centering variables
  z <- scale(o$x, center = o$xbar, scale = FALSE)

  ## extracting parameters and variables
  rhoc <- o$ccc
  mu1  <- o$xbar[1]
  mu2  <- o$xbar[2]
  sigma11 <- o$phi[1]
  sigma12 <- o$phi[2]
  sigma22 <- o$phi[3]
  z1 <- z[,1]
  z2 <- z[,2]

  ## preliminary definitions
  n <- o$dims[1]
  ones  <- rep(1, n)
  mu    <- mu1 - mu2
  sigma <- sqrt(sigma11 + sigma22 - 2. * sigma12)
  z     <- (cad - mu) / sigma
  diff  <- z1 - z2
  dn    <- dnorm(z)

  ## computation of 1st-order derivative of PoA
  mu.dot    <- diff / n
  sigma.dot <- ((n - 2) * diff * diff / n - sigma^2 * ones) / n
  unscaled  <- sigma * mu.dot + .5 * z * sigma.dot
  grad <- -2. * dn * unscaled / sigma^2

  if (show[1] && !both) # returning the gradient of PoA
    return (grad)

  ## computation of 2nd order derivatives of PoA
  mu.ddot    <- (outer(ones, diff) + outer(diff, ones)) / n
  sigma.ddot <- (2. * outer(ones, diff * diff) - outer(ones, sigma.dot) - outer(sigma.dot, ones)) / n^2
  term1 <- outer(sigma.dot, unscaled) / sigma^4
  term2 <- .5 * (outer(sigma.dot, mu.dot) - outer(mu.dot,sigma.dot)) / sigma
  term2 <- (term2 - sigma * mu.ddot + (cad - mu) * sigma.ddot) / sigma^2
  hess  <- 2. * dn * (term1 - term2 - z * outer(unscaled, unscaled) / sigma^4)
  hess  <- as.matrix(hess)

  if (show[2] && !both) # returning the hessian matrix of PoA
    return (hess)

  ## creating the output object
  out <- list(gradient = grad, hessian = hess)
  out
}

plot.influence.poa <- function(x, idn = 3, main = "")
{ # plotting direction of largest curvature
  hmax <- abs(x$hmax)
  cutoff <- mean(hmax) + 2 * sd(hmax)
  n <- length(hmax)

  if (is.null(idn))
    idn <- 0
  else {
    idn <- as.integer(idn)
    if(idn < 0 || idn > n)
    stop(paste("`idn' must be in { 1,..,",n,"}"))
  }
  if(idn > 0) {
    idx  <- 1:idn
    show <- order(-hmax)[idx]
    show <- intersect(show, (1:n)[hmax > cutoff])
  }

  plot(hmax, ylim = c(0,1), ylab = "hmax", main = main, font.main = 1)
  abline(h = cutoff, col = "red", lty = 2, lwd = 2)
  if (idn > 0) {
    which <- 1:n
    text(which[show], hmax[show], as.character(show), pos = 3)
  }
  invisible(x)
}
