## Id: ccc.influence.R
## Author: originally coded by Carla Leal, with contributions of Felipe Osorio
## Last update: 2018-04-12

influence.ccc <- function(object, method = "FI")
{ ## influence measures for the concordance correlation coefficient
  switch(method,
         "FI" = { # first-order influence measure
           grad <- derivatives.ccc(object, which = 1)
           curv <- outer(grad, grad)
         },
         "SI" = { # second-order influence measure
           derv <- derivatives.ccc(object, which = 1:2)
           curv <- derv$hessian + diag(derv$gradient)
         },
         "CF" = { # maximum slope influence measure
           grad <- derivatives.ccc(object, which = 1)
           hmax <- grad / sqrt(sum(grad^2))
         },
         "normal" = { # normal curvature
           derv <- derivatives.ccc(object, which = 1:2)
           grad <- derv$gradient
           n <- length(grad)
           scaling <- (1 + sum(grad^2)) * (diag(n) + outer(grad, grad))
           curv <- solve(scaling, derv$hessian)
         },
         "conformal" = { # normal conformal curvature
           derv <- derivatives.ccc(object, which = 1:2)
           grad <- derv$gradient
           hess <- derv$hessian
           n <- length(grad)
           scaling <- sqrt(sum(diag(crossprod(hess)))) * (diag(n) + outer(grad, grad))
           curv <- solve(scaling, derv$hessian)
         },
         stop("method =", method, " is not implemented.")
        )

  ## compute largest eigenvectors and create the output object
  if (method != "CF") {
    rs <- svd(curv, nv = 0)
    which <- abs(rs$d) < .Machine$double.eps
    hmax  <- rs$u[,1]
    total <- diag(curv)
    out <- list(hmax = hmax, total = total, vectors = rs$u[,!which])
  } else
    out <- list(hmax = hmax, gradient = grad)

  class(out) <- "influence.ccc"
  out
}

plot.influence.ccc <- function(x, idn = 3, main = "")
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

derivatives.ccc <- function(object, which =  1:2)
{ ## compute 1st- and 2nd-order derivatives of CCC objective function
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
  ones <- rep(1, n)
  diff <- mu1 - mu2
  s11  <- z1 * z1 - sigma11 * ones
  s12  <- z1 * z2 - sigma12 * ones
  s22  <- z2 * z2 - sigma22 * ones
  rel  <- rhoc / (n * sigma12)

  ## computation of 1st-order derivative of CCC
  star <- s11 + s22 + 2 * diff * (z1 - z2)
  grad <- rel * as.vector(s12 - .5 * rhoc * star)

  if (show[1] && !both) # returning the gradient of CCC
    return (grad)

  ## computation of 2nd order derivatives of CCC
  rel  <- rhoc^2 / sigma12
  z11  <- sigma11^2 * outer(ones, ones) - outer(ones, z1 * z1) - outer(z1, z1)
  z12  <- sigma12^2 * outer(ones, ones) - outer(ones, z1 * z2) - outer(z1, z2)
  z22  <- sigma22^2 * outer(ones, ones) - outer(ones, z2 * z2) - outer(z2, z2)
  term1 <- outer(grad / sigma12 - rel^2 * s12 / n, (s11 + s22) / n) + rel * (z11 + z22) / n^2
  term2 <- rel * (outer(z1 - z2, z1 - z2) - (n + 1) * diff * outer(z1 - z2, ones)) / n^2
  term3 <- 2 * rel * z12 / n^2 - .5 * rel * outer(s11 + s22 + 2 * diff * (z1 - z2), s12) / (sigma12 * n^2)
  hess  <- as.matrix(term3 - (term1 + term2))

  if (show[2] && !both) # returning the hessian matrix of CCC
    return (hess)

  ## creating the output object
  out <- list(gradient = grad, hessian = hess)
  out
}
