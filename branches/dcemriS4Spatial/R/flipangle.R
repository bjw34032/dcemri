##
##
## Copyright (c) 2009-2011 Brandon Whitcher and Volker Schmid
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
## $Id: flipangle.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#############################################################################
## dam() = double-angle method
#############################################################################

dam <- function(low, high, low.deg) {
  alpha <- acos(abs(high /(2*low)))
  (180/pi * alpha) / low.deg # radians to degrees
}

#############################################################################
## R10.lm() = estimate R1 using Levenburg-Marquardt
#############################################################################

R10.lm <- function(signal, alpha, TR, guess, control=nls.lm.control()) {
  func <- function(x, y) {
    R1 <- x[1]
    m0 <- x[2]
    signal <- y[[1]]
    theta <- pi/180 * y[[2]]    # degrees to radians
    TR <- y[[3]]
    signal -
      m0 * sin(theta) * (1 - exp(-TR*R1)) / (1 - cos(theta) * exp(-TR*R1))
  }
  require("minpack.lm") # Levenberg-Marquart fitting
  out <- nls.lm(par=guess, fn=func, control=control,
                y=list(signal, alpha, TR))
  list(R1=out$par[1], S0=out$par[2], info=out$info, message=out$message)
}

#############################################################################
## E10.lm() = estimate exp(-TR*R1) using Levenburg-Marquardt
#############################################################################

E10.lm <- function(signal, alpha, guess, control=nls.lm.control()) {
  func <- function(x, signal, alpha) {
    E1 <- x[1]
    m0 <- x[2]
    theta <- pi/180 * alpha    # degrees to radians
    signal - m0 * sin(theta) * (1 - E1) / (1 - cos(theta) * E1)
  }
  require("minpack.lm") # Levenberg-Marquart fitting
  out <- nls.lm(par=guess, fn=func, control=control, signal=signal,
                alpha=alpha)
  list(E10=out$par[1], m0=out$par[2], hessian=out$hessian, info=out$info,
       message=out$message)
}

#############################################################################
## setGeneric("R1.fast")
#############################################################################

setGeneric("R1.fast", function(flip, ...) standardGeneric("R1.fast"))
setMethod("R1.fast", signature(flip="array"),
          function(flip, flip.mask, fangles, TR, control=nls.lm.control(),
                   multicore=FALSE, verbose=FALSE) 
	    .dcemriWrapper("R1.fast", flip, flip.mask, fangles, TR, control,
                           multicore, verbose))

#############################################################################
## R1.fast()
#############################################################################

.R1.fast <- function(flip, flip.mask, fangles, TR, control=nls.lm.control(),
                     multicore=FALSE, verbose=FALSE) {

  if (length(dim(flip)) != 4) { # Check flip is a 4D array
    stop("Flip-angle data must be a 4D array.")
  }
  if (!is.logical(flip.mask)) { # Check flip.mask is logical
    stop("Mask must be logical.")
  }
    
  X <- nrow(flip)
  Y <- ncol(flip)
  Z <- nsli(flip)
  nvoxels <- sum(flip.mask)
  
  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  flip.mat <- matrix(flip[flip.mask], nrow=nvoxels)
  if (is.array(fangles)) {
    fangles.mat <- matrix(fangles[flip.mask], nrow=nvoxels)
  } else {
    fangles.mat <- matrix(fangles, nrow=nvoxels, ncol=length(fangles),
                          byrow=TRUE)
  }
  flip.list <- vector("list", nvoxels)
  for (k in 1:nvoxels) {
    flip.list[[k]] <- list(signal=flip.mat[k,], angles=fangles.mat[k,])
  }
  if (verbose) {
    cat("  Calculating R10 and M0...", fill=TRUE)
  }
  if (multicore && require("multicore")) {
    fit.list <- mclapply(flip.list, function(x) {
      E10.lm(x$signal, x$angles, guess=c(1, mean(x$signal)), control)
    })
  } else {
    fit.list <- lapply(flip.list, function(x) {
      E10.lm(x$signal, x$angles, guess=c(1, mean(x$signal)), control)
    })
  }
  rm(flip.list) ; gc()
  R10 <- M0 <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  for (k in 1:nvoxels) {
    if (fit.list[[k]]$info > 0 && fit.list[[k]]$info < 5) {
      R10$par[k] <- log(fit.list[[k]]$E10) / -TR
      M0$par[k] <- fit.list[[k]]$m0
      R10$error[k] <- sqrt(fit.list[[k]]$hessian[1,1])
      M0$error[k] <- sqrt(fit.list[[k]]$hessian[2,2])
    } else {
      R10$par[k] <- M0$par[k] <- R10$error[k] <- M0$error[k] <- NA
    }
  }

  ##R10 <- M0 <- numeric(nvoxels)
  ##for (k in 1:nvoxels) {
  ##  fit <- E10.lm(flip.mat[k,], fangles.mat[k,],
  ##                guess=c(1, mean(flip.mat[k,])), control)
  ##  if (fit$info == 1 || fit$info == 2 || fit$info == 3) {
  ##    R10[k] <- log(fit$E10) / -TR
  ##    M0[k] <- fit$m0
  ##  } else {
  ##    R10[k] <- M0[k] <- NA
  ##  }
  ##}

  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  R10.array <- M0.array <- R10error <- M0error <- array(NA, c(X,Y,Z))
  R10.array[flip.mask] <- R10$par
  M0.array[flip.mask] <- M0$par
  R10error[flip.mask] <- R10$error
  M0error[flip.mask] <- M0$error

  list(M0 = M0.array, R10 = R10.array, M0.error = NULL, R10.error = NULL)
}

#############################################################################
## setGeneric("CA.fast")
#############################################################################

setGeneric("CA.fast", function(dynamic, ...) standardGeneric("CA.fast"))
setMethod("CA.fast", signature(dynamic="array"),
	  function(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4,
                   control=nls.lm.control(maxiter=200), multicore=FALSE,
                   verbose=FALSE) 
	    .dcemriWrapper("CA.fast", dynamic, dyn.mask, dangle, flip,
                           fangles, TR, r1, control, multicore, verbose))

#############################################################################
## CA.fast() = estimate contrast-agent concentration and other stuff
#############################################################################

.CA.fast <- function(dynamic, dyn.mask, dangle, flip, fangles, TR,
                     r1=4, control=nls.lm.control(maxiter=200),
                     multicore=FALSE, verbose=FALSE) {

  if (length(dim(flip)) != 4) { # Check flip is a 4D array
    stop("Flip-angle data must be a 4D array.")
  }
  if (!is.logical(dyn.mask)) { # Check dyn.mask is logical
    stop("Mask must be logical.")
  }
  if (! identical(length(fangles), ntim(flip)) &&
      ! isTRUE(all.equal(dim(flip), dim(fangles)))) {
    ## Check that #(flip angles) are equal
    stop("Number of flip angles must agree with dimension of flip-angle data.")
  }
  
  R1est <- R1.fast(flip, dyn.mask, fangles, TR, control, multicore, verbose)
  
  if (verbose) {
    cat("  Calculating concentration...", fill=TRUE)
  }
  theta <- dangle * pi/180
  A <- sweep(sweep(dynamic, 1:3, dynamic[,,,1], "-"),
             1:3, R1est$M0, "/") / sin(theta)
  B <- (1 - exp(-TR * R1est$R10)) / (1 - cos(theta) * exp(-TR * R1est$R10))
  AB <- sweep(A, 1:3, B, "+")
  rm(A,B)
  R1t <- -(1/TR) * log((1 - AB) / (1 - cos(theta) * AB))
  rm(AB)
  conc <- sweep(R1t, 1:3, R1est$R10, "-") / r1

  list(M0 = R1est$M0, R10 = R1est$R10, R1t = R1t, conc = conc)
}

#############################################################################
## setGeneric("CA.fast2")
#############################################################################

setGeneric("CA.fast2", function(dynamic, ...) standardGeneric("CA.fast2"))
setMethod("CA.fast2", signature(dynamic="array"),
	  function(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4,
                   verbose=FALSE) 
          .dcemriWrapper("CA.fast2", dynamic, dyn.mask, dangle, flip,
                         fangles, TR, r1, verbose))

#############################################################################
## CA.fast2()
#############################################################################

.CA.fast2 <- function(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4,
                     verbose=FALSE) {
  
  if (length(dim(flip)) != 4) {  # Check flip is a 4D array
    stop("Flip-angle data must be a 4D array.")
  }
  if (! identical(length(fangles), ntim(flip)) &&
      ! isTRUE(all.equal(dim(flip), dim(fangles)))) {
    ## Check that #(flip angles) are equal
    stop("Number of flip angles must agree with dimension of flip-angle data.")
  }
  ##if (ntim(flip) != 2 || length(fangles) != 2) {
  ##  stop("Only two flip angles are allowed.")
  ##}
  if (!is.logical(dyn.mask)) { # Check dyn.mask is logical
    stop("Mask must be logical.")
  }
  nangles <- length(fangles)
  nvoxels <- sum(dyn.mask)
  M <- nrow(flip)
  N <- ncol(flip)
  Z <- nsli(dynamic)
  W <- ntim(dynamic)
  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  if (is.array(fangles)) {
    fangles.mat <- matrix(fangles[dyn.mask], nrow=nvoxels)
  } else {
    fangles.mat <- matrix(fangles, nrow=nvoxels, ncol=length(fangles),
                          byrow=TRUE)
  }
  dyn.mat <- matrix(dynamic[dyn.mask], nvoxels)
  flip.mat <- matrix(flip[dyn.mask], nvoxels)
  R10 <- M0 <- numeric(nvoxels)
  if (verbose) {
    cat("  Calculating R10 and M0...", fill=TRUE)
  }
  x <- flip.mat / tan(pi * fangles.mat / 180)
  x <- ifelse(is.finite(x), x, 0)
  y <- flip.mat / sin(pi * fangles.mat / 180)
  y <- ifelse(is.finite(y), y, 0)
  for (k in 1:nvoxels) {
    #x <- c(flip.mat[k,1] / tan(pi * fangles.mat[k,1] / 180),
    #       flip.mat[k,2] / tan(pi * fangles.mat[k,2] / 180))
    #x <- ifelse(is.finite(x), x, 0)
    #y <- c(flip.mat[k,1] / sin(pi * fangles.mat[k,1] / 180),
    #       flip.mat[k,2] / sin(pi * fangles.mat[k,2] / 180))
    #y <- ifelse(is.finite(y), y, 0)
    fit <- lsfit(x[k, ], y[k, ])$coefficients
    R10[k] <- log(fit[2]) / -TR
    M0[k] <- fit[1] / (1 - fit[2])
  }
  if (verbose) {
    cat("  Calculating concentration...", fill=TRUE)
  }
  theta <- dangle * pi/180
  CD <- conc <- matrix(NA, nvoxels, W)
  B <- (1 - exp(-TR * R10)) / (1 - cos(theta) * exp(-TR * R10))
  A <- (dyn.mat - dyn.mat[,1]) / M0 / sin(theta)
  R1t <- -(1/TR) * log((1 - (A+B)) / (1 - cos(theta) * (A+B)))
  conc <- (R1t - R10) / r1
  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  R10.array <- M0.array <- array(NA, c(M,N,Z))
  R10.array[dyn.mask] <- R10
  M0.array[dyn.mask] <- M0
  conc.array <- R1t.array <- array(NA, c(M,N,Z,W))
  mask4D <- array(dyn.mask, c(M,N,Z,W))
  conc.array[mask4D] <- unlist(conc)
  R1t.array[mask4D] <- unlist(R1t)
  list(M0 = M0.array, R10 = R10.array, R1t = R1t.array, conc = conc.array)
}
