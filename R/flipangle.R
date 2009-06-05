R10.lm <- function(signal, alpha, TR, guess, nprint=0) {
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
  out <- nls.lm(par=guess, fn=func, control=list(nprint=nprint, maxiter=150),
               y=list(signal, alpha, TR))
  list(R1=out$par[1], S0=out$par[2], info=out$info, message=out$message)
}

E10.lm <- function(signal, alpha, guess, nprint=0) {
  func <- function(x, signal, alpha) {
    E1 <- x[1]
    m0 <- x[2]
    theta <- pi/180 * alpha    # degrees to radians
    signal -
      m0 * sin(theta) * (1 - E1) / (1 - cos(theta) * E1)
  }
  require("minpack.lm") # Levenberg-Marquart fitting
  out <- nls.lm(par=guess, fn=func, control=list(nprint=nprint, maxiter=150),
               signal=signal, alpha=alpha)
  list(E10=out$par[1], m0=out$par[2], info=out$info, message=out$message)
}

CA.fast <- function(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4.39,
                    verbose=FALSE) {

  if (length(dim(flip)) != 4)  # Check flip is a 4D array
    stop("Flip-angle data must be a 4D array.")
  if (!is.logical(dyn.mask))  # Check dyn.mask is logical
    stop("Mask must be logical.")
  
  nangles <- length(fangles)
  nvoxels <- sum(dyn.mask)
  M <- nrow(flip)
  N <- ncol(flip)
  Z <- nsli(dynamic)
  W <- ntim(dynamic)
  
  if (verbose)
    cat("  Deconstructing data...", fill=TRUE)
  dyn.mat <- matrix(dynamic[dyn.mask], nvoxels)
  flip.mat <- matrix(flip[dyn.mask], nvoxels)
  R10 <- M0 <- numeric(nvoxels)

  if (verbose)
    cat("  Calculating R10 and M0...", fill=TRUE)
  for(k in 1:nvoxels) {
    fit <- E10.lm(flip.mat[k,], fangles, guess=c(1, mean(flip.mat[k,])))
    if(fit$info == 1 || fit$info == 2 || fit$info == 3) {
      R10[k] <- log(fit$E10) / -TR
      M0[k] <- fit$m0
    } else {
      R10[k] <- M0[k] <- NA
    }
  }

  if (verbose)
    cat("  Calculating concentration...", fill=TRUE)
  theta <- dangle * pi/180
  CD <- conc <- matrix(NA, nvoxels, W)
  B <- (1 - exp(-TR * R10)) / (1 - cos(theta) * exp(-TR * R10))
  A <- (dyn.mat - dyn.mat[,1]) / M0 / sin(theta)
  R1t <- -(1/TR) * log((1 - (A+B)) / (1 - cos(theta) * (A+B)))
  conc <- (R1t - R10) / r1
  rm(A,B,CD)

  if (verbose)
    cat("  Reconstructing results...", fill=TRUE)
  R10.array <- M0.array <- array(NA, c(M,N,Z,1))
  R10.array[dyn.mask] <- R10
  M0.array[dyn.mask] <- M0
  conc.array <- R1t.array <- array(NA, c(M,N,Z,W))
  mask4D <- array(dyn.mask, c(M,N,Z,W))
  conc.array[mask4D] <- unlist(conc)
  R1t.array[mask4D] <- unlist(R1t)
  
  list(M0 = M0.array, R10 = R10.array, R1t = R1t.array, conc = conc.array)
}

CA.fast2 <- function(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4.39,
                     verbose=FALSE) {
  
  if (length(dim(flip)) != 4)  # Check flip is a 4D array
    stop("Flip-angle data must be a 4D array.")
  if (!is.logical(dyn.mask))  # Check dyn.mask is logical
    stop("Mask must be logical.")
  
  nangles <- length(fangles)
  nvoxels <- sum(dyn.mask)
  M <- nrow(flip)
  N <- ncol(flip)
  if(ntim(flip) != 2 || length(fangles) != 2)
    stop("Only two flip angles are allowed.")
  Z <- nsli(dynamic)
  W <- ntim(dynamic)
  
  if (verbose)
    cat("  Deconstructing data...", fill=TRUE)
  dyn.mat <- matrix(dynamic[dyn.mask], nvoxels)
  flip.mat <- matrix(flip[dyn.mask], nvoxels)
  R10 <- M0 <- numeric(nvoxels)

  if (verbose)
    cat("  Calculating R10 and M0...", fill=TRUE)
  for(k in 1:nvoxels) {
    x <- c(flip.mat[k,1] / tan(pi*fangles[1]/180),
           flip.mat[k,2] / tan(pi*fangles[2]/180))
    y <- c(flip.mat[k,1] / sin(pi*fangles[1]/180),
           flip.mat[k,2] / sin(pi*fangles[2]/180))
    fit <- lsfit(x, y)$coefficients
    R10[k] <- log(fit[2]) / -TR
    M0[k] <- fit[1] / (1 - fit[2])
  }

  if (verbose)
    cat("  Calculating concentration...", fill=TRUE)
  theta <- dangle * pi/180
  CD <- conc <- matrix(NA, nvoxels, W)
  B <- (1 - exp(-TR * R10)) / (1 - cos(theta) * exp(-TR * R10))
  A <- (dyn.mat - dyn.mat[,1]) / M0 / sin(theta)
  R1t <- -(1/TR) * log((1 - (A+B)) / (1 - cos(theta) * (A+B)))
  conc <- (R1t - R10) / r1

  if (verbose)
    cat("  Reconstructing results...", fill=TRUE)
  R10.array <- M0.array <- array(NA, c(M,N,Z,1))
  R10.array[dyn.mask] <- R10
  M0.array[dyn.mask] <- M0
  conc.array <- R1t.array <- array(NA, c(M,N,Z,W))
  ## conc.array <- array(NA, c(M,N,Z,W))
  mask4D <- array(dyn.mask, c(M,N,Z,W))
  conc.array[mask4D] <- unlist(conc)
  R1t.array[mask4D] <- unlist(R1t)
  ## list(M0 = M0.array, R10 = R10.array, conc = conc.array)
  list(M0 = M0.array, R10 = R10.array, R1t = R1t.array, conc = conc.array)
}
