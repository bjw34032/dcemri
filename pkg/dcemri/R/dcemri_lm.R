##
##
## Copyright (c) 2009, Brandon Whitcher and Volker Schmid
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
##
## $Id$
##

dcemri.lm <- function(conc, time, img.mask, model="extended", aif=NULL,
                      nprint=0, user=NULL, verbose=FALSE, ...) {
  ## dcemri.lm - a function for fitting 1-compartment PK models to
  ## DCE-MRI images
  ##
  ## authors: Volker Schmid, Brandon Whitcher
  ##
  ## input:
  ##        conc: array of Gd concentration,
  ##        time: timepoints of aquisition,
  ##        img.mask: array of voxels to fit,
  ##        D(=0.1): Gd dose in mmol/kg,
  ##        model: AIF... "weinmann" or "parker",
  ##        update: re-do given parameter maps where parameters not fitted,
  ##        ktransmap, kepmap, ktranserror, keperror: given parameter maps.
  ##
  ## output: list with ktrans, kep, ve, std.error of ktrans and kep
  ##         (ktranserror and keperror)
  ##
  model <- switch(model,
                  weinmann = "weinmann",
                  extended = "extended",
                  orton.exp = "orton.exp",
                  orton.cos = "orton.cos",
                  stop("Unknown model ", model, call.=FALSE))
  aif <- switch(model,
                weinmann = ,
                extended = {
                  if (is.null(aif)) {
                    "tofts.kermode"
                  } else {
                    switch(aif,
                           tofts.kermode="tofts.kermode",
                           fritz.hansen="fritz.hansen",
                           stop("Only aif=\"tofts.kermode\" or aif=\"fritz.hansen\" acceptable aifs for model=\"weinmann\" or model=\"extended\"", call.=FALSE)
                           )
                  }
                },
                orton.exp = {
                  if (is.null(aif)) {
                    "orton.exp"
                  } else {
                    switch(aif,
                           orton.exp="orton.exp",
                           user="user",
                           stop("Only aif=\"orton.exp\" or aif=\"user\" acceptable aifs for model=\"orton.exp\""), call.=FALSE)
                  }
                },
                orton.cos = {
                  if (is.null(aif)) {
                    "orton.cos"
                  } else {
                    switch(aif,
                           orton.cos="orton.cos",
                           user="user",
                           stop("Only aif=\"orton.cos\" or aif=\"user\" acceptable aifs for model=\"orton.cos\""), call.=FALSE)
                  }
                },
                stop("Unknown model: " + model, call.=FALSE))
  
  require("minpack.lm")
  
  mod <- model
  nvoxels <- sum(img.mask)
  I <- nrow(conc)
  J <- ncol(conc)
  K <- nsli(conc)
  
  if (!is.numeric(dim(conc))) {
    I <- J <- K <- 1
  } else {
    if (length(dim(conc)) == 2)
      J <- K <- 1
  }

  if (verbose) cat("  Deconstructing data...", fill=TRUE)
  conc.mat <- matrix(conc[img.mask], nvoxels)
  conc.mat[is.na(conc.mat)] <- 0

  switch(aif,
         tofts.kermode = {
           D <- 0.1; a1 <- 3.99; a2 <- 4.78; m1 <- 0.144; m2 <- 0.0111
         },
         fritz.hansen= {
           D <- 1; a1 <- 2.4; a2 <- 0.62; m1 <- 3.0; m2 <- 0.016
         },
         orton.exp = {
           D <- 1; AB <- 323; muB <- 20.2; AG <- 1.07; muG <- 0.172
         },
         orton.cos = {
           D <- 1; aB <- 2.84; muB <- 22.8; aG <- 1.36; muG <- 0.171
         },
         user = {
           if (verbose) cat("  User-specified AIF parameters...", fill=TRUE);
           D <- try(user$D); AB <- try(user$AB); aB <- try(user$aB);
           muB <- try(user$muB); AG <- try(user$AG); aG <- try(user$aG); 
           muG <- try(user$muG)
         },
         print("WARNING: AIF parameters must be specified!"))
  
  model.weinmann <- function(time, th1, th3, ...) {
    ## Convolution of Tofts & Kermode AIF with single-compartment model
    erg <- D * exp(th1) * ((a1 / (m1 - exp(th3))) *
                           (exp(-(time * exp(th3))) - exp(-(time * m1))) +
                           (a2 / (m2 - exp(th3))) *
                           (exp(-(time * exp(th3))) - exp(-(time * m2))))
    erg[time <= 0] <- 0
    return(erg)
  }

  model.extended <- function(time, th0, th1, th3, ...) {
    ## Extended Tofts & Kermode model including the concentration of
    ## contrast agent in the blood plasma (vp)
    Cp <- function(tt, D, a1, a2, m1, m2)
      D * (a1 * exp(-m1 * tt) + a2 * exp(-m2 * tt))

    erg <- exp(th0) * Cp(time, D, a1, a2, m1, m2) +
      D * exp(th1) * ((a1 / (m1 - exp(th3))) *
                      (exp(-(time * exp(th3))) - exp(-(time * m1))) +
                      (a2 / (m2 - exp(th3))) *
                      (exp(-(time * exp(th3))) - exp(-(time * m2))))
    erg[time <= 0] <- 0
    return(erg)
  }
  
  model.orton.exp <- function(time, th0, th1, th3, ...) {
    ## Extended model using the exponential AIF from Matthew Orton (ICR)
    
    Cp <- function(tt, ...) 
      AB * tt * exp(-muB * tt) + AG * (exp(-muG * tt) - exp(-muB * tt))
    
    vp <- exp(th0)
    ktrans <- exp(th1)
    kep <- exp(th3)
    
    T1 <- AB * kep / (kep - muB)
    T2 <- time * exp(-muB * time) -
      (exp(-muB * time) - exp(-kep * time)) / (kep - muB)
    T3 <- AG * kep
    T4 <- (exp(-muG * time) - exp(-kep * time)) / (kep - muG) -
      (exp(-muB * time) - exp(-kep * time)) / (kep - muB)

    erg <- vp * Cp(time) + ktrans * (T1 * T2 + T3 * T4)
    erg[time <= 0] <- 0
    return(erg)
  }
  
  model.orton.cos <- function(time, th0, th1, th3, ...) {
    ## Extended model using the raised cosine AIF from Matthew Orton (ICR)

    A2 <- function(t, alpha, ...) {
      (1 - exp(-alpha*time)) / alpha - (alpha*cos(muB*time) + muB*sin(muB*time) - alpha*exp(-alpha*time)) / (alpha^2 + muB^2)
    }

    vp <- exp(th0)
    ktrans <- exp(th1)
    kep <- exp(th3)

    tB <- 2*pi/muB
    cp <- ifelse(time <= tB,
                 aB * (1 - cos(muB*time)) + aB * aG * A2(time, muG),
                 aB * aG * A2(time, muG) * exp(-muB * (time - tB)))
    erg <- ifelse(time <= tB,
                  vp * cp + aB * aG * ktrans / (kep - muG) * ((A2(time,muG) + (kep - muG) / aG - 1) * A2(time,kep)),
                  vp * cp + aB * aG * ktrans / (kep - muG) * (A2(tB,muG) * exp(-muB * (time - tB)) + ((kep - muG) / aG - 1) * A2(tB,kep) * exp(-kep * (time - tB))))
    erg[time <= 0] <- 0
    return(erg)
  }

  switch(mod,
         weinmann = {
           model <- model.weinmann
           func <- function(theta, signal, time, ...)
             signal - model(time, theta[1], theta[2])
           guess <- c("th1"=0, "th3"=0.1)
         },
         extended = {
           model <- model.extended
           func <- function(theta, signal, time, ...)
             signal - model.extended(time, theta[1], theta[2], theta[3])
           guess <- c("th0"=-1, "th1"=0, "th3"=0.1)
           Vp <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
        },
         orton.exp = {
           model <- model.orton.exp
           func <- function(theta, signal, time, ...)
             signal - model(time, theta[1], theta[2], theta[3], ...)
           guess <- c("th0"=-1, "th1"=0, "th3"=0.1)
           Vp <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
         },
         orton.cos = {
           model <- model.orton.cos
           func <- function(theta, signal, time, ...)
             signal - model(time, theta[1], theta[2], theta[3], ...)
           guess <- c("th0"=-1, "th1"=0, "th3"=0.1)
           Vp <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
         },
         stop("Model/AIF is not supported."))

  ktrans <- kep <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  sse <- rep(NA, nvoxels)

  if (verbose) cat("  Estimating the kinetic parameters...", fill=TRUE)
  for(k in 1:nvoxels) {
    fit <- nls.lm(par=guess, fn=func, control=list(nprint=nprint),
                  signal=conc.mat[k,], time=time)
    if(fit$info %in% 1:4) {
      ktrans$par[k] <- exp(fit$par["th1"])
      kep$par[k] <- exp(fit$par["th3"])
      ktrans$error[k] <- sqrt(fit$hessian["th1","th1"])
      kep$error[k] <- sqrt(fit$hessian["th3","th3"])
      if(mod %in% c("extended","orton.exp","orton.cos")) {
        Vp$par[k] <- exp(fit$par["th0"])
        Vp$error[k] <- sqrt(fit$hessian["th0","th0"])
        sse[k] <-
          sum((conc.mat[k,] - model(time, fit$par["th0"], fit$par["th1"],
                                    fit$par["th3"]))^2)
      } else {
        sse[k] <- sum((conc.mat[k,] - model(time, fit$par["th1"],
                                            fit$par["th3"]))^2)
      }
    }
  }

  if (verbose) cat("  Reconstructing results...", fill=TRUE)
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- ktrans$par
  B[img.mask] <- ktrans$error
  ktrans.out <- list(par=A, error=B)
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- kep$par
  B[img.mask] <- kep$error
  kep.out <- list(par=A, error=B)
  if(mod %in% c("extended","orton.exp","orton.cos")) {
    A <- B <- array(NA, c(I,J,K))
    A[img.mask] <- Vp$par
    B[img.mask] <- Vp$error
    Vp.out <- list(par=A, error=B)
  }
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- sse
  sse.out <- A
  
  if(mod %in% c("extended","orton.exp","orton.cos"))
    list(ktrans=ktrans.out$par, kep=kep.out$par, ktranserror=ktrans.out$error,
         keperror=kep.out$error, ve=ktrans.out$par/kep.out$par, vp=Vp.out$par,
         vperror=Vp.out$error, sse=sse.out, time=time)
  else
    list(ktrans=ktrans.out$par, kep=kep.out$par, ktranserror=ktrans.out$error,
         keperror=kep.out$error, ve=ktrans.out$par/kep.out$par, sse=sse.out,
         time=time)
}


