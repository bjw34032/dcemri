##
##
## Copyright (c) 2009,2010 Brandon Whitcher and Volker Schmid
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
## $Id:$
##

#############################################################################
## Regional Blood Volume (rCBV) using numeric integration
#############################################################################

rCBV <- function(Ct, Ca, time, Hf=1, rho=1) {
  trapezoid <- function(x, y) {
    ## Trapezoidal rule
    ## x values need not be equally spaced
    sum(diff(x) * (y[-1] + y[-length(y)])) / 2
  }
  ## Hf  = hematocrit factor = 0.45 (?)
  ## rho = density of brain tissue = 1.04 g/mL (?)
  AUC.Ct <- trapezoid(time, Ct)
  AUC.Ca <- trapezoid(time, Ca)
  (Hf / rho) * (AUC.Ct / AUC.Ca)
}

#############################################################################
## setGeneric("rCBV.fast")
#############################################################################

setGeneric("rCBV.fast", function(signal, ...) standardGeneric("rCBV.fast"))
setMethod("rCBV.fast", signature(signal="array"),
          function(signal, mask, aif, time, multicore=FALSE, verbose=FALSE) 
	    .dcemriWrapper("rCBV.fast", signal, mask, aif, time, multicore,
                           verbose))

#############################################################################
## rCBV.fast()
#############################################################################

.rCBV.fast <- function(signal, mask, aif, time, multicore=FALSE,
                       verbose=FALSE) {
  if (length(dim(signal)) != 4) { # Check signal is a 4D array
    stop("Data must be a 4D array.")
  }
  if (!is.logical(mask)) { # Check mask is logical
    stop("Mask must be logical.")
  }
  if (length(time) != ntim(signal)) { # Check signal and time match
    stop("Data must have the same length as the acquisition times.")
  }
  X <- nrow(signal)
  Y <- ncol(signal)
  Z <- nsli(signal)
  nvoxels <- sum(mask)
  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  signal.mat <- matrix(signal[mask], nrow=nvoxels)
  signal.list <- vector("list", nvoxels)
  for (k in 1:nvoxels) {
    signal.list[[k]] <- signal.mat[k,]
  }
  if (verbose) {
    cat("  Calculating rCBV...", fill=TRUE)
  }
  if (multicore && require("multicore")) {
    rcbv.list <- mclapply(signal.list, function(x) {
      rCBV(x, Ca=aif, time=time)
    })
  } else {
    rcbv.list <- lapply(signal.list, function(x) {
      rCBV(x, Ca=aif, time=time)
    })
  }
  rm(signal.list) ; gc()
  rcbv <- unlist(rcbv.list)
  rm(rcbv.list) ; gc()
  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  rcbv.array <- array(NA, dim(signal)[1:3])
  rcbv.array[mask] <- rcbv
  return(rcbv.array)
}

