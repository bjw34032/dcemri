##
##
## Copyright (c) 2010 Brandon Whitcher and Volker Schmid
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
## $Id:$
##

aifParameters <- function(type, user=NULL) {
  switch(type,
         tofts.kermode = list(D=0.1, a1=3.99, m1=0.144, a2=4.78, m2=0.0111),
         fritz.hansen = list(D=1.0, a1=2.4, m1=3.0, a2=0.62, m2=0.016),
         orton.exp = list(D=1.0, AB=323, muB=20.2, AG=1.07, muG=0.172),
         orton.cos = list(D=1.0, AB=2.84, muB=22.8, AG=1.36, muG=0.171),
         user = user,
         empirical = user,
         stop("AIF parameters must be specified!"))
}

compartmentalModel <- function(type) {
  switch(type,
         weinmann =
         function(time, theta, p) {
           ## Convolution of Tofts & Kermode AIF with single-compartment model
           th1 <- theta[1]
           th3 <- theta[2]
           erg <- p$D * exp(th1) * ((p$a1 / (p$m1 - exp(th3))) * (exp(-(time * exp(th3))) - exp(-(time * p$m1))) + (p$a2 / (p$m2 - exp(th3))) * (exp(-(time * exp(th3))) - exp(-(time * p$m2))))
           erg[time <= 0] <- 0
           return(erg)
         },
         extended =
         function(time, theta, p) {
           ## Extended Tofts & Kermode model including the concentration of
           ## contrast agent in the blood plasma (vp)
           Cp <- function(t, p) {
             p$D * (p$a1 * exp(-p$m1 * t) + p$a2 * exp(-p$m2 * t))
           }
           th0 <- theta[1]
           th1 <- theta[2]
           th3 <- theta[3]
           erg <- exp(th0) * Cp(time, p) + p$D * exp(th1) *
             ((p$a1 / (p$m1 - exp(th3))) *
              (exp(-(time * exp(th3))) - exp(-(time * p$m1))) +
              (p$a2 / (p$m2 - exp(th3))) *
              (exp(-(time * exp(th3))) - exp(-(time * p$m2))))
           erg[time <= 0] <- 0
           return(erg)
         },
         kety.orton.exp =
         function(time, theta, p) {
           ## Kety model using the exponential AIF from Matthew Orton (ICR)
           ktrans <- exp(theta[1])
           kep <- exp(theta[2])
           T1 <- p$AB * kep / (kep - p$muB)
           T2 <- time * exp(-p$muB * time) -
             (exp(-p$muB * time) - exp(-kep * time)) / (kep - p$muB)
           T3 <- p$AG * kep
           T4 <- (exp(-p$muG * time) - exp(-kep * time)) / (kep - p$muG) -
             (exp(-p$muB * time) - exp(-kep * time)) / (kep - p$muB)
           erg <- ktrans * (T1 * T2 + T3 * T4)
           erg[time <= 0] <- 0
           return(erg)
         },
         orton.exp =
         function(time, theta, p) {
           ## Extended model using the exponential AIF from Matthew Orton (ICR)
           Cp <- function(t, p) {
             p$AB * t * exp(-p$muB * t) + p$AG *
               (exp(-p$muG * t) - exp(-p$muB * t))
           }
           vp <- exp(theta[1])
           ktrans <- exp(theta[2])
           kep <- exp(theta[3])
           T1 <- p$AB * kep / (kep - p$muB)
           T2 <- time * exp(-p$muB * time) -
             (exp(-p$muB * time) - exp(-kep * time)) / (kep - p$muB)
           T3 <- p$AG * kep
           T4 <- (exp(-p$muG * time) - exp(-kep * time)) /
             (kep - p$muG) - (exp(-p$muB * time) - exp(-kep * time)) /
               (kep - p$muB)
           erg <- vp * Cp(time, p) + ktrans * (T1 * T2 + T3 * T4)
           erg[time <= 0] <- 0
           return(erg)
         },
         kety.orton.cos =
         function(time, theta, p) {
           ## Extended model with the raised cosine AIF from Matthew Orton (ICR)
           A2 <- function(time, alpha, p) {
             (1 - exp(-alpha * time)) / alpha -
               (alpha * cos(p$muB * time) + p$muB * sin(p$muB * time) -
                alpha * exp(-alpha * time)) / (alpha^2 + p$muB^2)
           }
           ktrans <- exp(theta[1])
           kep <- exp(theta[2])
           tB <- 2 * pi / p$muB
           cp <- ifelse(time <= tB,
                        p$aB * (1 - cos(p$muB*time)) + p$aB * p$aG * A2(time, p$muG, p),
                        p$aB * p$aG * A2(time, p$muG, p) * exp(-p$muB * (time - tB)))
           erg <- ifelse(time <= tB,
                         p$aB * p$aG * ktrans / (kep - p$muG) * ((A2(time, p$muG, p) + (kep - p$muG) / p$aG - 1) * A2(time, kep, p)),
                         p$aB * p$aG * ktrans / (kep - p$muG) * (A2(tB, p$muG, p) * exp(-p$muB * (time - tB)) + ((kep - p$muG) / p$aG - 1) * A2(tB, kep, p) * exp(-kep * (time - tB))))
           erg[time <= 0] <- 0
           return(erg)
         },
         orton.cos =
         function(time, theta, p) {
           ## Extended model with the raised cosine AIF from Matthew Orton (ICR)
           A2 <- function(time, alpha, p) {
             (1 - exp(-alpha * time)) / alpha -
               (alpha * cos(p$muB * time) + p$muB * sin(p$muB * time) -
                alpha * exp(-alpha * time)) / (alpha^2 + p$muB^2)
           }
           vp <- exp(theta[1])
           ktrans <- exp(theta[2])
           kep <- exp(theta[3])
           tB <- 2 * pi / p$muB
           cp <- ifelse(time <= tB,
                        p$aB * (1 - cos(p$muB * time)) + p$aB * p$aG * A2(time, p$muG, p),
                        p$aB * p$aG * A2(time, p$muG, p) * exp(-p$muB * (time - tB)))
           erg <- ifelse(time <= tB,
                         vp * cp + p$aB * p$aG * ktrans / (kep - p$muG) * ((A2(time, p$muG, p) + (kep - p$muG) / p$aG - 1) * A2(time, kep, p)),
                         vp * cp + p$aB * p$aG * ktrans / (kep - p$muG) * (A2(tB, p$muG, p) * exp(-p$muB * (time - tB)) + ((kep - p$muG) / p$aG - 1) * A2(tB, kep, p) * exp(-kep * (time - tB))))
           erg[time <= 0] <- 0
           return(erg)
         },
         weinmann.empirical =
         function(time, theta, aif) {
           th1 <- theta[1]
           th3 <- theta[2]
           tsec <- seq(min(time*60), ceiling(max(time*60)), by=1)
           aif.new <- approx(time*60, aif, tsec)$y
           erg <- approx(tsec, expConv(aif.new, exp(th1)/60, exp(th3)/60), time*60)$y
           erg[time <= 0] <- 0
           return(erg)
         },
         extended.empirical =
         function(time, theta, aif) {
           th0 <- theta[1]
           th1 <- theta[2]
           th3 <- theta[3]
           tsec <- seq(min(time*60), ceiling(max(time*60)), by=1)
           aif.new <- approx(time*60, aif, tsec)$y
           erg <- approx(tsec, exp(th0) * aif.new + expConv(aif.new, exp(th1)/60, exp(th3)/60), time*60)$y
           erg[time <= 0] <- 0
           return(erg)
         })
}

expConv <- function(input, k1, k2) {
  n <- length(input)
  convolution <- rep(0, n)
  ek2 <- exp(-k2)
  if (ek2 == 1) {
    convolution <- k1 * cumsum(input)
  } else {
    prev <- 0
    for (i in 1:n) {
      prev <- prev * ek2 + k1 * input[i] * (1 - ek2) / k2
      convolution[i] <- prev
    }
  }
  return(convolution)
}
  
