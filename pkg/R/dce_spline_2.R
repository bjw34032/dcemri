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
## Time-stamp: <>
## $Id$
##

dce.spline.onset <- function(data, time, input, time.input=time,
                             mask=array(1,dim(data)[1:(length(dim(data))-1)]),
                             k=4, knotpoints=25, knots=seq(min(time)-1e-3-(k-1)*(max(time)-min(time))/(knotpoints-k+1),1e-3+max(time)+k*(max(time)-min(time))/(knotpoints-k+1),(2e-3+max(time)-min(time))/(knotpoints-k+1)),
                             rw=2, nriters=500, thin=5, burnin=100,
                             nlr=TRUE, hyper.a=1e-5, hyper.b=1e-5,
                             epsilon.a=1, epsilon.b=1e-5, t0=0,
                             t0.compute=TRUE) {

  kk <- 0
  for (j in mask) {
    kk <- kk+1
    if (j==0) {
      for (i in 0:(length(time)-1))
        data[kk + i * prod(dim(mask))] <- 0
    }
  }
  data[is.na(data)] <- 0
  
  if (input=="FH" || input=="Weinmann") {
    if (input=="Weinmann") {
      D=0.1 ; a1=3.99 ; a2=4.78 ; m1=0.244 ; m2=0.0111
    }
    if (input=="FH") {
      D=1 ; a1=2.4 ; a2=0.62 ; m1=3 ; m2=0.016
    }
    input <- D * (a1*exp(-m1*time) + a2*exp(-m2*time))
    input.abl <- -D * (a1*m1*exp(-m1*time) + a2*m2*exp(-m2*time))
  } else {
    input.abl <- input[2:length(input)]-input[1:(length(input)-1)]
    input.abl <- input.abl / (time.input[2:length(input)] -
                              time.input[1:(length(input)-1)])
    input.abl <- c(0, input.abl)
  }

  ## function to make precision matrix for random walk
  R <- function(taux, rw) {
    RR <- matrix(0, nrow=length(taux)+rw, ncol=length(taux)+rw)
    if(rw==0) {
      for (i in 1:length(taux))
        RR[i,i] <- taux[i]
    }
    if(rw==1) {
      for (i in 1:length(taux)) {
        RR[i,i] <- RR[i,i] + taux[i]
        RR[i+1,i+1] <- RR[i+1,i+1] + taux[i]
        RR[i+1,i] <- RR[i+1,i] - taux[i]
        RR[i,i+1] <- RR[i,i+1] - taux[i]
      }
    }
    if(rw==2) {
      for (i in 1:length(taux)) {
        RR[i,i] <- RR[i,i]+taux[i]
        RR[i+1,i+1] <- RR[i+1,i+1] + 4*taux[i]
        RR[i+2,i+2] <- RR[i+2,i+2] + taux[i]
        RR[i+1,i] <- RR[i+1,i] - 2*taux[i]
        RR[i,i+1] <- RR[i,i+1] - 2*taux[i]
        RR[i+2,i+1] <- RR[i+2,i+1] - 2*taux[i]
        RR[i+1,i+2] <- RR[i+1,i+2] - 2*taux[i]
        RR[i+2,i] <- RR[i+2,i] + taux[i]
        RR[i,i+2] <- RR[i,i+2] + taux[i]
      }
    }
    return(RR)
  }
  
  ## redefinitions
  aifzeit <- time.input
  zeit <- time
  ## define B and A
  p <- length(knots) - k
  B <- splineDesign(knots, aifzeit, k, outer.ok=TRUE)
  if (sum(B[,dim(B)[2]] == 0) == dim(B)[1])
    B <- B[,-dim(B)[2]]
  if (sum(B[,1] == 0) == dim(B)[1])
    B <- B[,-1]
  p <- dim(B)[2]
  A <- matrix(0, nrow=length(zeit), ncol=length(aifzeit))
  ni <- zeit
  for (i in 1:length(zeit))
    for (j in 1:length(aifzeit))
      if (aifzeit[j] <= zeit[i])
        ni[i] <- j
  for (i in 1:length(zeit))
    for (j in 1:ni[i])
      A[i,j] <- input[1+ni[i]-j]
  A <- A * mean(diff(aifzeit))
  D <- A %*% B
  
  A.abl <- matrix(0, nrow=length(zeit), ncol=length(aifzeit))
  for (i in 1:length(zeit))
    for (j in 1:ni[i])
      A.abl[i,j] <- input.abl[1+ni[i]-j]
  A.abl <- A.abl * mean(diff(aifzeit))
  
  T <- length(zeit)
  dims <- dim(as.array(data))
  dims <- dims[-length(dims)]
  ddd <- 3
  if (length(dims) == 0) {
    dims <- c(1,1,1)
    ddd <- 0
  }
  if (length(dims) == 1) {
    dims <- c(1,1,dims)
    ddd <- 1
  }
  if (length(dims) == 2) {
    dims <- c(1,dims)
    ddd <- 2
  }
  data <- array(data, c(dims, T))
  tau <- array(1000, c(dims, p-rw, nriters))
  beta <- array(0, c(dims, p, nriters))
  MAX <- array(0, c(dims, nriters))
  D2 <- t(D) %*% D
  tauepsilon  <-  array(1000, c(dims, nriters))
  burnin <- min(burnin, nriters)
  TT <- max(T, p)
  
  result <- .C("dce_spline_run",
             as.integer(1),
             as.integer(burnin),
             as.integer(c(dims, T)),
             as.double(data),
             as.double(tau),
             as.double(tauepsilon),
             as.double(D),
             as.integer(rw),
             as.double(beta),
             as.double(c(hyper.a, hyper.b, epsilon.a, epsilon.b)),
             as.integer(p),
             as.double(1:TT),
             as.double(1:TT),
             as.double(1:TT),
             as.double(1:(TT^2)),
             as.double(1:(TT^2)),
             as.double(1:(TT^2)),
             as.double(1:(TT^2)),
             as.double(t(D)),
             PACKAGE="dcemri")
  
  tau <- array(result[[5]], c(dims, p-rw, nriters))
  beta <- array(result[[9]], c(dims, p, nriters))
  tauepsilon <- array(result[[6]], c(dims, nriters))
  
  ##print ("Burn in done")
  ##print(tau[,,,,1])
  ##print(beta[,,,,1])
  ##print(tauepsilon[,,,,1])

  result <- .C("dce_spline_run",
               as.integer(nriters),
               as.integer(thin),
               as.integer(c(dims, T)),
               as.double(data),
               as.double(tau),
               as.double(tauepsilon),
               as.double(D),
               as.integer(rw),
               as.double(beta),
               as.double(c(hyper.a, hyper.b, epsilon.a, epsilon.b)),
               as.integer(p),
               as.double(1:TT),
               as.double(1:TT),
               as.double(1:TT),
               as.double(1:(TT^2)),
               as.double(1:(TT^2)),
               as.double(1:(TT^2)),
               as.double(1:(TT^2)),
               as.double(t(D)),
               PACKAGE="dcemri")

  tau <- array(result[[5]], c(dims, p-rw, nriters))
  beta <- array(result[[9]], c(dims, p, nriters))
  tauepsilon <- array(result[[6]], c(dims, nriters))
  
  t0 <- rep(t0, nriters)
  
  if (t0.compute) {
    print ("Computing t0")
    
    q05 <- function(x)
      quantile(x, .005, na.rm=TRUE)
    med.na <- function(x)
      median(x, na.rm=TRUE)
    
    d <- array(NA, c(dims[1:3], T, nriters))
    du <- array(NA, dims[1:3])
    du2 <- array(NA, dims[1:3])
    for (x in 1:dims[1]) {
      for (y in 1:dims[2]) {
        for (z in 1:dims[3]) {
          for (j in 1:nriters)
            d[x,y,z,,j] <- D %*% beta[x,y,z,,j]
        }
      }
    }
    
    d1 <- apply(d, 1:4, q05)
    d2 <- apply(d, 1:4, med.na)
    
    for (x in 1:dims[1]) {
      for (y in 1:dims[2]) {
        for (z in 1:dims[3]) {
            du[x,y,z] <- min(which(d1[x,y,z,] > 0))
            if (is.na(du[x,y,z]))
              du[x,y,z] <- 1
          }
      }
    }
    
    beta.abl <- beta.abl2 <- rep(0, p)
    
    B2 <- splineDesign(knots, aifzeit, k-2)
    B2 <- B2[,(1:p)+1]
    
    for (j in 1:nriters) {
      for (x in 1:dims[1]) {
        for (y in 1:dims[2]) {
          for (z in 1:dims[3]) {
            beta.abl <- 1:p
            for (q in 1:(p-1))
              beta.abl[q] <- (beta[x,y,z,q+1,j] - beta[x,y,z,q,j]) * k / (knots[q+k+1] - knots[q+1])
            beta.abl[p] <- 0
            ABL2 <- A %*% B2 %*% beta.abl
            du2[x,y,z] <- zeit[du[x,y,z]] - d2[x,y,z,du[x,y,z]] / ABL2[du[x,y,z]]
          }
        }
      }
      t0[j] <- median(du2, na.rm=TRUE)
    }
  
    t0[t0 < 0] <- 0
  }
  
  print ("Done.")
  ## print(t0)
  t0 <- median(t0)
  return(t0)
}


