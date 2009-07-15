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

dce.spline <- function(data,time,input,time.input=time,mask=array(1,dim(data)[1:(length(dim(data))-1)]),k=4,knotpoints=25,
                       knots=seq(min(time)-1e-3-(k-1)*(max(time)-min(time))/(knotpoints-k+1),1e-3+max(time)+k*(max(time)-min(time))/(knotpoints-k+1),(2e-3+max(time)-min(time))/(knotpoints-k+1)),
                       rw=2,nriters=500,thin=5,burnin=100,nlr=TRUE,
                       hyper.a=1e-5,hyper.b=1e-5,epsilon.a=1,epsilon.b=1e-5,silent=FALSE) {
  
  kk<-0
  for (j in mask){
    kk<-kk+1
    if (j==0) {
      for (i in 0:(length(time)-1))
        data[kk+i*prod(dim(mask))] <- 0
    }
  }
  data[is.na(data)]<-0
  
  if (input=="FH"||input=="Weinmann") {
    if (input=="Weinmann") { D=0.1; a1=3.99; a2=4.78; m1=0.244; m2=0.0111 }
    if (input=="FH") { D=1; a1=2.4; a2=0.62; m1=3; m2=0.016 }
    input <- D*(a1*exp(-m1*time)+a2*exp(-m2*time))
  }
  
  ##function to make precision matrix for random walk
  R <- function(taux,rw) {
    RR <- matrix(0,nrow=length(taux)+rw,ncol=length(taux)+rw)
    if (rw==0) {
      for (i in 1:length(taux))
        RR[i,i] <- taux[i]
    }
    if (rw==1) {
      for (i in 1:length(taux)) {
        RR[i,i] <- RR[i,i]+taux[i]
        RR[i+1,i+1] <- RR[i+1,i+1]+taux[i]
        RR[i+1,i] <- RR[i+1,i]-taux[i]
        RR[i,i+1] <- RR[i,i+1]-taux[i]
      }
    }
    if (rw==2) {
      for (i in 1:length(taux)) {
        RR[i,i] <- RR[i,i]+taux[i]
        RR[i+1,i+1] <- RR[i+1,i+1]+4*taux[i]
        RR[i+2,i+2] <- RR[i+2,i+2]+taux[i]
        RR[i+1,i] <- RR[i+1,i]-2*taux[i]
        RR[i,i+1] <- RR[i,i+1]-2*taux[i]
        RR[i+2,i+1] <- RR[i+2,i+1]-2*taux[i]
        RR[i+1,i+2] <- RR[i+1,i+2]-2*taux[i]
        RR[i+2,i] <- RR[i+2,i]+taux[i]
        RR[i,i+2] <-RR [i,i+2]+taux[i]
      }
    }
    return(RR)
  }
  
  ##redefinitions
  aifzeit <- time.input
  zeit <- time
  
  ##define B and A
  p <- length(knots)-k
  B <- splineDesign(knots, aifzeit, k, outer.ok=TRUE)
  if (sum(B[,dim(B)[2]]==0) == dim(B)[1])
    B <- B[,-dim(B)[2]]
  if (sum(B[,1]==0) == dim(B)[1])
    B <- B[,-1]
  p <- dim(B)[2]
  A <- matrix(0, nrow=length(zeit), ncol=length(aifzeit))
  ni <- zeit
  for (i in 1:length(zeit)) {
    for (j in 1:length(aifzeit)) {
      if (aifzeit[j]<=zeit[i])
        ni[i] <- j
    }
  }
  for (i in 1:length(zeit)) {
    for (j in 1:ni[i])
      A[i,j] <- input[1+ni[i]-j]
  }
  A <- A*mean(diff(aifzeit))
  D <- A%*%B
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
  data <- array(data, c(dims,T))
  
  tau <- array(1000, c(dims,p-rw,nriters))
  beta <- array(0, c(dims,p,nriters))
  MAX <- array(0, c(dims,nriters))
  D2 <- t(D)%*%D
  tauepsilon <- array(1000, c(dims,nriters))
  burnin <- min(burnin, nriters)
  
  if(!silent) { print(dims) }
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
               as.double(1:T),
               as.double(1:T),
               as.double(1:T),
               as.double(1:(T^2)),
               as.double(1:(T^2)),
               as.double(1:(T^2)),
               as.double(1:(T^2)),
               as.double(t(D)),
               as.integer(silent),
               PACKAGE="dcemri")
  
  tau <- array(result[[5]], c(dims,p-rw,nriters))
  beta <- array(result[[9]], c(dims,p,nriters))
  tauepsilon <- array(result[[6]], c(dims,nriters))
  
  if(!silent) { print("Burn in done") }
  
  ##print(tau[,,,,1])
  ##print(beta[,,,,1])
  ##print(tauepsilon[,,,1])
  
  result <- .C("dce_spline_run",
               as.integer(nriters),
               as.integer(thin),
               as.integer(c(dims,T)),
               as.double(data),
               as.double(tau),
               as.double(tauepsilon),
               as.double(D),
               as.integer(rw),
               as.double(beta),
               as.double(c(hyper.a, hyper.b, epsilon.a, epsilon.b)),
               as.integer(p),
               as.double(1:T),
               as.double(1:T),
               as.double(1:T),
               as.double(1:(T^2)),
               as.double(1:(T^2)),
               as.double(1:(T^2)),
               as.double(1:(T^2)),
               as.double(t(D)),
               as.integer(silent),
               PACKAGE="dcemri")
  
  tau <- array(result[[5]], c(dims,p-rw,nriters))
  beta <- array(result[[9]], c(dims,p,nriters))
  tauepsilon <- array(result[[6]], c(dims,nriters))
  
  if(!silent) { print("Iterations done") }
  
  MAX0 <- apply(beta, c(1:3,5), max)
  MAX <- apply(MAX0, c(1:3), median)
  
  parhat <- NULL
  
  modelTC0 <- function(time, F, E, ve) {
    kep <- E*F/ve
    erg <- E*exp(-kep*(time))
    erg <- erg*F
    eval(erg)
  }
  
  model <- function(time, F, E, TC, ve) {
    TC2<-2*TC
    if (TC < 0) { TC2 <- 0 }
    kep <- E*F/ve
    erg <- E*exp(-kep*(time-TC))
    ## erg[time<TC]<-1
    erg[time<TC] <- 1 - time[time<TC2]*(1-E) / TC
    erg <- erg*F
    if (TC < 0)
      erg <- rep(-10^16, length(time))
    ## erg<-expression(erg)
    eval(erg)
  }
  
  fcn <- function(p, time, x, N.Err, fcall, jcall)
    (x - do.call("fcall", c(list(time = time), as.list(p))))
  fcn.jac <- function(p, time, x, N.Err, fcall, jcall) {
    N.Err <- rep(N.Err, length(p))
    -do.call("jcall", c(list(time = time), as.list(p)))/N.Err
  }
  
  if (nlr=="med"||nlr=="mean") {
    if(!silent) { print("Starting non-linear regression") }
    
    X <- dims[1]
    Y <- dims[2]
    Z <- dims[3]
    
    parhat <- list(Fp = array(0,c(X,Y,Z)),
                   E = array(0,c(X,Y,Z)),
                   ktrans = array(0,c(X,Y,Z)),
                   kep = array(0,c(X,Y,Z)),
                   ve = array(0,c(X,Y,Z)),
                   TC = array(0,c(X,Y,Z)),
                   vp = array(0,c(X,Y,Z))
                   )
    for (i in 1:X) {
      for (j in 1:Y) {
        for (m in 1:Z) {
          if (nlr=="med")
            me5 <- apply(beta[i,j,m,,],1,median) %*% t(B)
          if (nlr=="mean")
            me5 <- apply(beta[i,j,m,,],1,mean) %*% t(B)
          kt <- max(me5)
          start <- which(kt==me5)
          try(tmp <- nls.lm(par=list("E"=.6,"F"=MAX[i,j,m,],"TC"=0,"ve"=.05), 
                            fn=fcn, fcall=model, #jcall= modelJ,
                            time=zeit, x=me5, N.Err=sqrt(300),
                            control=list(nprint=1,ftol=10^-20)))
          
          if (!(parhat$TC>0))
            try(tmp<-nls.lm(par=list("E"=.6,"F"=MAX[i,j,m,],"TC","ve"=.05), 
                            fn=fcn, fcall=modelTC0, #jcall= modelJ,
                            time=zeit, x=me5, N.Err=sqrt(300),
                            control=list(nprint=1,ftol=10^-20)))
          
          parhat$E[i,j,m] <- tmp$par$E
          parhat$Fp[i,j,m] <- tmp$par$F
          parhat$ve[i,j,m] <- tmp$par$ve
          parhat$TC[i,j,m] <- tmp$par$TC
          parhat$ktrans[i,j,m] <- tmp$par$E * tmp$par$F
          parhat$kep[i,j,m] <- parhat$ktrans[i,j,m] / tmp$par$ve
          parhat$vp[i,j,m] <- tmp$par$TC * tmp$par$F
        }
      }
    }
    if (ddd==0) {
      parhat$E <- parhat$E[1,1,1]
      parhat$Fp <- parhat$Fp[1,1,1]
      parhat$TC <- parhat$TC[1,1,1]
      parhat$ve <- parhat$ve[1,1,1]
      parhat$vp <- parhat$vp[1,1,1]
      parhat$ktrans <- parhat$ktrans[1,1,1]
      parhat$kep <- parhat$kep[1,1,1]
    }
    if (ddd==1) {
      parhat$E <- parhat$E[1,1,]
      parhat$Fp <- parhat$Fp[1,1,]
      parhat$TC <- parhat$TC[1,1,]
      parhat$ve <- parhat$ve[1,1,]
      parhat$vp <- parhat$vp[1,1,]
      parhat$ktrans <- parhat$ktrans[1,1,]
      parhat$kep <- parhat$kep[1,1,]
    }
    if (ddd==2) {
      parhat$E <- parhat$E[1,,]
      parhat$Fp <- parhat$Fp[1,,]
      parhat$TC <- parhat$TC[1,,]
      parhat$ve <- parhat$ve[1,,]
      parhat$vp <- parhat$vp[1,,]
      parhat$ktrans <- parhat$ktrans[1,,]
      parhat$kep <- parhat$kep[1,,]
    }
    nlr <- FALSE
  }
  if (nlr) {
    if(!silent)print ("Starting non-linear regression")
    X <- dims[1]
    Y <- dims[2]
    Z <- dims[3]

    parhat <- list(Fp = array(0,c(X,Y,Z,nriters)),
                   E = array(0,c(X,Y,Z,nriters)),
                   ktrans = array(0,c(X,Y,Z,nriters)),
                   kep = array(0,c(X,Y,Z,nriters)),
                   ve = array(0,c(X,Y,Z,nriters)),
                   TC = array(0,c(X,Y,Z,nriters)),
                   vp = array(0,c(X,Y,Z,nriters))
                   )

    for (i in 1:X)
      for (j in 1:Y)
        for (m in 1:Z) {
          for (k in 1:nriters) {
            me5 <- beta[i,j,m,,k]%*%t(B)
            kt <- max(me5)
            start <- which(kt==me5)
            try(tmp <- nls.lm(par=list("E"=1,"F"=1,"TC"=0.1,"ve"=.1), 
                              fn=fcn, fcall = model, #jcall= modelJ,
                              time=zeit, x=me5, N.Err=sqrt(300),
                              control=list(nprint=-1,ftol=10^-20)))
            parhat$E[i,j,m,k] <- tmp$par$E
            parhat$Fp[i,j,m,k] <- tmp$par$F
            parhat$ve[i,j,m,k] <- tmp$par$ve
            parhat$TC[i,j,m,k] <- tmp$par$TC
            parhat$ktrans[i,j,m,k] <- tmp$par$E * tmp$par$F
            parhat$kep[i,j,m,k] <- parhat$ktrans[i,j,m,k] / tmp$par$ve
            parhat$vp[i,j,m,k] <- tmp$par$TC * tmp$par$F
          }
        }
    if (ddd==0) {
      parhat$E <- parhat$E[1,1,1,]
      parhat$Fp <- parhat$Fp[1,1,1,]
      parhat$TC <- parhat$TC[1,1,1,]
      parhat$ve <- parhat$ve[1,1,1,]
      parhat$vp <- parhat$vp[1,1,1,]
      parhat$ktrans <- parhat$ktrans[1,1,1,]
      parhat$kep <- parhat$kep[1,1,1,]
    }
    if (ddd==1) {
      parhat$E <- parhat$E[1,1,,]
      parhat$Fp <- parhat$Fp[1,1,,]
      parhat$TC <- parhat$TC[1,1,,]
      parhat$ve <- parhat$ve[1,1,,]
      parhat$vp <- parhat$vp[1,1,,]
      parhat$ktrans <- parhat$ktrans[1,1,,]
      parhat$kep <- parhat$kep[1,1,,]
    }
    if (ddd==2) {
      parhat$E <- parhat$E[1,,,]
      parhat$Fp <- parhat$Fp[1,,,]
      parhat$TC <- parhat$TC[1,,,]
      parhat$ve <- parhat$ve[1,,,]
      parhat$vp <- parhat$vp[1,,,]
      parhat$ktrans <- parhat$ktrans[1,,,]
      parhat$kep <- parhat$kep[1,,,]
    }
  }
  if (ddd==0) {
    ##tau <- tau[1,1,1,,]
    ##beta <- beta[1,1,1,,]
    ##tauepsilon <- tauepsilon[1,1,1,,]
    MAX <- MAX[1,1,1]
  }
  if (ddd==2) {
    ##tau <- tau[1,1,,,]
    ##beta <- beta[1,1,,,]
    ##tauepsilon <- tauepsilon[1,1,,,]
    MAX <- MAX[1,1,]
  }
  if (ddd==3) {
    ##tau <- tau[1,,,,]
    ##beta <- beta[1,,,,]
    ##tauepsilon <- tauepsilon[1,,,,]
    MAX <- MAX[1,,]
  }

  if (!silent) { print ("Done.") }
  if (nlr) {
    ktrans <- apply(parhat$ktrans, 1:ddd, median, na.rm=TRUE)
    kep <- apply(parhat$kep, 1:ddd, median, na.rm=TRUE)
    vp <- apply(parhat$vp, 1:ddd, median, na.rm=TRUE)
    sigma2 <- 1/apply(tauepsilon, 1:3, median, na.rm=TRUE)
  }
  result <- list(fp.hat = MAX,
                 tau = tau,
                 beta = beta,
                 tauepsilon = tauepsilon,
                 A = A,
                 B = B,
                 D = D,
                 knots = knots,
                 estimates = parhat,
                 fp.hat.samples = MAX0)
  return(result)
}

perfSpline <- function(data,  time,aif, time.input, searchint=.25,
                       do.plot=TRUE) {
  require("splines")
  interpSpline(aif ~ time.input, bSpline=TRUE) -> XX
  time.input2 <- seq(0, max(time.input), 0.25)
  B <- splineDesign(XX$knots, time.input2, XX$order, outer.ok=TRUE)
  aif2 <- B %*% XX$coefficients

  bb <- c()
  for (t0 in seq(0, 5*searchint, searchint)) {
    tau <- min(which(time[1:25]>t0))
    rest <- dce.spline(data[1:25+tau-1], time[1:25+tau-1]-t0, aif2,
                       time.input2,nlr=FALSE, knotpoints=30, rw=2,
                       nriters=500, hyper.a=1e-10, hyper.b=1e-10,
                       epsilon.b=1e-10, silent=TRUE)
    beta <- apply(rest$beta, c(3,4), median)
    bb <- c(bb, (max((rest$B %*% beta[1,])[1:50])))
    ## print(bb)
  }
  bc <- bb[6]
  while(bc < 3*mean(bb[1:3])) {
    t0=t0+searchint
    tau <- min(which(time[1:25] > t0))
    rest <- dce.spline(data[1:25+tau-1], time[1:25+tau-1]-t0, aif2,
                       time.input2, nlr=FALSE, knotpoints=30, rw=2,
                       nriters=500, hyper.a=1e-10, hyper.b=1e-10,
                       epsilon.b=1e-10, silent=TRUE)
    beta <- apply(rest$beta,c(3,4),median)
    bc <- max((rest$B%*%beta[1,])[1:50])
    bb <- c(bb,bc)
    ## print(bb)
  }

  tt <- seq(0, (length(bb)-1)*searchint, searchint)
  a <- lm(bb[1:3] ~ tt[1:3])
  b <- lm(bb[length(bb)-(2:0)] ~ tt[length(bb)-(2:0)])
  t0 <- ((b$coefficients[1] - a$coefficients[1]) /
         (a$coefficients[2]-b$coefficients[2]))

  if (do.plot) {
    plot(tt, bb, xlab="Delay (seconds)", ylab="MBF")
    lines(tt, a$coefficients[1] + a$coefficients[2]*tt)
    lines(tt, b$coefficients[1] + b$coefficients[2]*tt)
  }
  if (t0 < 0)
    t0 <- 0
  tau <- min(which(time[1:25] > t0))
  rest <- dce.spline(t(data[1:25+tau-1]), time[1:25+tau-1]-t0, aif2,
                     time.input2, nlr=FALSE, knotpoints=30, rw=2,
                     nriters=500, hyper.a=1e-10, hyper.b=1e-10,
                     epsilon.b=1e-10, silent=TRUE)
  beta <- apply(rest$beta, c(3,4), median)
  bc=max(rest$B %*% beta[1,])
  if (do.plot)
    plot(time.input2, rest$B %*% beta[1,], type="l",
         xlab="Time (seconds)", ylab="Impulse Response Function")
  
  list("mbf"=bc, "t0"=t0, "fit"=rest)
}
