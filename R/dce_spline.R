dce.spline <- function(data, time, input, time.input=time,
                       mask=array(1,dim(data)[1:(length(dim(data))-1)]),
                       k=4, knotpoints=25,
                       knots=seq(min(time)-1e-3-(k-1)*(max(time)-min(time))/(knotpoints-k+1),1e-3+max(time)+k*(max(time)-min(time))/(knotpoints-k+1),(2e-3+max(time)-min(time))/(knotpoints-k+1)),
                       rw=2, nriters=500, thin=5, burnin=100, nlr=TRUE,
                       hyper.a=1e-5, hyper.b=1e-5, epsilon.a=1,
                       epsilon.b=1e-5) {
  kk <- 0
  for (j in mask) {
    kk <- kk+1
    if (j==0)
      for (i in 0:(length(time)-1)) {
        data[kk+i*prod(dim(mask))] <- 0
      }
  }
  data[is.na(data)] <- 0

  require("splines") ; require("MASS") ; require("minpack.lm")

  if (input=="FH"||input=="Weinmann") {
    if (input=="Weinmann"){D=0.1;a1=3.99;a2=4.78;m1=0.244;m2=0.0111;}
    if (input=="FH"){D=1;a1=2.4;a2=0.62;m1=3;m2=0.016;}
    input <- D*(a1*exp(-m1*time)+a2*exp(-m2*time))
  }

  ## function to make precision matrix for random walk
  R <- function(taux,rw) {
    RR <- matrix(0,nrow=length(taux)+rw,ncol=length(taux)+rw)
    if(rw==0) {
      for (i in 1:length(taux)) {
        RR[i,i] <- taux[i]
      }
    }
    if(rw==1) {
      for (i in 1:length(taux)) {
        RR[i,i] <- RR[i,i]+taux[i]
        RR[i+1,i+1] <- RR[i+1,i+1]+taux[i]
        RR[i+1,i] <- RR[i+1,i]-taux[i]
        RR[i,i+1] <- RR[i,i+1]-taux[i]
      }
    }
    if(rw==2) {
      for (i in 1:length(taux)) {
        RR[i,i] <- RR[i,i]+taux[i]
        RR[i+1,i+1] <- RR[i+1,i+1]+4*taux[i]
        RR[i+2,i+2] <- RR[i+2,i+2]+taux[i]
        RR[i+1,i] <- RR[i+1,i]-2*taux[i]
        RR[i,i+1] <- RR[i,i+1]-2*taux[i]
        RR[i+2,i+1] <- RR[i+2,i+1]-2*taux[i]
        RR[i+1,i+2] <- RR[i+1,i+2]-2*taux[i]
        RR[i+2,i] <- RR[i+2,i]+taux[i]
        RR[i,i+2] <- RR[i,i+2]+taux[i]
      }
    }
    return(RR)
  }

  ## redefinitions
  aifzeit <- time.input
  zeit <- time

  ## define B and A
  p <- length(knots)-k
  B <- splineDesign(knots,aifzeit,k)        ## require("splines")
  if (sum(B[,dim(B)[2]]==0)==dim(B)[1])
    B <- B[,-dim(B)[2]]
  if (sum(B[,1]==0)==dim(B)[1])
    B <- B[,-1]
  p <- dim(B)[2]
  A <- matrix(0,nrow=length(zeit),ncol=length(aifzeit))
  ni <- zeit
  for (i in 1:length(zeit))
    for (j in 1:length(aifzeit))
      if (aifzeit[j]<=zeit[i])
        ni[i] <- j
  for (i in 1:length(zeit))
    for (j in 1:ni[i]) {
      A[i,j] <- input[1+ni[i]-j]
    }
  A <- A*mean(diff(aifzeit))
  D <- A%*%B
  
  T <- length(zeit)
  dims <- dim(as.array(data))
  dims <- dims[-length(dims)]
  ddd <- 3
  if (length(dims)==0){
    dims <- c(1,1,1)
    ddd <- 0
  }
  if (length(dims)==1) {
    dims <- c(1,1,dims)
    ddd <- 1
  }
  if (length(dims)==2){
    dims <- c(1,dims)
    ddd <- 2
  }
  data <- array(data,c(dims,T))

  tau  <-  array(1000,c(dims,p-rw,nriters))
  beta  <-  array(0,c(dims,p,nriters))
  MAX  <-  array(0,c(dims,nriters))
  D2  <-  t(D)%*%D
  tauepsilon  <-  array(1000,c(dims,nriters))
  burnin <- min(burnin,nriters)

  print(dims)
  result <- .C("dce_spline_run",
               as.integer(1),
               as.integer(burnin),
               as.integer(c(dims,T)),
               as.double(data),
               as.double(tau),
               as.double(tauepsilon),
               as.double(D),
               as.integer(rw),
               as.double(beta),
               as.double(c(hyper.a,hyper.b,epsilon.a,epsilon.b)),
               as.integer(p),
               as.double(1:T),
               as.double(1:T),
               as.double(1:T),
               as.double(1:(T^2)),
               as.double(1:(T^2)),
               as.double(1:(T^2)),
               as.double(1:(T^2)),
               as.double(t(D)),
               PACKAGE="dcemri")
  tau  <-  array(result[[5]],c(dims,p-rw,nriters))
  beta  <-  array(result[[9]],c(dims,p,nriters))
  tauepsilon  <-  array(result[[6]],c(dims,nriters))
  
  print ("  Burn in done")

  print(tau[,,,,1])
  print(beta[,,,,1])
  print(tauepsilon[,,,1])

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
               as.double(c(hyper.a,hyper.b,epsilon.a,epsilon.b)),
               as.integer(p),
               as.double(1:T),
               as.double(1:T),
               as.double(1:T),
               as.double(1:(T^2)),
               as.double(1:(T^2)),
               as.double(1:(T^2)),
               as.double(1:(T^2)),
               as.double(t(D)),
               PACKAGE="dcemri")
  tau  <-  array(result[[5]],c(dims,p-rw,nriters))
  beta  <-  array(result[[9]],c(dims,p,nriters))
  tauepsilon  <-  array(result[[6]],c(dims,nriters))

  print ("  Iterations done")

  MAX <- apply(beta,c(1:3,5),max)
  MAX <- apply(beta,c(1:3),median)

  parhat <- NULL

  modelTC0 <- function(time, F, E, ve) {
    kep <- E*F/ve
    ######################################################################
    erg  <-  E*exp(-kep*(time-TC)) ## no visible binding for global variable 'TC'
    ######################################################################
    erg <- erg*F
    eval(erg)
  }
  model <- function(time, F, E, TC, ve) {
    TC2 <- 2*TC
    if (TC<0)
      TC2 <- 0
    kep <- E*F/ve
    erg  <-  E*exp(-kep*(time-TC)) # erg[time<TC] <- 1
    erg[time<TC] <- 1-time[time<TC2]*(1-E)/TC
    erg <- erg*F
    if (TC< 0)
      erg <- rep(-10^16,length(time)) # erg <- expression(erg)
    eval(erg)
  }
  fcn <- function(p, time, x, N.Err, fcall, jcall)
    (x - do.call("fcall", c(list(time = time), as.list(p))))
  fcn.jac <- function(p, time, x, N.Err, fcall, jcall) {
    N.Err  <-  rep(N.Err, length(p))
    -do.call("jcall", c(list(time = time), as.list(p)))/N.Err
  }
  
  if (nlr=="med"||nlr=="mean") {
    print ("  Starting non-linear regression")
    X <- dims[1]
    Y <- dims[2]
    Z <- dims[3]

    parhat <- list(Fp=array(0,c(X,Y,Z)),
                   E=array(0,c(X,Y,Z)),
                   ktrans=array(0,c(X,Y,Z)),
                   kep=array(0,c(X,Y,Z)),
                   ve=array(0,c(X,Y,Z)),
                   TC=array(0,c(X,Y,Z)),
                   vp=array(0,c(X,Y,Z)))

    for (i in 1:X)
      for (j in 1:Y)
        for (m in 1:Z) {
          if (nlr=="med")
            me5 <- apply(beta[i,j,m,,],1,median)%*%t(B)
          if (nlr=="mean")
            me5 <- apply(beta[i,j,m,,],1,mean)%*%t(B)
          kt <- max(me5)
          start <- which(kt==me5)
          try(tmp <- nls.lm(par=list("E"=.6,"F"=MAX[i,j,m,],"TC"=0,"ve"=.05), 
                            fn = fcn, fcall = model, #jcall= modelJ,
                            time = zeit, x = me5 ,N.Err=sqrt(300),
                            control = list(nprint=1,ftol=10^-20)))
          
          if (!(parhat$TC>0))
            try(tmp <- nls.lm(par=list("E"=.6,"F"=MAX[i,j,m,],"TC","ve"=.05), 
                              fn = fcn, fcall = modelTC0, #jcall= modelJ,
                              time = zeit, x = me5 ,N.Err=sqrt(300),
                              control = list(nprint=1,ftol=10^-20)))

          parhat$E[i,j,m] <- tmp$par$E
          parhat$Fp[i,j,m] <- tmp$par$F
          parhat$ve[i,j,m] <- tmp$par$ve
          parhat$TC[i,j,m] <- tmp$par$TC
          parhat$ktrans[i,j,m] <- tmp$par$E*tmp$par$F
          parhat$kep[i,j,m] <- parhat$ktrans[i,j,m]/tmp$par$ve
          parhat$vp[i,j,m] <- tmp$par$TC*tmp$par$F
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
    print ("  Starting non-linear regression")
    X <- dims[1]
    Y <- dims[2]
    Z <- dims[3]
    parhat <- list(Fp=array(0,c(X,Y,Z,nriters)),
                   E=array(0,c(X,Y,Z,nriters)),
                   ktrans=array(0,c(X,Y,Z,nriters)),
                   kep=array(0,c(X,Y,Z,nriters)),
                   ve=array(0,c(X,Y,Z,nriters)),
                   TC=array(0,c(X,Y,Z,nriters)),
                   vp=array(0,c(X,Y,Z,nriters)))

    for (i in 1:X)
      for (j in 1:Y)
        for (m in 1:Z) {
          for (k in 1:nriters) {
            me5 <- beta[i,j,m,,k]%*%t(B)
            kt <- max(me5)
            start <- which(kt==me5)
            try(tmp <- nls.lm(par=list("E"=1,"F"=1,"TC"=0.1,"ve"=.1), 
                              fn = fcn, fcall = model, #jcall= modelJ,
                              time = zeit, x = me5 ,N.Err=sqrt(300),
                              control = list(nprint=-1,ftol=10^-20)))

            parhat$E[i,j,m,k] <- tmp$par$E
            parhat$Fp[i,j,m,k] <- tmp$par$F
            parhat$ve[i,j,m,k] <- tmp$par$ve
            parhat$TC[i,j,m,k] <- tmp$par$TC
            parhat$ktrans[i,j,m,k] <- tmp$par$E*tmp$par$F
            parhat$kep[i,j,m,k] <- parhat$ktrans[i,j,m,k]/tmp$par$ve
            parhat$vp[i,j,m,k] <- tmp$par$TC*tmp$par$F
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
    MAX <- MAX[1,1,1]
  }
  if (ddd==2) {
    MAX <- MAX[1,1,]
  }
  if (ddd==3) {
    MAX <- MAX[1,,]
  }

  print ("  Done.")

  ktrans <- apply(parhat$ktrans, 1:ddd, median, na.rm=TRUE)
  kep <- apply(parhat$kep, 1:ddd, median, na.rm=TRUE)
  vp <- apply(parhat$vp, 1:ddd, median, na.rm=TRUE)
  sigma2 <- 1/apply(tauepsilon, 1:3, median, na.rm=TRUE)
  result <- list(fp.hat=MAX, tau=tau, beta=beta, tauepsilon=tauepsilon,
                 A=A, B=B, D=D, knots=knots, estimates=parhat)
  return(result)
}
