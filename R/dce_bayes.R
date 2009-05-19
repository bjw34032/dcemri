dce.bayes.single <- function(conc, time, nriters=9500, thin=30,
                             burnin=2000, tune=267, tau.gamma=1, tau.theta=1,
                             ab.vp=c(1,19), ab.tauepsilon=c(1,1/1000),
                             aif.model=0, aif.parameter=c(2.4,0.62,3,0.016),
                             vp=1) {
  if (sum(is.na(conc)) > 0)
    return(NA)
  else {
    n <- floor((nriters-burnin)/thin)
    if (tune > (0.5*nriters))
      tune <- floor(nriters/2)
    test <- .C("dce_bayes_run_single",
               as.integer(c(nriters, thin, burnin, tune)),
               as.double(conc),
               as.double(tau.gamma),
               as.double(tau.theta),
               as.double(ab.vp),
               as.double(ab.tauepsilon),
               as.double(c(aif.model, aif.parameter)),
               as.integer(vp),
               as.double(time),
               as.integer(length(time)),
               as.double(rep(0,n)),
               as.double(rep(0,n)),
               as.double(rep(0,n)),
               as.double(rep(0,n)),
               PACKAGE="dcemri")
    return(list("ktrans"=test[[11]], "kep"=test[[12]], "vp"=test[[13]],
                "sigma2"=1/test[[14]]))
  }
}

dce.bayes.wrapper <- function(conc, time, nriters=9500, thin=30,
                              burnin=2000, tune=267, ab.vp=c(1,19),
                              ab.tauepsilon=c(1,1/1000), aif.model=0,
                              aif.parameter=c(2.4,0.62,3,0.016), vp=1) {
  if (aif.model=="tofts-kermode")
    aif.model <- 0
  if (aif.model=="orton")
    aif.model <- 1
  D <- dim(conc)
  D <- c(D[1:(length(D)-1)], rep(1,4-length(D)))
  conc[is.na(conc)] <- 0
  n <- floor((nriters-burnin)/thin)
  N <- prod(D)
  test <- .C("dce_bayes_run",
             as.integer(c(nriters,thin,burnin,tune)),
             as.double(conc),
             as.double(rep(0,length(time))),
             as.double(1),
             as.double(1),
             as.double(ab.vp),
             as.double(ab.tauepsilon),
             as.double(c(aif.model,aif.parameter)),
             as.integer(vp),
             as.double(time),
             as.integer(D),
             as.integer(length(time)),
             as.double(rep(0,n*N)),
             as.double(rep(0,n*N)),
             as.double(rep(0,n*N)),
             as.double(rep(0,n*N)),
             PACKAGE="dcemri")
  return(list("ktrans"=test[[13]], "kep"=test[[14]], "vp"=test[[15]],
              "sigma2"=1/test[[16]]))
}

dce.bayes.single.temp <- function(conc, time, nriters, thin, burnin, tune,
                                  tau.gamma, tau.theta, ab.vp, ab.tauepsilon,
                                  aif.model, aif.parameter, vp) {
  temp <- dce.bayes.single(conc, time, nriters, thin, burnin, tune, tau.gamma,
                           tau.theta, ab.vp, ab.tauepsilon, aif.model,
                           aif.parameter, vp)
  if (is.na(temp[1])) {
    #######################################################################
    temp <- rep(NA, 4*n) ## no visible binding for global variable 'n'
    #######################################################################
    print("Fit fail")
  } else {
    temp <- c(temp$ktrans, temp$kep, temp$vp, temp$sigma2)
  }
  return(temp)
}

dce.bayes <- function(conc, time,
                      mask=array(1,dim(conc)[1:(length(dim(conc))-1)]),
                      nriters=9500, thin=30, burnin=2000, tune=267,
                      tau.ktrans=1, tau.kep=tau.ktrans, ab.vp=c(1,19),
                      ab.tauepsilon=c(1,1/1000), aif.model=0,
                      aif.parameter=c(2.4,3,0.62,0.016), vp.do=TRUE,
                      samples=TRUE) {
  if (aif.model=="tofts-kermode")
    aif.model <- 0
  if (aif.model=="orton")
    aif.model <- 1
  
  d <- length(dim(conc))
  n <- floor((nriters-burnin)/thin)
  
  D <- dim(conc)
  T <- D[length(D)]
  D0 <- D <- D[-length(D)]
  D <- c(D, rep(1,5-length(D)))
  conc <- array(conc,c(D,T))
  
  temp <- conc[,,,,,1]
  temp <- temp[mask]
  data <- array(NA,c(length(temp),T))
  data[,1] <- temp
  
  for (i in 2:T) {
    temp <- conc[,,,,,i]
    temp <- temp[mask]
    data[,i] <- temp
  }
  
  t0 <- max(which(time<0))
  time <- time[-(1:t0)]
  data <- data[,-(1:t0)]
  
  N <- dim(data)[1]
  
  temp <- apply(data, 1, dce.bayes.single.temp, time, nriters, thin, burnin,
                tune, tau.ktrans, tau.kep, ab.vp, ab.tauepsilon, aif.model,
                aif.parameter, vp.do)
  
  D1 <- D2 <- D3 <- D4 <- D5 <- array(1,D)
  for (i in 1:D[1])
    D1[i,,,,] <- i
  for (i in 1:D[2])
    D2[,i,,,] <- i
  for (i in 1:D[3])
    D3[,,i,,] <- i
  for (i in 1:D[4])
    D4[,,,i,] <- i
  for (i in 1:D[5])
    D5[,,,,i] <- i
  
  mask <- array(mask,D)
  
  D1 <- D1[mask]
  D2 <- D2[mask]
  D3 <- D3[mask]
  D4 <- D4[mask]
  D5 <- D5[mask]
  
  coord <- cbind(D1,D2,D3,D4,D5)
    
  ktrans <- temp[1:n,]
  kep <- temp[n+(1:n),]
  if (vp.do)
    vp <- temp[2*n+(1:n),]
  sigma2 <- temp[3*n+(1:n),]
  
  ktrans.mx <- kep.mx <- vp.mx <- sigma2.mx <- array(NA,c(D,n))
  for (i in 1:N) {
    ktrans.mx[D1[i], D2[i], D3[i], D4[i], D5[i],] <- ktrans[,i]
    kep.mx[D1[i], D2[i], D3[i], D4[i], D5[i],] <- kep[,i]
    if (vp.do)vp.mx[D1[i], D2[i], D3[i], D4[i], D5[i],] <- vp[,i]
    sigma2.mx[D1[i], D2[i], D3[i], D4[i], D5[i],] <- sigma2[,i]
  }
  
  ktrans.med <- apply(ktrans.mx, 1:5, median)
  kep.med <- apply(kep.mx, 1:5, median)
  if (vp.do)vp.med <- apply(vp.mx, 1:5, median)
  sigma2.med <- apply(sigma2.mx, 1:5, median)
  
  ktrans <- array(ktrans.mx, c(D0,n))
  kep <- array(kep.mx, c(D0,n))
  if (vp.do)vp <- array(vp.mx, c(D0,n))
  sigma2 <- array(sigma2.mx, c(D0,n))
  
  ktrans.med <- array(ktrans.med, D0)
  kep.med <- array(kep.med, D0)
  if (vp.do)vp.med <- array(vp.med, D0)
  sigma2.med <- array(sigma2.med, D0)

  print("  done")
  if (vp.do && samples)
    res <- list("ktrans"=ktrans.med, "kep"=kep.med, "vp"=vp.med,
                "sigma2"=sigma2.med, "ktrans.sample"=ktrans, "kep.sample"=kep,
                "vp.sample"=vp, "sigma2.sample"=sigma2)
  if (!vp.do && samples)
    res <- list("ktrans"=ktrans.med, "kep"=kep.med, "sigma2"=sigma2.med,
                "ktrans.sample"=ktrans, "kep.sample"=kep,
                "sigma2.sample"=sigma2)
  if (vp.do && !samples)
    res <- list("ktrans"=ktrans.med, "kep"=kep.med, "vp"=vp.med,
                "sigma2"=sigma2.med)
  if (!vp.do && !samples)
    res <- list("ktrans"=ktrans.med, "kep"=kep.med, "sigma2"=sigma2.med)
  
  return(res)
}


