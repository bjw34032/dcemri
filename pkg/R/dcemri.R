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
## Time-stamp: <>
## $Id$
##

dcemri.lm <- function(conc, time, img.mask, model="extended",
                      aif="tofts.kermode", nprint=0, user=NULL, ...) {
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
  ##        model: AIF... "weinman" or "parker",
  ##        update: re-do given parameter maps where parameters not fitted,
  ##        ktransmap, kepmap, ktranserror, keperror: given parameter maps.
  ##
  ## output: list with ktrans, kep, ve, std.error of ktrans and kep
  ##         (ktranserror and keperror)
  ##

  require("minpack.lm")

  mod <- model
  nvoxels <- sum(img.mask)
  I <- nrow(conc)
  J <- ncol(conc)
  K <- nsli(conc)

  if (!is.numeric(dim(conc))) {I <- J <- K <- 1} 
	else if (length(dim(conc))==2) {J <- K <-1}

  cat("  Deconstructing data...", fill=TRUE)
  conc.mat <- matrix(conc[img.mask], nvoxels)
  conc.mat[is.na(conc.mat)] <- 0

  switch(aif,
         tofts.kermode = {
           D <- 0.1; a1 <- 3.99; a2 <- 4.78; m1 <- 0.144; m2 <- 0.0111
         },
         fritz.hansen = {
           D <- 1; a1 <- 2.4; a2 <- 0.62; m1 <- 3.0; m2 <- 0.016
         },
         orton.exp = {
           D <- 1; AB <- 323; muB <- 20.2; AG <- 1.07; muG <- 0.172
         },
         orton.cos = {
           D <- 1; aB <- 2.84; muB <- 22.8; aG <- 1.36; muG <- 0.171
         },
         user = {
           cat("  User-specified AIF parameters...", fill=TRUE);
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
           func <- function(theta, signal, time)
             signal - model(time, theta[1], theta[2])
           guess <- c("th1"=0, "th3"=0.1)
         },
         extended = {
           model <- model.extended
           func <- function(theta, signal, time)
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

  cat("  Estimating the kinetic parameters...", fill=TRUE)
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

  cat("  Reconstructing results...", fill=TRUE)
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- ktrans$par
  B[img.mask] <- ktrans$error
  ktrans.out <- list(par = A, error = B)
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- kep$par
  B[img.mask] <- kep$error
  kep.out <- list(par = A, error = B)
  if(mod %in% c("extended","orton.exp","orton.cos")) {
    A <- B <- array(NA, c(I,J,K))
    A[img.mask] <- Vp$par
    B[img.mask] <- Vp$error
    Vp.out <- list(par = A, error = B)
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

dcemri.bayes <- function(conc, time, img.mask, model="extended",
                      aif="tofts.kermode", user=NULL, 
	              nriters=9500, thin=30, burnin=2000, tune=267, 
	              tau.ktrans=1, tau.kep=tau.ktrans, ab.vp=c(1,19),
	              ab.tauepsilon=c(1,1/1000), samples=FALSE, multicore=FALSE,
	              ...) {

  ## dcemri.bayes - a function for fitting 1-compartment PK models to
  ## DCE-MRI images using Bayes inference
  ##
  ## authors: Volker Schmid, Brandon Whitcher
  ##
  ## input:
  ##        conc: array of Gd concentration,
  ##        time: timepoints of aquisition,
  ##        img.mask: array of voxels to fit,
  ##        D(=0.1): Gd dose in mmol/kg,
  ##        model: AIF... "weinman" or "parker",
  ##
  ## output: list with ktrans, kep, ve, std.error of ktrans and kep
  ##         (ktranserror and keperror), samples if samples=TRUE
  ##


dce.bayes.single<-function(conc,time,nriters=7000,thin=10,burnin=2000,tune=267,tau.gamma=1,tau.theta=1,ab.vp=c(1,19),ab.tauepsilon=c(1,1/1000),
aif.model=0,aif.parameter=c(2.4,0.62,3,0.016),vp=1)
  {
  if (sum(is.na(conc))>0)return(NA)
  else
    {
  
 n<-floor((nriters-burnin)/thin)
if (tune>(0.5*nriters))tune=floor(nriters/2);

  singlerun<-.C("dce_bayes_run_single",as.integer(c(nriters,thin,burnin,tune)),as.double(conc),
    as.double(tau.gamma),as.double(tau.theta),as.double(ab.vp),as.double(ab.tauepsilon),
    as.double(c(aif.model,aif.parameter)),as.integer(vp),as.double(time),as.integer(length(time)),
    as.double(rep(0,n)),as.double(rep(0,n)),as.double(rep(0,n)),as.double(rep(0,n),as.integer(0)), PACKAGE="dcemri")
    
    return(list("ktrans"=singlerun[[11]],"kep"=singlerun[[12]],"vp"=singlerun[[13]],"sigma2"=1/singlerun[[14]]))
    }
}

  mod <- model
  nvoxels <- sum(img.mask)
  I <- nrow(conc)
  J <- ncol(conc)
  K <- nsli(conc)

  if (!is.numeric(dim(conc))) {I <- J <- K <- 1} 
	else if (length(dim(conc))==2) {J <- K <-1}

  cat("  Deconstructing data...", fill=TRUE)
  conc.mat <- matrix(conc[img.mask], nvoxels)
  conc.mat[is.na(conc.mat)] <- 0

  switch(aif,
         tofts.kermode = {
           D <- 0.1; a1 <- 3.99; a2 <- 4.78; m1 <- 0.144; m2 <- 0.0111
         },
         fritz.hansen = {
           D <- 1; a1 <- 2.4; a2 <- 0.62; m1 <- 3.0; m2 <- 0.016
         },
         orton.exp = {
           D <- 1; AB <- 323; muB <- 20.2; AG <- 1.07; muG <- 0.172
         },
         orton.cos = {
           D <- 1; aB <- 2.84; muB <- 22.8; aG <- 1.36; muG <- 0.171
         },
         user = {
           cat("  User-specified AIF parameters...", fill=TRUE);
           D <- try(user$D); AB <- try(user$AB); aB <- try(user$aB);
           muB <- try(user$muB); AG <- try(user$AG); aG <- try(user$aG); 
           muG <- try(user$muG)
         },
         print("WARNING: AIF parameters must be specified!"))
  
	# translate "model" to "aif.model" and "vp.do"
	switch(model,
	weinmann={aif.model=0
	vp.do=FALSE},
	extended={aif.model=0
	vp.do=TRUE},
	orton.exp={aif.model=1
	vp.do=TRUE},
	stop("Model is not supported."))


  ktrans <- kep <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  sigma2 <- rep(NA, nvoxels)
  Vp <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))

  if (mod %in% c("extended","weinmann")){aif.parameter=c(D*a1,D*a2,m1,m2)}
	else{ aif.parameter=c(D*aB,D*aG,muB,muG)}

  if (samples)
	{
	sigma2.samples <- ktrans.samples <- kep.samples <- c()
        if(mod %in% c("extended","orton.exp","orton.cos")) {
	     Vp.samples <- c()
             }
	}

  cat("  Estimating the kinetic parameters...", fill=TRUE)


   conc.list<-list()
   for (i in 1:nvoxels)
    conc.list[[i]]=conc.mat[i,]

  if (!multicore)
  {
   fit <- lapply(conc.list,dce.bayes.single,time,nriters,thin,burnin,tune,tau.ktrans,
   tau.kep,ab.vp,ab.tauepsilon,aif.model,aif.parameter,vp.do) 
  }
  else
  {
   require(multicore)
   fit <- mclapply(conc.list,dce.bayes.single,time,nriters,thin,burnin,tune,tau.ktrans,
   tau.kep,ab.vp,ab.tauepsilon,aif.model,aif.parameter,vp.do) 
  }

  cat("  Reconstructing results...", fill=TRUE)

   for(k in 1:nvoxels) {
      ktrans$par[k] <- median(fit[[k]]$ktrans)
      kep$par[k] <- median(fit[[k]]$kep)
      ktrans$error[k] <- sqrt(var(fit[[k]]$ktrans))
      kep$error[k] <- sqrt(var(fit[[k]]$kep))
      if(mod %in% c("extended","orton.exp","orton.cos")) {
        Vp$par[k] <- median(fit[[k]]$vp)
        Vp$error[k] <- sqrt(var(fit[[k]]$vp))
	if (samples){Vp.samples<-c(Vp.samples,fit[[k]]$v)}
       }
     sigma2[k] <- median(fit[[k]]$sigma2)
     if (samples)
	{
	   ktrans.samples<-c(ktrans.samples,fit[[k]]$ktrans)
	   kep.samples<-c(kep.samples,fit[[k]]$kep)
	   sigma2.samples<-c(sigma2.samples,fit[[k]]$sigma2)
 	}
      }
   
  

  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- ktrans$par
  B[img.mask] <- ktrans$error
  ktrans.out <- list(par = A, error = B)
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- kep$par
  B[img.mask] <- kep$error
  kep.out <- list(par = A, error = B)
  if(mod %in% c("extended","orton.exp","orton.cos")) {
    A <- B <- array(NA, c(I,J,K))
    A[img.mask] <- Vp$par
    B[img.mask] <- Vp$error
    Vp.out <- list(par = A, error = B)
  }
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- sigma2
  sigma2.out <- A
  

  if (samples)
	{
	  extract.samples <- function(sample,I,J,K,NRI)
		{
		A <- array(NA, c(I,J,K,NRI))
		for (i in 1:I)
		for (j in 1:J)
		for (k in 1:K)
		{
		voxelcount = i+(j-1)*I+(k-1)*I*J-1
		A[i,j,k,] <- sample[(1:NRI) + voxelcount*NRI]
		}
		return(A)
		}

            NRI <- length(ktrans.samples)/length(ktrans$par)
     	    ktrans.out <- list(par = ktrans.out$par, error=ktrans.out$error, samples = extract.samples(ktrans.samples,I,J,K,NRI))
     	    kep.out <- list(par = kep.out$par, error=kep.out$error, samples = extract.samples(kep.samples,I,J,K,NRI))
 	  if(mod %in% c("extended","orton.exp","orton.cos")) {
             	Vp.out <- list(par = Vp.out$par, error=Vp.out$error, samples = extract.samples(Vp.samples,I,J,K,NRI))
	   }
             sigma2.samples <- extract.samples(sigma2.samples,I,J,K,NRI)
	}

  if(mod %in% c("extended","orton.exp","orton.cos"))
     {
    if (samples)
	{
	list(ktrans=ktrans.out$par, kep=kep.out$par, ktranserror=ktrans.out$error,
         keperror=kep.out$error, ve=ktrans.out$par/kep.out$par, vp=Vp.out$par,
         vperror=Vp.out$error, sigma2=sigma2.out, ktrans.samples=ktrans.out$samples, 
         kep.samples = kep.out$samples, vp.samples = Vp.out$samples,
	 sigma2.samples = sigma2.samples, time=time)
	}
	else
	{
	list(ktrans=ktrans.out$par, kep=kep.out$par, ktranserror=ktrans.out$error,
         keperror=kep.out$error, ve=ktrans.out$par/kep.out$par, vp=Vp.out$par,
         vperror=Vp.out$error, sigma2=sigma2.out, time=time)
	}
      }
  else
	{
	if (samples)
	{
	    list(ktrans=ktrans.out$par, kep=kep.out$par, ktranserror=ktrans.out$error,
        	 keperror=kep.out$error, ve=ktrans.out$par/kep.out$par, sigma2=sigma2.out,
         	 ktrans.samples=ktrans.out$samples, kep.samples = kep.out$samples, 
		 sigma2.samples = sigma2.samples, time=time)
	}
	else
	{
	    list(ktrans=ktrans.out$par, kep=kep.out$par, ktranserror=ktrans.out$error,
        	 keperror=kep.out$error, ve=ktrans.out$par/kep.out$par, sigma2=sigma2.out,
         	time=time)
	}
	}
}


dcemri.spline <- function(conc, time, img.mask, time.input=time, model="weinmann", aif="tofts.kermode", user=NULL, aif.observed=NULL,
	       nriters=500, thin=5, burnin=100, 
	       ab.hyper=c(1e-5,1e-5), ab.tauepsilon=c(1,1/1000), 
	       k=4,p=25,rw=2,
               knots=seq(min(time)-1e-3-(k-1)*(max(time)-min(time))/(knotpoints-k+1),1e-3+max(time)+k*(max(time)-min(time))/(knotpoints-k+1),(2e-3+max(time)-min(time))/(knotpoints-k+1)), 
	       nlr = FALSE, t0.compute=FALSE, samples=FALSE, multicore=FALSE,
                       
	              ...) {

  ## dcemri.spline - a function for fitting Bayesian Penalty Splines to 
  ## DCE-MRI images and computing kinetic parameters
  ##
  ## authors: Volker Schmid, Brandon Whitcher
  ##
  ## input:
  ##        conc: array of Gd concentration,
  ##        time: timepoints of aquisition,
  ##        img.mask: array of voxels to fit,
  ##        D(=0.1): Gd dose in mmol/kg,
  ##        model: AIF... "weinman" or "parker",
  ##
  ## output: list with ktrans, kep, ve, std.error of ktrans and kep
  ##         (ktranserror and keperror)
  ##

##internal function
    q05 <- function(x)
      quantile(x, .005, na.rm=TRUE)
    med.na <- function(x)
      median(x, na.rm=TRUE)

  fcn <- function(p, time, x, N.Err, fcall, jcall)
    (x - do.call("fcall", c(list(time = time), as.list(p))))

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

nls.lm.single <- function(fitted, par,
                         fn, fcall,model,time)
	{
	fcall2=fcall
	if (length(fcall)>1) fcall2=fcall[[1]]
	fit=nls.lm(par=par, fn=fn, fcall=fcall2, time=time, x=fitted, N.Err=sqrt(300),
                            control=list(nprint=0,ftol=10^-20))
	if (model=="AATH")
	if(fit$par$TC<1e-6)
	{
	fit=nls.lm(par=par[-3], fn=fn, fcall=fcall[[2]], time=time, x=fitted, N.Err=sqrt(300),
                            control=list(nprint=0,ftol=10^-20))
	fit$par$TC=0
	}
	return(fit)
}

dce.spline.single<-function(conc, time, D, inputtime, p, rw, knots, t0.compute=FALSE, nlr=FALSE, nriters=500,thin=5, burnin=100,ab.hyper=c(1e-5,1e-5),ab.tauepsilon=c(1,1/1000),silent=0, multicore=FALSE,model=NULL,model.func=NULL,model.guess=NULL,samples=FALSE)
  {
  if (sum(is.na(conc))>0)return(NA)
  else
    {
      T=length(time)

  tau <- array(1000, c(p-rw,nriters))
  beta <- array(0, c(p,nriters))
  MAX <- array(0, c(nriters))
  tauepsilon <- array(1000, c(nriters))
  burnin <- min(burnin, nriters)

  result <- .C("dce_spline_run",
               as.integer(1),
               as.integer(burnin),
               as.integer(c(1,1,1,T)),
               as.double(conc),
               as.double(tau),
               as.double(tauepsilon),
               as.double(D),
               as.integer(rw),
               as.double(beta),
               as.double(c(ab.hyper[1], ab.hyper[2], ab.tauepsilon[1], ab.tauepsilon[2])),
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
	}

  tau <- array(result[[5]], c(p-rw,nriters))
  beta <- array(result[[9]], c(p,nriters))
  tauepsilon <- array(result[[6]], c(nriters))

  result <- .C("dce_spline_run",
               as.integer(nriters),
               as.integer(thin),
               as.integer(c(1,1,1,T)),
               as.double(conc),
               as.double(tau),
               as.double(tauepsilon),
               as.double(D),
               as.integer(rw),
               as.double(beta),
               as.double(c(ab.hyper[1], ab.hyper[2], ab.tauepsilon[1], ab.tauepsilon[2])),
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
	
  tau <- array(result[[5]], c(p-rw,nriters))
  beta <- array(result[[9]], c(p,nriters))
  tauepsilon <- array(result[[6]], c(nriters))


  t0=time[1]
  if(t0.compute)
        {
	    d <- array(NA, c(T, nriters))
    for (j in 1:nriters)
        d[,j] <- D %*% beta[,j]
    d1 <- apply(d, 1, q05)
    d2 <- apply(d, 1, med.na)

    du <- min(which(d1 > 0))
                
    beta.abl <- beta.abl2 <- rep(0, p)
    
    B2 <- splineDesign(knots, inputtime, k-2)
    B2 <- B2[,(1:p)+1]
    
    for (j in 1:nriters) {
            beta.abl <- 1:p
            for (q in 1:(p-1))
              beta.abl[q] <- (beta[q+1,j] - beta[q,j]) * k / (knots[q+k+1] - knots[q+1])
            beta.abl[p] <- 0
            ABL2 <- A %*% B2 %*% beta.abl
            du2 <- time[du] - d2[du] / ABL2[du]  
            t0[j] <- du2
    	}
  	if (t0<0) t0<-0
  	if (t0>max(time)) t0<-0
    }

	fitted=list()
	for (i in 1:nriters)
	fitted[[i]] <- B%*% beta[,i]
	if (multicore)
		MAX0= mclapply(fitted,max)
	else
		MAX0= lapply(fitted,max)
	MAX = c()
	for (i in 1:nriters)
	   MAX=c(MAX,MAX0[[i]])

parameters=list()

if (nlr)
	{
	if(model=="AATH") model.guess[2]=median(MAX)
	if (multicore)
	response <- mclapply(fitted,nls.lm.single, par=model.guess,
                         fn=fcn, fcall = model.func, model=model,
                         time=time-t0)
	else
	response <- lapply(fitted,nls.lm.single, par=model.guess,
                         fn=fcn, fcall = model.func, model=model,
                         time=time-t0)
	
	if (model=="AATH")
	{
	E=F=TC=ve=c()
	for (i in 1:nriters)
	{
		E=c(E,response[[i]]$par$E)
		F=c(F,response[[i]]$par$F)
		TC=c(TC,response[[i]]$par$TC)
		ve=c(ve,response[[i]]$par$ve)
	}
	parameters=list("E"=median(E),"F"=median(F),"TC"=median(TC),"ve"=median(ve))
	if(samples)parameters=list("E.samples"=E,"F.samples"=F,"TC.samples"=TC,"ve.samples"=ve)
	}
	if (model=="weinmann")
	{
	ktrans=kep=c()
	for (i in 1:nriters)
	{
		ktrans=c(ktrans,response[[i]]$par$logktrans)
		kep=c(kep,response[[i]]$par$logkep)
	}
	parameters=list("ktrans"=median(exp(ktrans)),"kep"=median(exp(kep)))
	if (samples) parameters=list("ktrans.sample"=exp(ktrans),"kep.sample"=exp(kep))
	}

	}

  return(list("beta"=beta,"tau"=tau,"tauepsilon"=tauepsilon,"t0"=t0,"Fp"=median(MAX),"Fp.samples"=MAX,"fitted"=fitted,"par"=parameters))
}

# main function

  knotpoints=p
  mod <- model
  nvoxels <- sum(img.mask)
  I <- nrow(conc)
  J <- ncol(conc)
  K <- nsli(conc)
  T <- length(time)

  if (!is.numeric(dim(conc))) {I <- J <- K <- 1} 
	else if (length(dim(conc))==2) {J <- K <-1}

  cat("  Deconstructing data...", fill=TRUE)
  conc.mat <- matrix(conc[img.mask], nvoxels)
  conc.mat[is.na(conc.mat)] <- 0
   conc.list<-list()
   for (i in 1:nvoxels)
    conc.list[[i]]=conc.mat[i,]

  switch(aif,
         tofts.kermode = {
           D <- 0.1; a1 <- 3.99; a2 <- 4.78; m1 <- 0.144; m2 <- 0.0111
           input <- D*(a1*exp(-m1*time)+a2*exp(-m2*time))
         },
         fritz.hansen = {
           D <- 1; a1 <- 2.4; a2 <- 0.62; m1 <- 3.0; m2 <- 0.016
           input <- D*(a1*exp(-m1*time)+a2*exp(-m2*time))	 
	},
	observed = {
	   input = aif.observed
       },
         print("WARNING: AIF parameters must be specified!"))
  
  model.func <- model.guess <- NULL

  if (model=="weinmann")
	{
	  ktrans <- kep <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  	  sigma2 <- rep(NA, nvoxels)
	model.func <- function(time, logktrans, logkep) {
             ktrans=exp(logktrans)
	     kep=exp(logkep)
             erg <- ktrans*exp(-kep*(time))
             eval(erg)
                }
 
	model.guess <- list("logktrans"=-1,"logkep"=0)
	}

  if (model=="AATH")
	{
	  E <- F <- TC <- ve <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  	  sigma2 <- rep(NA, nvoxels)

  	model.func <- list()
	model.func[[1]] <- function(time, F, E, TC, ve) {
    		TC2<-2*TC
    		if (TC < 0) { TC2 <- 0 }
    		kep <- E*F/ve
    		erg <- E*exp(-kep*(time-TC))
     		erg[time<TC] <- 1 - time[time<TC2]*(1-E) / TC
    		erg <- erg*F
    		if (TC < 0)
      			erg <- rep(-10^16, length(time))
	        eval(erg)
  		}
	model.func[[2]] <- function(time, F, E, ve) {
             kep <- E*F/ve
             erg <- E*exp(-kep*(time))
             erg <- erg*F
             eval(erg)
                }
 
	model.guess <- list("E"=.6,"F"=2,"TC"=0,"ve"=.05)
	}


  ##define B and A
  p <- length(knots)-k
  B <- splineDesign(knots, time.input, k, outer.ok=TRUE)
  if (sum(B[,dim(B)[2]]==0) == dim(B)[1])
    B <- B[,-dim(B)[2]]
  if (sum(B[,1]==0) == dim(B)[1])
    B <- B[,-1]
  p <- dim(B)[2]
  A <- matrix(0, nrow=length(time), ncol=length(time.input))
  ni <- time
  for (i in 1:length(time)) {
    for (j in 1:length(time.input)) {
      if (time.input[j]<=time[i])
        ni[i] <- j
    }
  }
  for (i in 1:length(time)) {
    for (j in 1:ni[i])
      A[i,j] <- input[1+ni[i]-j]
  }
  A <- A*mean(diff(time.input))
  D <- A%*%B
  T <- length(time)



  cat("  Estimating the parameters...", fill=TRUE)

  if (!multicore)
  {
   fit <- lapply(conc.list,dce.spline.single, time, D, time.input, p, rw, knots, nriters=500,thin=5, burnin=100,ab.hyper=c(1e-5,1e-5),ab.tauepsilon=c(1,1/1000),t0.compute=t0.compute,nlr=nlr,multicore=multicore,model=model,model.func=model.func,model.guess=model.guess,samples=samples)
  }
  else
  {
   require(multicore)
   fit <- mclapply(conc.list,dce.spline.single, time, D, time.input, p, rw, knots, nriters=500,thin=5, burnin=100,ab.hyper=c(1e-5,1e-5),ab.tauepsilon=c(1,1/1000),t0.compute=t0.compute,nlr=nlr,multicore=multicore,model=model,model.func=model.func,model.guess=model.guess,samples=samples)
  }

  cat("  Reconstructing results...", fill=TRUE)

  t0 <- c()
  for (k in 1:nvoxels)t0<-c(t0,fit[[k]]$t0)
  t0.img <- array(NA,c(I,J,K))
  t0.img[img.mask]=t0
  t0=t0.img

  Fp <- c()
  for (k in 1:nvoxels)Fp<-c(Fp,fit[[k]]$Fp)
  Fp.img <- array(NA,c(I,J,K))
  Fp.img[img.mask]=Fp
  Fp.samples <- array(NA,c(nvoxels,nriters))
  for (i in 1:nvoxels)Fp.samples[i,]=fit[[k]]$Fp.samples
  if(samples)Fp=array(NA,c(I,J,K,nriters))
  if(samples)for (j in 1:nriters)Fp[,,,j][img.mask]=Fp.samples[,j]

	
  if (nlr)
  {
    ktrans<-ve<-c()
    if (model=="weinmann")
    {
	kep=c()
    	if (samples)ktrans.sample <- array(NA,c(nvoxels,nriters))
    	for (k in 1:nvoxels)
    	{
    	ktrans<-c(ktrans,fit[[k]]$par$ktrans)
    	if (samples)ktrans[k,]=fit[[k]]$par$ktrans.samples
    	}
    	if (samples)kep.sample <- array(NA,c(nvoxels,nriters))
    	for (k in 1:nvoxels)
    	{
    	kep<-c(kep,fit[[k]]$par$kep)
    	if (samples)kep[k,]=fit[[k]]$par$kep.samples
    	}
	ve = ktrans/kep
	if (samples) ve.samples=ktrans.samples/kep.samples
    }
    if (model=="AATH")
    {
	E<-F<-TC<-c()
   	if (samples)E.sample <- array(NA,c(nvoxels,nriters))
    	for (k in 1:nvoxels)
    	{
    	E<-c(E,fit[[k]]$par$E)
    	if (samples)E[k,]=fit[[k]]$par$E.samples
    	}
   	if (samples)F.sample <- array(NA,c(nvoxels,nriters))
    	for (k in 1:nvoxels)
    	{
    	F<-c(F,fit[[k]]$par$F)
    	if (samples)F[k,]=fit[[k]]$par$F.samples
    	}
   	if (samples)TC.sample <- array(NA,c(nvoxels,nriters))
    	for (k in 1:nvoxels)
    	{
    	TC<-c(TC,fit[[k]]$par$TC)
    	if (samples)TC[k,]=fit[[k]]$par$TC.samples
    	}
   	if (samples)ve.sample <- array(NA,c(nvoxels,nriters))
    	for (k in 1:nvoxels)
    	{
    	ve<-c(ve,fit[[k]]$par$ve)
    	if (samples)ve[k,]=fit[[k]]$par$ve.samples
    	}
	ktrans = E*F
	if (samples) ktrans.samples=E.samples/F.samples
    }
  }


  beta.sample=array(NA,c(nvoxels,p,nriters))
  for(k in 1:nvoxels) {
		beta.sample[k,,]<-fit[[k]]$beta
 	}
  response.sample=array(NA,c(nvoxels,T,nriters))
  for(k in 1:nvoxels) 
  for(j in 1:nriters) {
		response.sample[k,,j]<-fit[[k]]$fitted[[j]]
 	}

  fitted.sample = array(NA,c(nvoxels,T,nriters))
  for (i in 1:nvoxels)
  for (j in 1:nriters)
	fitted.sample[i,,j]=D%*%beta.sample[i,,j]

  beta <- array(NA,c(I,J,K,p,nriters))
  beta.med <- array(NA,c(I,J,K,p))
  fitted <- array(NA,c(I,J,K,T,nriters))
  fitted.med <- array(NA,c(I,J,K,T))
  response <- array(NA,c(I,J,K,T,nriters))
  response.med <- array(NA,c(I,J,K,T))

  if (nlr)
	{
         ktrans.med <- ve.med <- array(NA,c(I,J,K))
         ktrans.med[img.mask]=ktrans
	 ve.med[img.mask]=ve
	if (model=="weinmann")
	{
        kep.med <- array(NA,c(I,J,K))
	kep.med[img.mask]=kep
	}
	if (model=="AATH")
	{
        E.med <- F.med <- TC.med <- array(NA,c(I,J,K))
	E.med[img.mask]=E
	F.med[img.mask]=F
	TC.med[img.mask]=TC
	} 
	if (samples)
	  {
	  ktrans <- ve <- array(NA,c(I,J,K,nriters))
  	  for (i in 1:nriters)
		{
		  ktrans[,,,i][img.mask]<-ktrans.sample[,i]
		  ve[,,,i][img.mask]<-ve.sample[,i]
		}
	   if (model=="weinmann")
		{
		  kep <- array(NA,c(I,J,K,nriters))
	  	  for (i in 1:nriters)
			{
			  kep[,,,i][img.mask]<-kep.sample[,i]
			}
		}
	   if (model=="AATH")
		{
		  E <- F <- TC <- array(NA,c(I,J,K,nriters))
	  	  for (i in 1:nriters)
			{
			  E[,,,i][img.mask]<-E.sample[,i]
			  F[,,,i][img.mask]<-F.sample[,i]
			  TC[,,,i][img.mask]<-TC.sample[,i]
			}
		}
          }
       }
	
  for (j in 1:p)
   {
    beta.med[,,,j][img.mask]<-apply(beta.sample[,j,],1,median)
	for (i in 1:nriters)
	beta[,,,j,i][img.mask]<-beta.sample[,j,i]
   }

  for (j in 1:T)
   {
      response.med[,,,j][img.mask]<-apply(response.sample[,j,],1,median)
	for (i in 1:nriters)
	response[,,,j,i][img.mask]<-response.sample[,j,i]
     fitted.med[,,,j][img.mask]<-apply(fitted.sample[,j,],1,median)
	for (i in 1:nriters)
	fitted[,,,j,i][img.mask]<-fitted.sample[,j,i]
   }

return.list<-list("beta"=beta.med)
return.list <- c(return.list,list("beta.sample"=beta))
return.list <- c(return.list,list("beta.test"=beta.sample))
return.list <- c(return.list,list("fit"=fitted.med))
return.list <- c(return.list,list("fit.sample"=fitted))
return.list <- c(return.list,list("response"=response.med))
return.list <- c(return.list,list("response.sample"=response))
return.list <- c(return.list,list("Fp"=Fp.img))
if(samples)return.list <- c(return.list,list("Fp.samples"=Fp))
return.list <- c(return.list,list("A"=A))
return.list <- c(return.list,list("B"=B))
return.list <- c(return.list,list("D"=D))
if (nlr)return.list <- c(return.list,list("ktrans"=ktrans.med))
if (nlr)return.list <- c(return.list,list("ve"=ve.med))
if (nlr&(model=="weinmann"))return.list <- c(return.list,list("kep"=kep.med))
if (nlr&(model=="AATH"))return.list <- c(return.list,list("E"=E.med))
if (nlr&(model=="AATH"))return.list <- c(return.list,list("F"=F.med))
if (nlr&(model=="AATH"))return.list <- c(return.list,list("TC"=TC.med))
if (nlr&samples)return.list <- c(return.list,list("ktrans.sample"=ktrans))
if (nlr&samples)return.list <- c(return.list,list("ve.sample"=ve))
if (nlr&samples&(model=="weinmann"))return.list <- c(return.list,list("kep.samples"=kep))
if (nlr&samples&(model=="AATH"))return.list <- c(return.list,list("E.samples"=E))
if (nlr&samples&(model=="AATH"))return.list <- c(return.list,list("F.samples"=F))
if (nlr&samples&(model=="AATH"))return.list <- c(return.list,list("TC.samples"=TC))
return.list <- c(return.list,list("t0"=t0))

return(return.list)
}
