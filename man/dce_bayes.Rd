\name{dcemri.bayes}
\alias{dcemri.bayes}
\title{Bayesian Methods for Pharmacokinetic Modeling of Dynamic
  Contrast-Enhanced MRI Data}
\description{
  Bayesian analysis of contrast agent concentration time curves from
  DCE-MRI.
}
\usage{
dcemri.bayes(conc, time, img.mask, model="extended",
                      aif="tofts.kermode", user=NULL, 
	              nriters=9500, thin=30, burnin=2000, tune=267, 
	              tau.ktrans=1, tau.kep=tau.ktrans, ab.vp=c(1,19),
	              ab.tauepsilon=c(1,1/1000), samples=TRUE,
	              ...)
}
\arguments{
  \item{conc}{Matrix or array of concentration time series.
    Last dimension must be time.}
  \item{time}{Time in minutes.}
  \item{img.mask}{Mask matrix or array. Voxels with mask=0 will be excluded.}
  \item{nriters}{Total number of iterations.}
  \item{thin}{Thining factor.}
  \item{burnin}{Number of iterations for burn-in.}
  \item{tune}{Number for iterations for tuning. Algorithm will be 
    tuned to an acceptance rate of 0.3 to 0.6.}
  \item{tau.ktrans}{Variance parameter for prior on $\log(K^{trans})$.}
  \item{tau.kep}{Variance parameter for prior on $\log(k_{ep})$.}
  \item{tau.gamma}{}
  \item{tau.theta}{}
  \item{ab.vp}{Hyper prior parameters for $v_p$ Beta prior.}
  \item{ab.tauepsilon}{Hyper prior parameters for observation error Gamma Prior.}
  \item{aif.model}{Model for arterial input funktion. So far implemented:
    0 or "tofts-kermode" double exponential model by Tofts,Kermode;
    1 or "orton" exponential model by Orton et al.}
  \item{aif.parameter}{Vector of AIF parameters. 
    For Tofts/Kermode: $a_1$,$m_1$,$a_2$,$m_2$; for Orton et al.:
    $A_b$,$\mu_b$,$A_g$,$\mu_g$.}
  \item{vp}{}
  \item{approx}{}
  \item{vp.do}{If TRUE, model includes vascular term.}
  \item{samples}{If TRUE, output includes samples.}
}
\value{
}
\details{
}
\references{
}
\seealso{
}
\examples{

data("buckley")
img <- t(breast$data)[,(1:60)*5]
time <- seq(5,300,5)/60
mask <- rep(TRUE,13)

fit <- dcemri.bayes(img, time, mask,aif="fritz.hansen")
plot(breast$ktrans,fit$ktrans, xlab="True Ktrans",ylab="Estimated Ktrans (Bayes)",xlim=c(0,0.5),ylim=c(0,0.5))
lines(c(-0.1,0.6),c(-0.1,0.6),lty=3, col="grey")

fit.lm <- dcemri.lm(img, time, mask,aif="fritz.hansen")
plot(breast$ktrans, fit$ktrans, xlab="True Ktrans",ylab="Estimated Ktrans",xlim=c(0,0.75),ylim=c(0,0.75))
points(breast$ktrans, fit.lm$ktrans, pch=2)
legend(0.6,0.2,pch=c(2,1),legend=c("LM","Bayes"))
lines(c(-0.1,0.8),c(-0.1,0.8),lty=3, col="grey")

}
\author{Volker Schmid}
\keyword{}
