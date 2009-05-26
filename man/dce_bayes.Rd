\name{dce.bayes}
\alias{dce.bayes.single}
\alias{dce.bayes.wrapper}
\alias{dce.bayes.single.temp}
\alias{dce.bayes}
\title{Bayesian Methods for Pharmacokinetic Modeling of Dynamic
  Contrast-Enhanced MRI Data}
\description{
  Bayesian analysis of contrast agent concentration time curves from
  DCE-MRI.
}
\usage{
dce.bayes.single(conc, time, nriters=9500, thin=30,
                             burnin=2000, tune=267, tau.gamma=1, tau.theta=1,
                             ab.vp=c(1,19), ab.tauepsilon=c(1,1/1000),
                             aif.model=0, aif.parameter=c(2.4,0.62,3,0.016),
                             vp=1, approx=FALSE)
dce.bayes.wrapper(conc, time, nriters=9500, thin=30,
                              burnin=2000, tune=267, ab.vp=c(1,19),
                              ab.tauepsilon=c(1,1/1000), aif.model=0,
                              aif.parameter=c(2.4,0.62,3,0.016), vp=1)
dce.bayes.single.temp(conc, time, nriters, thin, burnin, tune,
                                  tau.gamma, tau.theta, ab.vp, ab.tauepsilon,
                                  aif.model, aif.parameter, vp, approx)
dce.bayes(conc, time, mask=array(TRUE,dim(conc)[1:(length(dim(conc))-1)]),
                      nriters=9500, thin=30, burnin=2000, tune=267,
                      tau.ktrans=1, tau.kep=tau.ktrans, ab.vp=c(1,19),
                      ab.tauepsilon=c(1,1/1000), aif.model=0,
                      aif.parameter=c(2.4,3,0.62,0.016), vp.do=TRUE,
                      samples=TRUE, approx=FALSE)
}
\arguments{
  \item{conc}{Matrix or array of concentration time series.
    Last dimension must be time.}
  \item{time}{Time in minutes.}
  \item{mask}{Mask matrix or array. Voxels with mask=0 will be excluded.}
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
}
\author{Volker Schmid}
\keyword{}
