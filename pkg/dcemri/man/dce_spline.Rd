% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-08-06 14:22:21 (bjw34032)>
% $Id$

\name{dcemri.spline}
\alias{dcemri.spline}
\title{Bayesian P-Splines for Dynamic Contrasat-Enhanced MRI Data}
\description{
  Quantitative analysis of DCE-MRI typically involves the convolution of
  an arterial input function (AIF) with a nonlinear pharmacokinetic
  model of the contrast agent concentration.  This function takes a
  semi-parametric penalized spline smoothing approach, with which the
  AIF is convolved with a set of B-splines to produce a design matrix
  using locally adaptive smoothing parameters based on Bayesian
  penalized spline models (P-splines).
} 
\usage{
dcemri.spline(conc, time, img.mask, time.input=time, model="weinmann", 
              aif="tofts.kermode", user=NULL, aif.observed=NULL,
              nriters=500, thin=5, burnin=100, 
              ab.hyper=c(1e-5,1e-5), ab.tauepsilon=c(1,1/1000), 
              k=4,p=25,rw=2,
              knots=NULL, 
              nlr = FALSE, t0.compute=FALSE, samples=FALSE, multicore=FALSE,
\dots) 
}
\arguments{
  \item{conc}{An array of Gd concentration}
  \item{time}{Time points of aquisition of Gd concentration}
  \item{img.mask}{Array of voxels to fit}
  \item{time.input}{Time points of observed Arterial Gd concentration, defaults to time}
  \item{model}{The type of compartmental model to be used. Acceptable models
    include: \dQuote{AATH} or \dQuote{weinmann} (default).}
  \item{aif}{Arterial input function to use. Values include:
    \dQuote{tofts.kermode}, \dQuote{fritz.hansen} or \dQuote{observed}.
    If \dQuote{observed} you must provide the observed concentrations in
    aif.observed.}
  \item{aif.observed}{Arterial concentrations observed at timepoints
    time.input}
  \item{multicore}{Use multicore library}
  \item{nlr}{Return the generated samples}
  \item{user}{}
  \item{ab.hyper}{}
  \item{ab.tauepsilon}{}
  \item{p}{}
  \item{t0.compute}{}
  \item{samples}{}
  \item{k}{}
  \item{knots}{}
  \item{rw}{}
  \item{nriters}{}
  \item{thin}{}
  \item{burnin}{}
  \item{...}{}
}
\details{
  See Schmid \emph{et al.} (2009) for more details.
}
\value{
  
}
\references{
  Schmid, V., Whitcher, B., Padhani, A.R. and G.-Z. Yang (2009) A
  semi-parametric technique for the quantitative analysis of dynamic
  contrast-enhanced MR images based on Bayesian P-splines, \emph{IEEE
    Transactions on Medical Imaging}, \bold{28} (6), 789-798.
}
\seealso{
  \code{\link{dcemri.lm}}, \code{\link{dcemri.bayes}}
}
\examples{
data(buckley)
img <- array(t(breast$data)[,(1:60)*5], c(13,1,1,60))
time.min <- buckley$time.min[(1:60)*5]
mask <- array(rep(TRUE, 13), c(13,1,1))
## Generate AIF params using the orton.exp function from Buckley's AIF
aifparams <- with(buckley, orton.exp.lm(time.min,input))
aifparams$D <- 1 
pseudoobservedaif <- aif.orton.exp(tt=time.min, AB=aifparams$AB,
                                   AG=aifparams$AG, muB=aifparams$muB,
                                   muG=aifparams$muG)

fit.spline <- dcemri.spline(img, time.min, mask, aif="fritz.hansen", nlr=TRUE)
fit.spline.aif <- dcemri.spline(img, time.min, mask, aif="observed",
                                aif.observed=pseudoobservedaif, nlr=TRUE)

plot(breast$ktrans, fit.spline$ktrans, xlim=c(0,0.75), ylim=c(0,0.75),
     xlab="True Ktrans", ylab="Estimated Ktrans (Spline)")
points(breast$ktrans, fit.spline.aif$ktrans, pch=2)
legend(0.4, 0.2, pch=c(2,1), legend=c("Modelled AIF","Fritz Hansen"))
lines(c(-0.1,0.8), c(-0.1,0.8), lty=3, col="grey")
cbind(breast$ktrans, fit.spline$ktrans[,,1], fit.spline.aif$ktrans[,,1])
}
\author{Volker Schmid}
\keyword{models}
