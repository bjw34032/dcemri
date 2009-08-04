% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-07-15 09:44:00 (bjw34032)>
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
dcemri.spline(data, time, input, time.input=time,
mask=array(1,dim(data)[1:(length(dim(data))-1)]), k=4, knotpoints=25, 
knots=seq(min(time)-1e-3-(k-1)*(max(time)-min(time))/(knotpoints-k+1),1e-3+max(time)+k*(max(time)-min(time))/(knotpoints-k+1),(2e-3+max(time)-min(time))/(knotpoints-k+1)),
rw=2, nriters=500, thin=5, burnin=100, nlr=TRUE, hyper.a=1e-5,
hyper.b=1e-5, epsilon.a=1, epsilon.b=1e-5, silent=FALSE)
}
\arguments{
  \item{data}{}
  \item{time}{}
  \item{input}{}
  \item{time.input}{}
  \item{mask}{}
  \item{k}{}
  \item{knotpoints}{}
  \item{knots}{}
  \item{rw}{}
  \item{nriters}{}
  \item{thin}{}
  \item{burnin}{}
  \item{nlr}{}
  \item{hyper.a}{}
  \item{hyper.b}{}
  \item{epsilon.a}{}
  \item{epsilon.b}{}
  \item{t0}{}
  \item{t0.compute}{}
  \item{silent}{}
  \item{aif}{}
  \item{searchint}{}
  \item{do.plot}{}
}
\details{
  See Schmid \emph{et al.} for more details.
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
data("buckley")
img <- array(t(breast$data)[,(1:60)*5], c(13,1,1,60))
time <- seq(5,300,5) / 60
mask <- array(rep(TRUE, 13), c(13,1,1))
#Generate AIF params using the orton.exp function from Buckley's AIF
aifparams <- with(buckley, orton.exp.lm(time.min,input))
aifparams$D<-1 
pseudoobservedaif=aif.orton.exp(tt=time,AB=aifparams$AB,AG=aifparams$AG,muB=aifparams$muB,muG=aifparams$muG)

fit.spline <- dcemri.spline(img, time, mask, aif="fritz.hansen", nlr=TRUE)
fit.spline.aif <- dcemri.spline(img, time, mask, aif="observed", aif.observed=pseudoobservedaif, nlr=TRUE)

plot(breast$ktrans, fit.spline$ktrans, xlim=c(0,0.75), ylim=c(0,0.75),
     xlab="True Ktrans", ylab="Estimated Ktrans (Spline)")
points(breast$ktrans, fit.spline.aif$ktrans, pch=2)
legend(0.4, 0.2, pch=c(2,1), legend=c("Modelled AIF","Fritz Hansen"))
lines(c(-0.1,0.8), c(-0.1,0.8), lty=3, col="grey")

# not for windows, yet
system.time(fit <- dcemri.spline(img, time, mask, aif="fritz.hansen",nlr=TRUE,multicore=TRUE)

}
\author{Volker Schmid}
\keyword{models}
