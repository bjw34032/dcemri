% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-06-19 16:50:13 (bjw34032)>

\name{dce.spline}
\alias{dce.spline}
\alias{dce.spline.onset}
\alias{perfSpline}
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
dce.spline(data, time, input, time.input=time,
mask=array(1,dim(data)[1:(length(dim(data))-1)]), k=4, knotpoints=25, 
knots=seq(min(time)-1e-3-(k-1)*(max(time)-min(time))/(knotpoints-k+1),1e-3+max(time)+k*(max(time)-min(time))/(knotpoints-k+1),(2e-3+max(time)-min(time))/(knotpoints-k+1)),
rw=2, nriters=500, thin=5, burnin=100, nlr=TRUE, hyper.a=1e-5,
hyper.b=1e-5, epsilon.a=1, epsilon.b=1e-5, silent=FALSE)
dce.spline.onset(data, time, input, time.input=time,
mask=array(1,dim(data)[1:(length(dim(data))-1)]), k=4, knotpoints=25,
knots=seq(min(time)-1e-3-(k-1)*(max(time)-min(time))/(knotpoints-k+1),1e-3+max(time)+k*(max(time)-min(time))/(knotpoints-k+1),(2e-3+max(time)-min(time))/(knotpoints-k+1)),
rw=2, nriters=500, thin=5, burnin=100, nlr=TRUE, hyper.a=1e-5,
hyper.b=1e-5, epsilon.a=1, epsilon.b=1e-5, t0=0, t0.compute=TRUE)
perfSpline(data, time, aif, time.input, searchint=.25, do.plot=TRUE)
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

}
\author{Volker Schmid}
\keyword{models}
