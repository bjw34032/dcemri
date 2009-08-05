% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-07-29 15:59:19 (bjw34032)>
% $Id$

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
	              nriters=7000, thin=10, burnin=2000, tune=267, 
	              tau.ktrans=1, tau.kep=tau.ktrans, ab.vp=c(1,19),
	              ab.tauepsilon=c(1,1/1000), samples=FALSE, multicore=FALSE
	              ...)
}
\arguments{
  \item{conc}{Matrix or array of concentration time series.
    Last dimension must be time.}
  \item{time}{Time in minutes.}
  \item{img.mask}{Mask matrix or array. Voxels with mask=0 will be
    excluded.}
  \item{model}{is a character string that identifies the type of
    compartmental model to be used.  Acceptable models include:
    \dQuote{weinmann} Tofts & Kermode AIF convolved with single compartment
	 model 
    \dQuote{extended} (default) Weinmann model extended with additional
	vascular compartment, ...
    }
  \item{aif}{Model for arterial input function. So far implemented:
    0 or \dQuote{tofts-kermode} double exponential model by Tofts and
    Kermode; 1 or \dQuote{orton} exponential model by Orton \emph{et
      al}.}
  \item{user}{Vector of AIF parameters. 
    For Tofts and Kermode: \eqn{a_1},\eqn{m_1},\eqn{a_2},\eqn{m_2}; for
    Orton \emph{et al.}: \eqn{A_b},\eqn{\mu_b},\eqn{A_g},\eqn{\mu_g}.}
  \item{nriters}{Total number of iterations.}
  \item{thin}{Thining factor.}
  \item{burnin}{Number of iterations for burn-in.}
  \item{tune}{Number for iterations for tuning. Algorithm will be 
    tuned to an acceptance rate of 0.3 to 0.6.}
  \item{tau.ktrans}{Variance parameter for prior on \eqn{\log(K^{trans})}.}
  \item{tau.kep}{Variance parameter for prior on \eqn{\log(k_{ep})}.}
  \item{ab.vp}{Hyper-prior parameters for \eqn{v_p}{vp} Beta prior.}
  \item{ab.tauepsilon}{Hyper-prior parameters for observation error
    Gamma prior.} 
  \item{samples}{If \code{TRUE} output includes samples drawn
    from the posterior distribution for all parameters.}
  \item{multicore}{If \code{TRUE} algorithm is parallelised using multicore (requires multicore package and R > 2.9.0).}
}
\details{
  See Schmid \emph{et al.} (2006) for details.
}
\value{
  Parameter estimates and their standard errors are provided for the
  masked region of the multidimensional array.  They include
  \item{ktrans}{Transfer rate from plasma to the extracellular,
    extravascular space (EES).}
  \item{ktranserror}{Error on \eqn{K^{trans}}{ktrans}.}
  \item{kep}{Rate parameter for transport from the EES to plasma.}
  \item{keperror}{Error on \eqn{k_{ep}}{kep}.}
  \item{ve}{Fractional occupancy by EES (the ratio between ktrans and
    kep).}
  \item{vperror}{Error on \eqn{v_e}{ve}.}
  \item{vp}{Fractional occupancy by plasma.}
  \item{sigma2}{The residual sum-of-squares from the model fit.}
  \item{time}{Acquisition times (for plotting purposes).}
  Note, not all parameters are available under all models choices.
}
\references{
  Schmid, V., Whitcher, B., Padhani, A.R., Taylor, N.J. and Yang, G.-Z. 
  (2006) Bayesian methods for pharmacokinetic models in dynamic
  contrast-enhanced magnetic resonance imaging, \emph{IEEE Transactions
    on Medical Imaging}, \bold{25} (12), 1627-1636.
}
\seealso{
  \code{\link{dcemri.lm}}, \code{\link{dcemri.spline}}
}
\examples{
data("buckley")
img <- array(t(breast$data)[, (1:60)*5], c(13, 1, 1, 60))
time.min <- buckley$time.min[(1:60)*5]
mask <- array(rep(TRUE, 13), c(13,1,1))
#Generate AIF params using the orton.exp function from Buckley's AIF
aifparams <- with(buckley, orton.exp.lm(time.min,input))
aifparams$D<-1

fit.bayes <- dcemri.bayes(img, time.min, mask, aif="fritz.hansen")
fit.bayes.aif <- dcemri.bayes(img, time.min, mask, model="orton.exp", aif="user", user=aifparams)
plot(breast$ktrans, fit.bayes$ktrans, xlim=c(0,0.75), ylim=c(0,0.75),
     xlab="True Ktrans", ylab="Estimated Ktrans (Bayes)")
points(breast$ktrans, fit.bayes.aif$ktrans, pch=2)
legend(0.4, 0.2, pch=c(2,1), legend=c("Orton Exp (modelled)","Fritz Hansen"))
lines(c(-0.1,0.8), c(-0.1,0.8), lty=3, col="grey")


fit.lm <- dcemri.lm(img, time.min, mask, aif="fritz.hansen")
fit.lm.aif <- dcemri.lm(img, time.min, mask, model="orton.exp", aif="user", user=aifparams)
plot(breast$ktrans, fit.bayes$ktrans, xlim=c(0,0.75), ylim=c(0,0.75),
     xlab="True Ktrans", ylab="Estimated Ktrans")
points(breast$ktrans, fit.bayes.aif$ktrans, pch=3)
points(breast$ktrans, fit.lm$ktrans, pch=2)
points(breast$ktrans, fit.lm.aif$ktrans, pch=4)
legend(0.6, 0.2, pch=c(4,3,2,1), legend=c("LM OE","Bayes OE","LM FH","Bayes FH"))
lines(c(-0.1,0.8), c(-0.1,0.8), lty=3, col="grey")
cbind(breast$ktrans, fit.bayes$ktrans[,,1], fit.bayes.aif$ktrans[,,1], fit.lm$ktrans[,,1], fit.lm.aif$ktrans[,,1])

# not for windows, yet
\dontrun{
system.time(fit.bayes <- dcemri.bayes(img, time.min, mask, aif="fritz.hansen", samples=TRUE, multicore=TRUE))
}
}
\author{Volker Schmid}
\keyword{models}
