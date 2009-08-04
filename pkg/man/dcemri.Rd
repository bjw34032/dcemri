% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-07-14 18:11:40 (bjw34032)>
% $Id$

\name{dcemri}
\alias{dcemri.lm}
\title{Pharmacokinetic Models for Dynamic Contrast-Enhanced MRI Data}
\description{
  Parameter estimation for single compartment models is performed using
  literature-based or user-specified arterial input functions.  The
  Levenburg-Marquardt algorithm does the heavy lifting.
}
\usage{
dcemri.lm(conc, time, img.mask, model="extended", aif="tofts.kermode",
          nprint=0, user=NULL, ...)
}
\arguments{
  \item{conc}{is a multidimensional (1d to 4d) array of contrast agent
    concentrations.  The last dimension is assumed to be temporal, while
    the previous dimensions are assued to be spatial.}
  \item{time}{is a vector of acquisition times (in minutes) relative to 
    injection of the contrast agent.  Negative values should be used
    prior to the injection.}
  \item{img.mask}{is a (logical) multidimensional array that identifies
    the voxels to be analyzed. Has to have same dimension as \code{conc}
    minus temporal dimension.}
  \item{model}{is a character string that identifies the type of
    compartmental model to be used.  Acceptable models include:
    \dQuote{weinmann} Tofts & Kermode AIF convolved with single compartment
	 model 
    \dQuote{extended} (default) Weinmann model extended with additional
	vascular compartment, ...
    \dQuote{orton.exp} extended model using Orton's exponential AIF
    \dQuote{orton.cos} extended model using Orton's raised cosine AIF
    }
  \item{aif}{is a character string that identifies the type of arterial
    input function (AIF) to be used.  Acceptable AIF models include:
    \code{tofts.kermode}, \code{fritz.hansen}, \code{orton.exp},
    \code{orton.cos}, and \code{user}.  All AIF models set the
    parametric form and parameter values; except \code{user} where a set
    of user-defined parameter values are allowed for the
    \code{orton.exp} model only!}
  \item{nprint}{is an integer, that enables controlled printing of
    iterates if it is positive.  In this case, estimates of \code{par}
    are printed at the beginning of the first iteration and every
    \code{nprint} iterations thereafter and immediately prior to return.
    If \code{nprint} is not positive, no tracing information on the
    progress of the optimization is produced.}
  \item{user}{is a list with the following paramters required: D, AB,
    muB, AG, muG.}
  \item{...}{Additional parameters to the function.}
}
\details{
  Compartmental models are the solution to the modified general rate
  equation (Kety 1951).  The specific parametric models considered here
  include the basic Kety model
  \deqn{C_t(t)=K^{trans}\left[C_p(t)\otimes\exp(-k_{ep}t)\right],}
  where \eqn{\otimes} is the convoluation operator, and the so-called
  extended Kety model
  \deqn{C_t(t)=v_pC_p(t)+K^{trans}\left[C_p(t)\otimes\exp(-k_{ep}t)\right].}
  The arterial input function must be either literature-based (with
  fixed parameters) or the exponential AIF from Orton \emph{et al.}
  (2008) with user-defined parameters.
}
\value{
  Parameter estimates and their standard errors are provided for the
  masked region of the multidimensional array.  They include
  \item{ktrans}{Transfer rate from plasma to the extracellular,
    extravascular space (EES).}
  \item{kep}{Rate parameter for transport from the EES to plasma.}
  \item{ve}{Fractional occupancy by EES (the ratio between ktrans and
    kep).}
  \item{vp}{Fractional occupancy by plasma.}
  \item{ktranserror}{Standard error for ktrans.}
  \item{keperror}{Standard error for kep.}
  \item{vperror}{Standard error for vp.}
  The residual sum-of-squares is also provided, along with the original
  acquisition times (for plotting purposes).
}
\references{
  Ahearn, T.S., Staff, R.T., Redpath, T.W. and Semple, S.I.K. (2005) The
  use of the Levenburg-Marquardt curve-fitting algorithm in
  pharmacokinetic modelling of DCE-MRI data, \emph{Physics in Medicine
    and Biology}, \bold{50}, N85-N92.

Fritz-Hansen, T., Rostrup, E., Larsson, H.B.W, Sondergaard, L., Ring,
  P. and Henriksen, O. (1993) Measurement of the arterial concentration
  of Gd-DTPA using MRI: A step toward quantitative perfusion imaging,
  \emph{Magnetic Resonance in Medicine}, \bold{36}, 225-231.

  Orton, M.R., Collins, D.J., Walker-Samuel, S., d'Arcy, J.A.,
  Hawkes, D.J., Atkinson, D. and Leach, M.O. (2007) Bayesian
  estimation of pharmacokinetic parameters for DCE-MRI with a robust
  treatment of enhancement onset time, \emph{Physics in Medicine and
  Biology} \bold{52}, 2393-2408.
  
  Orton, M.R., d'Arcy, J.A., Walker-Samuel, S., Hawkes, D.J.,
  Atkinson, D., Collins, D.J. and Leach, M.O. (2008) Computationally
  efficient vascular input function models for quantitative kinetic
  modelling using DCE-MRI, \emph{Physics in Medicine and Biology}
  \bold{53}, 1225-1239.

  Tofts, P.S., Brix, G, Buckley, D.L., Evelhoch, J.L., Henderson, E.,
  Knopp, M.V., Larsson, H.B.W., Lee, T.-Y., Mayr, N.A., Parker, G.J.M.,
  Port, R.E., Taylor, J. and Weiskoff, R. (1999) Estimating kinetic
  parameters from dynamic contrast-enhanced \eqn{T_1}-weighted MRI of a
  diffusable tracer: Standardized quantities and symbols, \emph{Journal
    of Magnetic Resonance}, \bold{10}, 223-232.

  Tofts, P.S. and Kermode, A.G. (1984) Measurement of the blood-brain
  barrier permeability and leakage space using dynamic MR
  imaging. 1. Fundamental concepts, \emph{Magnetic Resonance in
    Medicine}, \bold{17}, 357-367.

  Weinmann, H.J., Laniado, M. and Mutzel, W. (1984) Pharmacokinetics of
  Gd-DTPA/dimeglumine after intraveneous injection into healthy
  volunteers, \emph{Physiological Chemistry and Physics and Medical
    NMR}, \bold{16}, 167-172. 
}
\seealso{
  \code{\link{dcemri.bayes}}, \code{\link{dcemri.spline}},
  \code{\link[minpack.lm]{nls.lm}}
}
\examples{
data("buckley")
#Set up breast data for dcemri
img <- array(t(breast$data)[,(1:60)*5], c(13,1,1,60))
time <- seq(5, 300, by=5) / 60
mask <- array(rep(TRUE, 13), c(13,1,1))
#Generate AIF params using the orton.exp function from Buckley's AIF
aifparams <- with(buckley, orton.exp.lm(time.min,input))
aifparams$D<-1
fit <- dcemri.lm(img, time, mask, model="orton.exp", aif="user", user=aifparams)
}
\author{Brandon Whitcher, Volker Schmid}
\keyword{models}
