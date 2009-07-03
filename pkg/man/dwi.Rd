% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-06-19 17:04:28 (bjw34032)>

\name{Apparent Diffusion Coefficient}
\alias{adc.lm}
\alias{ADC.fast}
\title{Estimate the Apparent Diffusion Coefficient (ADC)}
\description{
  Estimation of apparent diffusion coefficient (ADC) values, using a
  single exponential function, is achieved through nonlinear
  optimization.
}
\usage{
adc.lm(signal, b, guess, nprint=0)
ADC.fast(dwi, bvalues, dwi.mask, verbose=FALSE)
}
\arguments{
  \item{signal}{Signal intensity vector as a function of b-values.}
  \item{b,bvalues}{Diffusion weightings (b-values).}
  \item{guess}{Initial values of \eqn{S_0} and \eqn{D}.}
  \item{nprint}{is an integer, that enables controlled printing of
    iterates if it is positive.  In this case, estimates of \code{par}
    are printed at the beginning of the first iteration and every
    \code{nprint} iterations thereafter and immediately prior to return.
    If \code{nprint} is not positive, no tracing information on the
    progress of the optimization is produced.}
  \item{dwi}{Multidimensional array of diffusion-weighted images.}
  \item{dwi.mask}{Logical array that defines the voxels to be analyzed.}
  \item{verbose}{Additional information will be printed when
    \code{verbose=TRUE}.}
}
\details{
  The \code{adc.lm} function estimates parameters for a vector of
  observed MR signal intensities using the following relationship
  \deqn{S(b) = S_0 \exp{-bD},}
  where \eqn{S_0} is the baseline signal intensity and \eqn{D} is the
  apparent diffusion coefficient (ADC).  It requires the routine
  \code{nls.lm} that applies the Levenberg-Marquardt algorithm.  Note,
  low b-values (\eqn{<50} or \eqn{<100} depending on who you read)
  should be avoided in the parameter estimation because they do not
  represent information about the diffusion of water in tissue.

  The \code{ADC.fast} function rearranges the assumed multidimensional
  (2D or 3D) structure of the DWI data into a single matrix to take
  advantage of internal R functions instead of loops, and called
  \code{adc.lm}.
}
\value{
  A list structure is produced with estimates of \eqn{S_0}, \eqn{D} and
  information about the convergence of the nonlinear optimization
  routine.
}
\references{
  Buxton, R.B. (2002) \emph{Introduction to Functional Magnetic
    Resonance Imaging: Principles & Techniques}, Cambridge University
  Press: Cambridge, UK.

  Callahan, P.T. (2006) \emph{Principles of Nuclear Magnetic Resonance
    Microscopy}, Clarendon Press: Oxford, UK.
  
  Koh, D.-M. and Collins, D.J. (2007) Diffusion-Weighted MRI in the
  Body: Applications and Challenges in Oncology, \emph{American Journal
    of Roentgenology}, \bold{188}, 1622-1635.
}
\seealso{
  \code{\link[minpack.lm]{nls.lm}}
}
\examples{

}
\author{Brandon Whitcher}
\keyword{models}
