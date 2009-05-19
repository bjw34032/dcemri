\name{Relaxivity}
\alias{R10.lm}
\alias{E10.lm}
\alias{CA.fast}
\alias{CA.fast2}
\title{Estimate Intrinsic Tissue Relaxivity}
\description{
  Estimation of the intrinsic tissue relaxivity is achieved through
  nonlinear optimization.
}
\usage{
R10.lm(signal, alpha, TR, guess, nprint=0)
E10.lm(signal, alpha, guess, nprint=0)
CA.fast(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4.39, verbose=FALSE)
CA.fast2(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4.39, verbose=FALSE)
}
\arguments{
  \item{signal}{Signal intensity vector as a function of flip angles.}
  \item{alpha}{Vector of angles (in degrees).}
  \item{TR}{Relaxation time (in seconds) from the magnetic resonance
    scan.}
  \item{guess}{Initial values of $m_0$ and $R_{10}$.}
  \item{nprint}{is an integer, that enables controlled printing of
    iterates if it is positive.  In this case, estimates of \code{par}
    are printed at the beginning of the first iteration and every
    \code{nprint} iterations thereafter and immediately prior to return.
    If \code{nprint} is not positive, no tracing information on the
    progress of the optimization is produced.}
  \item{dynamic}{}
  \item{dyn.mask}{}
  \item{dangle}{}
  \item{flip}{}
  \item{fangles}{}
  \item{r1}{}
  \item{verbose}{}
}
\value{
  A list structure is produced with estimates of $m_0$, $R_{10}$ and
  information about the convergence of the nonlinear optimization
  routine.
}
\details{
  The ".lm" functions require the function \code{nls.lm} that uses the
  Levenberg-Marquardt algorithm.
}
\references{
  K.-L. Li, X. P. Zhu, J. Waterton and A. Jackson (2000)
  Improved {3D} quantiative mapping of blood volume and endothelial
  permeability in brain tumors,
  \emph{Journal of Magnetic Resonance Imaging},
  \bold{12}, 347-357.

  K.-L. Li, X. P. Zhu, I. D. Kamaly-Asl, D. R. Checkley,
  J. J. L. Tessier, J. C. Waterton and A. Jackson (2000)
  Quantification of endothelial permeability, leakage space, and blood
  volume in brain tumors using combined T1 and T2* contrast-enhanced
  dynamic {MR} imaging,
  \emph{Journal of Magnetic Resonance Imaging},
  \bold{11}, 575-585.
}
%\seealso{}
%\examples{}
\author{B. Whitcher}
\keyword{misc}
