\name{Relaxivity}
\alias{R10.lm}
\alias{E10.lm}
\alias{CA.fast}
\alias{CA.fast2}
\title{Estimate Intrinsic Tissue Relaxivity}
\description{
  Estimation of the intrinsic tissue relaxivity is achieved through
  nonlinear optimization and the dynamic signal intensities are
  converted into contrast agent concentration.
}
\usage{
R10.lm(signal, alpha, TR, guess, nprint=0)
E10.lm(signal, alpha, guess, nprint=0)
CA.fast(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4.39, verbose=FALSE)
CA.fast2(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4.39, verbose=FALSE)
}
\arguments{
  \item{signal}{is the vector of signal intensities as a function of
    flip angles.}
  \item{alpha}{is the vector of flip angles (in degrees).}
  \item{TR}{is the relaxation time (in seconds) used in the acquisition
    of the MRI data.}
  \item{guess}{is the vector of initial values for the parameters of
    interest: \eqn{m_0} and \eqn{R_{10}}.}
  \item{nprint}{is an integer, that enables controlled printing of
    iterates if it is positive.  In this case, estimates of \code{par}
    are printed at the beginning of the first iteration and every
    \code{nprint} iterations thereafter and immediately prior to return.
    If \code{nprint} is not positive, no tracing information on the
    progress of the optimization is produced.}
  \item{dynamic}{a multidimensional array of contrast agent
    concentrations.  The last dimension is assumed to be temporal, while
    the previous dimenions are assued to be spatial.}
  \item{dyn.mask}{is a (logical) multidimensional array that identifies
    the voxels to be analyzed.}
  \item{dangle}{is the flip angle used to acquire the dynamic MRI data.}
  \item{flip}{a multidimensional array of contrast agent
    concentrations.  The last dimension is assumed to be a function of
    the flip angles, while the previous dimenions are assued to be
    spatial.}
  \item{fangles}{is the vector of flip angles (in degrees).}
  \item{r1}{is the spin-lattice relaxivity constant (default = 4.39 for
    1.5T).  For 3T data it may be necessary to adjust this value.} 
  \item{verbose}{is a logical variable (default = \code{FALSE}) that
    allows text-based feedback during execution of the function.}
}
\details{
  The \code{E10.lm} and \code{R10.lm} functions estimate parameters for
  a vector of observed MR signal intensities, as a function of flip
  angle, using the following relationship
  \deqn{S(\alpha) = \frac{\sin(\alpha) \left(1 -
      \exp{-\text{TR}/\text{T}_1}\right)}{\left(1 - \cos(\alpha)
      \exp{-\text{TR}/\text{T}_1}\right)}.}
  The only difference between the two functions is exactly what is being
  estimated in the nonlinear least squares formulation.  They both
  require the function \code{\link[minpack.lm]{nls.lm}} that uses the
  Levenberg-Marquardt algorithm.

  The \code{CA.fast} function rearranges the assumed multidimensional
  (2D or 3D) structure of the multiple flip-angle data into a single
  matrix to take advantage of internal R functions instead of loops, and
  called \code{E10.lm}.  Conversion of the dynamic signal intensities to
  contrast agent concentration is performed via
  \deqn{[Gd] = \frac{1}{r_1}\left(\frac{1}{\text{T}_1} -
    \frac{1}{\text{T}_{10}}\right).} 

  The \code{CA2.fast} function assumes only two flip angles have been
  acquired and uses an approximation to the nonlinear relationship
  between signal intensity and flip angle enable to conversion from
  signal intensity to contrast agent concentration.
}
\value{
  A list structure is produced with estimates of \eqn{m_0}, \eqn{R_{10}}
  and information about the convergence of the nonlinear optimization
  routine.
}
\references{
  Buxton, R.B. (2002) \emph{Introduction to Functional Magnetic
    Resonance Imaging: Principles & Techniques}, Cambridge University
  Press: Cambridge, UK.

  Li, K.-L., Zhu, X.P., Waterton, J. and Jackson, A. (2000)
  Improved 3D quantiative mapping of blood volume and endothelial
  permeability in brain tumors, \emph{Journal of Magnetic Resonance
    Imaging}, \bold{12}, 347-357.

  Li, K.-L., Zhu, X.P., Kamaly-Asl, I.D., Checkley, D.R., Tessier,
  J.J.L., Waterton, J.C. and Jackson, A. (2000) Quantification of
  endothelial permeability, leakage space, and blood volume in brain
  tumors using combined T1 and T2* contrast-enhanced dynamic MR imaging,
  \emph{Journal of Magnetic Resonance Imaging}, \bold{11}, 575-585.
}
%\seealso{}
\examples{

}
\author{B. Whitcher}
\keyword{misc}
