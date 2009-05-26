\name{Apparent Diffusion Coefficient}
\alias{adc.lm}
\alias{ADC.fast}
\title{Estimate the Apparent Diffusion Coefficient (ADC)}
\description{
  Estimation of the apparent diffusion coefficient is achieved through
  nonlinear optimization.
}
\usage{
adc.lm(signal, b, guess, nprint=0)
ADC.fast(dwi, bvalues, dwi.mask, verbose=FALSE)
}
\arguments{
  \item{signal}{Signal intensity vector as a function of b-values.}
  \item{b,bvalues}{Diffusion weightings (b-values).}
  \item{guess}{Initial values of $S_0$ and $D$.}
  \item{nprint}{is an integer, that enables controlled printing of
    iterates if it is positive.  In this case, estimates of \code{par}
    are printed at the beginning of the first iteration and every
    \code{nprint} iterations thereafter and immediately prior to return.
    If \code{nprint} is not positive, no tracing information on the
    progress of the optimization is produced.}
  \item{dwi}{}
  \item{dwi.mask}{}
  \item{verbose}{}
}
\value{
  A list structure is produced with estimates of $S_0$, $D$ and
  information about the convergence of the nonlinear optimization
  routine.
}
\details{
  The \code{adc.lm} function requires the routine \code{nls.lm} that
  applies the Levenberg-Marquardt algorithm.
}
%\references{}
%\seealso{}
%\examples{}
\author{Brandon Whitcher}
\keyword{misc}
