% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-07-14 17:59:01 (bjw34032)>
% $Id$

\name{conv.fft}
\alias{conv.fft}
\title{Convolution of 3D Arrays using the Fourier Transform}
\description{
  Convolve a 
}
\usage{
conv.fft(A, B, C, FFTA=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{is a three-dimensional array ("the template").}
  \item{B}{is a three-dimensional array ("the target").}
  \item{C}{is a three dimensional array.}
  \item{FFTA}{is the three-dimensional Fourier transform of \code{A},
    this may save time when looping over multiple \dQuote{targets}.}
}
\details{
  The arrays \eqn{A} and \eqn{B} are transformed into the Fourier domain
  and multiplied together (equivalent to a convolution in the image
  domain across all spatial locations simultaneously).  
}
\value{
  A three-dimensional array, the same dimension as the input arrays,
  that is the convolution of the \dQuote{target} to the
  \dQuote{template} at all spatial locations.
}
\references{
  Briggs, W.L. and Henson, V.E. (1995) \emph{The DFT: An Owner's Manual
    for the Discrete Fourier Transform}, SIAM: Philadelphia.
}
\author{Brandon Whitcher}
\seealso{
  \code{\link{fft}}, \code{\link{ftm}}, \code{\link{shift3D}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
