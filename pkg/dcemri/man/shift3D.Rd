% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-07-14 17:49:45 (bjw34032)>
% $Id$

\name{shift3D}
\alias{shift3D}
\title{Shift a 3D Array in One Dimension}
\description{
  One axis of the three-dimensional array is translated by an integer
  amount.  This is useful when applying convolution operators in the
  Fourier domain.
}
\usage{
shift3D(A, s, type, fill=0)
}
\arguments{
  \item{A}{is a three-dimensional array.}
  \item{s}{is the integer number of translation steps.}
  \item{type}{is a character string using anatomical coordinates
    assuming a transverse acquisition scheme (\dQuote{LR} = left-right =
    x-axis, \dQuote{AP} = anterior-posterior = y-axis, \dQuote{SI} =
    superior-inferior = z-axis).}
  \item{fill}{is the quantity used to fill gaps induced by the
    translations (circular boundary conditions are NOT used).}
}
%\details{}
\value{
  A three-dimensional array is returned, the same dimension as the
  original array, with one dimension translated.
}
%\references{}
\author{Brandon Whitcher}
\seealso{
  \code{\link{conv.fft}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
