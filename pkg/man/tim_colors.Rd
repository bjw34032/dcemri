% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-06-19 17:10:27 (bjw34032)>

\name{tim.colors}
\alias{tim.colors}
\title{Tim's Useful Color Table}
\description{
  A pleasing rainbow style color table patterned after that used in 
  Matlab by Tim.
}
\usage{
tim.colors(n=64)
}
\arguments{
  \item{n}{is the number of color levels (default = 64).}
}
\details{
  Based on the \code{\link[fields]{tim.colors}} function in the
  \pkg{fields} package.  The \code{tim.colors} function here has been
  modified to break any dependence on code in the \pkg{fields} package.
  Spline interpolation (\code{interpSpline}) is used when the number of
  requested colors is not the default.
}
\value{
  A vector of character strings giving the colors in hexadecimal
  format.
}
\seealso{
  \code{\link{topo.colors}}, \code{\link{terrain.colors}}
}
\examples{
tim.colors(10) 
image(outer(1:20,1:20,"+"), col=tim.colors(75), main="tim.colors")
}
\keyword{aplot}
