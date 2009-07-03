% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-06-19 17:05:47 (bjw34032)>

\name{Miscellaneous Functions for Analyze/NIfTI Data}
\alias{nsli}
\alias{NSLI}
\alias{ntim}
\alias{NTIM}
\title{Miscellaneous Functions for Analyze/NIfTI Data}
\description{
  Functions to extract the dimensions from Analyze/NIfTI data.
}
\usage{
nsli(x)
ntim(x)
NSLI(x)
NTIM(x)
}
\arguments{
  \item{x}{Three- or four-dimensional array (e.g., read in from an
    Analyze/NIfTI file).}
}
\details{
  Simple calls to \code{dim} to replicate the functionality of
  \code{nrow} and \code{ncol} for higher dimensions of an array that are
  commonly required when manipulating medical imaging data.
}
\value{
  Third (slice) or fourth (time) dimension of the array.
}
\references{
  
}
\seealso{
  \code{\link{read.hdr}}, \code{\link{read.img}}
}
\examples{

}
\author{B. Whitcher}
\keyword{misc}
