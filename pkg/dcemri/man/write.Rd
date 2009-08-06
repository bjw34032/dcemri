% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-07-14 17:53:16 (bjw34032)>
% $Id$

\name{Write Output as Analyze/NIfTI}
\alias{write.analyze.img}
\alias{write.nifti.img}
\title{Write Analyze/NIfTI Format Binary Image Files}
\description{
  This function saves both the header information and multidimensional
  image array into a pair of binary files in Analyze format or a single
  binary file in NIfTI format.
}
\usage{
write.analyze.img(fname, hdr, img, type, gzipped=TRUE, warn=-1)
write.nifti.img(fname, hdr, img, type, gzipped=TRUE, warn=-1)
}
\arguments{
  \item{fname}{is the pathname to save the Analzye pair (.img and .hdr)
    or NIfTI file (.nii) \bold{without} the suffix.}
  \item{hdr}{is the appropriate header object information.}   
  \item{img}{is the multidimensional array of image data.}
  \item{type}{is a character string describing the image format.  Valid
    character strings include:\dQuote{uint1}, \dQuote{uint8},
    \dQuote{int16}, \dQuote{int32}, \dQuote{float32}, \dQuote{float64}.}
  \item{gzipped}{is a character string that enables exportation of
    compressed (.gz) files (default = \code{TRUE}).}
  \item{warn}{is a number to regulate the display of warnings (default =
    -1).  See \code{\link{options}} for more details.}
}
\details{
  Both functions utilize the internal \code{writeBin} and
  \code{writeChar} command to write infromation to a binary file.

  Current acceptable data types include
  \itemize{
    \item{\dQuote{uint1}}{DT BINARY (1 bit per voxel)}
    \item{\dQuote{uint8}}{DT UNSIGNED CHAR (8 bits per voxel)}
    \item{\dQuote{int16}}{DT SIGNED SHORT (16 bits per voxel)}
    \item{\dQuote{int32}}{DT SINGED INT (32 bits per voxel)}
    \item{\dQuote{float}}{DT FLOAT (32 bits per voxel)}
    \item{\dQuote{double}}{DT DOUBLE (64 bits per voxel)}
    }
}
\value{
  \code{write.analyze.img} and \code{write.nifti.img} return nothing.
}
\references{
  Analyze 7.5
  \url{http://www.mayo.edu/bir/PDF/ANALYZE75.pdf}
  
  NIfTI-1
  \url{http://nifti.nimh.nih.gov/}
}
\seealso{
    \code{\link{read.hdr}}, \code{\link{read.img}}
}
\examples{
norm<-dnorm(seq(-5,5,length=32),sd=1)
norm<-(norm-min(norm))/max(norm-min(norm))
img<-norm%o%norm%o%norm
img<-round(img*255)
X<-dim(img)[1]
Y<-dim(img)[2]
Z<-dim(img)[3]

## Analyze
hdr <- make.hdr(X, Y, Z, 1, "int16")
write.analyze.img("xyz64", hdr, img, "int16")
## NIfTI
nhdr <- make.hdr(X, Y, Z, 1, "INT", "nifti")
write.nifti.img("xyz64", nhdr, img, "int16")
## These files should be viewable in, for example, FSLview.
## Make sure you adjust the min/max values for proper visualization.
data <- read.img("xyz64.nii.gz")
}
\author{Brandon Whitcher, Volker Schmid}
\keyword{file}
