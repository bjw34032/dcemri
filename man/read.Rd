\name{Read Analyze}
\alias{read.analyze.hdr}
\alias{read.analyze.img}
\alias{make.hdr}
\alias{read.hdr}
\alias{read.img}
\title{Read Analyze Format Headers and Binary Image Files}
\description{
  These functions read in the header information and multidimensional
  array from a binary file in Analyze 7.5 format.
}
\usage{
read.analyze.hdr(fname, gzipped=TRUE, verbose=FALSE, warn=-1)
read.analyze.img(fname, gzipped=TRUE, signed=FALSE, verbose=FALSE, warn=-1)
make.hdr(X, Y, Z, T, datatype, type)
read.hdr(fname, verbose=FALSE, warn=-1)
read.img(fname, verbose=FALSE, warn=-1)
}
\arguments{
  \item{fname}{Pathname of the Analyze pair of files .img and .hdr
    without the suffix.}
  \item{signed}{Only changed to \code{FALSE} when reading in unsigned
    single-byte integers.}
  \item{gzipped}{Allows the importation of compressed (.gz) files.}
  \item{verbose}{}
  \item{warn}{}
  \item{X,Y,Z,T}{}
  \item{datatype}{}
  \item{type}{}
}
\details{
  Both functions utilize the internal \code{readBin} and \code{rawToChar}
  command to efficiently extract information from a binary file.  The
  types of data are limited to 1- and 2-byte integers, 4-byte floats and
  8-byte doubles.

  The function \code{read.hdr} decides whether the input file should be
  read using \code{read.analyze.hdr} or \code{read.nifti.hdr} and if
  compression is required.  A similar description is valid for
  \code{read.img}. 
}
\value{
  \code{read.analyze.hdr} returns a list containing all the fields from
  the header file.  \code{read.analyze.img} returns a multidimensional
  array depending on the field read in from the header file.
}
\references{
  Analyze 7.5
  \url{http://www.mayo.edu/bir/PDF/ANALYZE75.pdf}
}
\seealso{
  \code{\link{read.nifti.img}},
  \code{\link{write.analyze.img}},
  \code{\link{write.nifti.img}}
}
\author{Brandon Whitcher, Volker Schmid}
\keyword{file}
