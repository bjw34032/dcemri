\name{Read NIfTI}
\alias{read.nifti.hdr}
\alias{read.nifti.img}
\title{Read NIfTI Format Headers and Image Files}
\description{
  These functions read in the header information and multidimensional
  array from a binary file in NIfTI-1 format.
}
\usage{
read.nifti.hdr(fname, onefile=TRUE, gzipped=TRUE, verbose=FALSE, warn=-1)
read.nifti.img(fname, onefile=TRUE, gzipped=TRUE, verbose=FALSE, warn=-1)
}
\arguments{
  \item{fname}{Pathname of the NIfTI file(s).}
  \item{onefile}{Allows the scanning of single-file format (.nii) or
    dual-file format (.hdr and .img) NIfTI files.}
  \item{gzipped}{Allows the importation of compressed (.gz) files.}
  \item{verbose}{}
  \item{warn}{}
}
\value{
  \code{read.nifti.hdr} returns a list containing all the fields from
  the header file.  \code{read.nifti.img} returns a multidimensional
  array depending on the fields obtained from the header file.
}
\details{
  Both functions utilize the internal \code{readBin} and \code{readChar}
  command to efficiently extract information from the binary file(s).
  The types of data are limited to 1- and 2-byte integers, 4-byte floats
  and 8-byte doubles.
}
\references{
  NIfTI-1
  http://nifti.nimh.nih.gov/
}
\seealso{
  \code{\link{read.analyze.hdr}},
  \code{\link{read.analyze.img}},
  \code{\link{write.nifti.img}}
}
\examples{
}
\author{Brandon Whitcher}
\keyword{file}
