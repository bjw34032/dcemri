\name{Write NIfTI}
\alias{write.nifti.img}
\title{Save NIfTI Format Header/Image Files}
\description{
  This function writes the appropriate header and multidimensional array
  information to binary file in NIfTI-1 format.
}
\usage{
write.nifti.img(fname, hdr, img, type, gzipped=TRUE)
}
\arguments{
  \item{fname}{Pathname of the NIfTI file(s).}
  \item{hdr}{Valid NIfTI header object.}
  \item{img}{Multidimensional array containing the image information.}
  \item{type}{Output type for the image data.  Valid character strings
    include: "UINT1", "UINT8", "INT16", "INT32", "FLOAT32", "FLOAT64".}
  \item{gzipped}{Allows the export of compressed (.gz) file(s).}
}
\value{
  \code{write.nifti.img} returns nothing, but writes a file.
}
\details{
  The function utilizes internal \code{writeBin} and \code{writeChar}
  commands to export information to the binary file(s).  The types of
  data are limited to 1-bit integers, 1-, 2- and 4-byte integers, 4-byte
  floats and 8-byte doubles.
}
\references{
  NIfTI-1
  http://nifti.nimh.nih.gov/
}
\seealso{
  \code{\link{write.analyze.img}}, \code{\link{read.nifti.img}}
}
\examples{
}
\author{Brandon Whitcher}
\keyword{file}
