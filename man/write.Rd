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
write.analyze.img(fname, hdr, img, type, gzipped=TRUE)
write.nifti.img(fname, hdr, img, type, gzipped=TRUE)
}
\arguments{
  \item{fname}{is the pathname to save the Analzye pair (.img and .hdr)
    or NIfTI file (.nii) \bold{without} the suffix.}
  \item{hdr}{is the appropriate header object information.}   
  \item{img}{is the multidimensional array of image data.}
  \item{type}{is a character string describing the image format.  Valid
    character strings include: "uint1", "uint8", "int16", "int32",
    "float32", "float64".}
  \item{gzipped}{is a character string that enables exportation of
    compressed (.gz) files (default = \code{TRUE}).}
}
\details{
  Both functions utilize the internal \code{writeBin} and
  \code{writeChar} command to write infromation to a binary file.

  Current acceptable data types include
  \item{"uint1"}{DT BINARY (1 bit per voxel)}
  \item{"uint8"}{DT UNSIGNED CHAR (8 bits per voxel)}
  \item{"int16"}{DT SIGNED SHORT (16 bits per voxel)}
  \item{"int32"}{DT SINGED INT (32 bits per voxel)}
  \item{"float"}{DT FLOAT (32 bits per voxel)}
  \item{"double"}{DT DOUBLE (64 bits per voxel)}
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
X <- Y <- Z <- 64
img <- array(sample(c(0,1), X*Y*Z, replace=TRUE), c(X,Y,Z))
## Analyze
hdr <- make.hdr(X, Y, Z, 1, "int16")
write.analyze.img("xyz64", hdr, img, "int16")
## NIfTI
nhdr <- make.hdr(X, Y, Z, 1, "INT", "nifti")
write.nifti.img("XYZ64", nhdr, img, "int16")
## These files should be viewable in, for example, FSLview.
## Make sure you adjust the min/max values for proper visualization.
data <- read.img("XYZ64.nii.gz")
}
\author{Brandon Whitcher, Volker Schmid}
\keyword{file}
