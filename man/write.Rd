\name{Write Analyze}
\alias{write.analyze.img}
\title{Save Analyze Format Binary Image Files}
\description{
  This function saves both the header information and multidimensional
  image array into a pair of binary files in Analyze format.
}
\usage{
write.analyze.img(fname, hdr, img, type, gzipped=TRUE)
}
\arguments{
  \item{fname}{Pathname to save the Analzye pair of files .img and .hdr 
    without the suffix.}
  \item{hdr}{Header information as read in with \code{read.hdr}.}   
  \item{img}{Multidimensional array of image data.}
  \item{type}{Character string describing the image format.}
  \item{gzipped}{Allows the exportation of compressed (.gz) files.}
}
\value{
  \code{write.img} returns nothing.
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
\references{
  Analyze 7.5
  http://www.mayo.edu/bir/PDF/ANALYZE75.pdf
}
\seealso{
    \code{\link{read.hdr}}, \code{\link{read.img}}
}
%\examples{}
\author{Volker Schmid, Brandon Whitcher}
\keyword{file}
