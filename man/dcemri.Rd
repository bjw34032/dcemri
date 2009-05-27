\name{dcemri}
\alias{dcemri.lm}
\title{Pharmacokinetic Models for Dynamic Contrast-Enhanced MRI Data}
\description{
  TBA
}
\usage{
dcemri.lm(conc, time, img.mask, model="extended", aif="tofts.kermode",
          nprint=0, user=NULL, ...)
}
\arguments{
  \item{conc}{}
  \item{time}{}
  \item{img.mask}{}
  \item{model}{}
  \item{aif}{}
  \item{nprint}{}
  \item{user}{}
  \item{...}{}
}
\value{
}
\details{
}
\references{
}
\seealso{
  dce.bayes, dce.spline
}
\examples{

data("buckley")
img <- t(breast$data)[,(1:60)*5]
time <- seq(5,300,5)/60
mask <- rep(TRUE,13)
fit <- dcemri.lm(img,time,mask)
fit$ktrans
}
\author{Brandon Whitcher, Volker Schmid}
\keyword{}
