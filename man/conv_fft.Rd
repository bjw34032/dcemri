\name{conv.fft}
\alias{conv.fft}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
conv.fft(A, B, C)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ ~~Describe \code{A} here~~ }
  \item{B}{ ~~Describe \code{B} here~~ }
  \item{C}{ ~~Describe \code{C} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(A, B, C) {
  if (length(dim(A)) == 3) {
    if (length(dim(A)) == length(dim(B)) && length(dim(A)) == length(C)) {
      X <- nrow(A)
      Y <- ncol(A)
      Z <- nsli(A)
      out <- Re(fft(fft(A) * Conj(fft(B)), inv=TRUE))[X:1,Y:1,Z:1] / (X*Y*Z)
      out <- out[c((X-C[1]+1):X, 1:(X-C[1])),,]
      out <- out[,c((Y-C[2]+1):Y, 1:(Y-C[2])),]
      out <- out[,,c((Z-C[3]+1):Z, 1:(Z-C[3]))]
      return(out)
    } else {
      stop("Objects are not all the same dimension!")
    }
  } else {
    stop("Only three dimensional objects!")
  }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
