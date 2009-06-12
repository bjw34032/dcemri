\name{find.center}
\alias{find.center}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find.center(M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{ ~~Describe \code{M} here~~ }
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
function(M) {
  X <- nrow(M)
  Y <- ncol(M)
  Z <- nsli(M)
  xx <- array(1:X, dim(M))
  yy <- array(rep(1:Y, each=X), dim(M))
  zz <- array(rep(1:Z, each=X*Y), dim(M))
  center <- c(mean(xx[M]), mean(yy[M]), mean(zz[M]))
  trunc(center)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
