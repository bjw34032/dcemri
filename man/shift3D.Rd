\name{shift3D}
\alias{shift3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
shift3D(A, s, type, fill = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ ~~Describe \code{A} here~~ }
  \item{s}{ ~~Describe \code{s} here~~ }
  \item{type}{ ~~Describe \code{type} here~~ }
  \item{fill}{ ~~Describe \code{fill} here~~ }
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
function(A, s, type, fill=0) {
  X <- nrow(A)
  Y <- ncol(A)
  Z <- nsli(A)
  if (s != 0) {
    if (type == "LR") {
      if (s > 0) {
        A <- A[c((X-s+1):X,1:(X-s)),,]
        A[1:s,,] <- fill
      } else {
        A <- A[c((abs(s)+1):X,1:abs(s)),,]
        A[(X-abs(s)+1):X,,] <- fill
      }
    } else {
      if (type == "AP") {
        if (s > 0) {
          A <- A[,c((Y-s+1):Y,1:(Y-s)),]
          A[,1:s,] <- fill
        } else {
          A <- A[,c((abs(s)+1):Y,1:abs(s)),]
          A[,(Y-abs(s)+1):Y,] <- fill
        }
      } else {
        if (type == "SI") {
          if (s > 0) {
            A <- A[,,c((Z-s+1):Z,1:(Z-s))]
            A[,,1:s] <- 0
          } else {
            A <- A[,,c((abs(s)+1):Z,1:abs(s))]
            A[,,(Z-abs(s)+1):Z] <- fill
          }
        } else {
          stop("Type of translation not recognized!")
        }
      }
    }
  }
  return(A)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
