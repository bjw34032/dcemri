##
##
## Copyright (c) 2009, Brandon Whitcher and Volker Schmid
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
##

conv.fft <- function(A, B, C) {
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

find.center <- function(M) {
  X <- nrow(M)
  Y <- ncol(M)
  Z <- nsli(M)
  xx <- array(1:X, dim(M))
  yy <- array(rep(1:Y, each=X), dim(M))
  zz <- array(rep(1:Z, each=X*Y), dim(M))
  center <- c(mean(xx[M]), mean(yy[M]), mean(zz[M]))
  trunc(center)
}

shift3D <- function(A, s, type, fill=0) {
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

