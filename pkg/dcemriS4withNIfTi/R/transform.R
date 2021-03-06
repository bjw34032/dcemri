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
## $Id$
##

############################################################################
## performPermutation
############################################################################

performPermutation <- function(T, real.dimensions, data, verbose=FALSE) {
  workingdims <- (
                  function(r) {
                    lr <- length(r)
                    if (lr <= 5)
                      c(r, rep(1, 5 - lr))
                    else
                      stop("array has dim > 5")
                  }
                  )(real.dimensions) # An anonymous function

  if (sum(T != 0) == 3 && det(T) != 0) {
    ## Now ensure T is descaled and work out the permutation
    trans <- sign(T)
    perms <- abs(trans %*% 1:3)
    if (length(perms) != length(workingdims))
      perms <- (c(perms, (length(perms)+1):length(workingdims)))

    reverselist <- c(trans %*% rep(1,3) < 0,
                     rep(FALSE, length(workingdims)-3))

    if (any(reverselist[2:length(reverselist)]) || 
	any(perms != 1:length(perms))) {

      if (verbose) cat("need to permute", fill=TRUE)
      translatedData <- array(data, workingdims)
      ## Now if we have to do a permutation or reverse more than the first
      ## index we will be going slow anyway, so...
      prs <- (
              function(reverse, dims) {
                function(x) { 
                  if (reverse[x])
                    rev(1:dims[x]) 
                  else
                    1:dims[x] 
                }
              }
              )(reverselist, workingdims) # An anonymous function
      translatedData <- translatedData[prs(1), prs(2), prs(3), prs(4),
                                       prs(5), drop=FALSE]
      return(array(aperm(translatedData, perms), dim=real.dimensions))
    } else {
      if (reverselist[1]) {
        ## We just need to reverse the first index.
        return(array(array(data, dim=workingdims)[workingdims[1]:1,,,,],
                     dim=real.dimensions))
      } else {
        return(array(data, dim=real.dimensions))
      }
    }
  } else {
    stop("Transformation is not simple, cannot reorient!")
  }
}

############################################################################
## reorient
############################################################################

reorient <- function(nim, data, verbose=FALSE, invert=FALSE) {
  ## from nifti1.h there are three different methods of orienting the
  ## i,j,k data into x,y,z space.

  ## Method 1 is the default for ANALYZE 7.5 files, and will be dealt
  ## with last

  ## This function will try to reorient the data into an i, j, k space
  ## where increasing (+) (i,j,k) is correlated with
  ## (LEFT,ANTERIOR,SUPERIOR)
  real.dimensions <- nim@"dim_"[2:(1+nim@"dim_"[1])]

  if (nim@"qform_code" > 0) {
    ## Method 2. 
    ## This method determines [x] coordinates by the pixdim[] scales
    ## (the values of which we don't care about just the signs):
    S <- diag(nim@"pixdim"[2:4])
    ## A scaling factor, (because quaternions have to have det(R)=1):
    scalingFactor <- nim@"pixdim"[1]
    ## which is either = 1 or -1
    if (abs(scalingFactor) != 1) { 
      if (verbose)
	cat("ScalingFactor (nim@\"pixdim\"[1]) <-", scalingFactor,
	    "!= -1 or 1. Defaulting to 1", fill=TRUE)
      scalingFactor <- 1 
    }

    ## and applied to S[3,3]
    S[3,3] <- S[3,3] * scalingFactor

    ## a quaternion rotation matrix, R:
    R <- quaternion2rotation(nim@"quatern_b", nim@"quatern_c", nim@"quatern_d")
    ## And a shift, which we'll have to update if we do any reversals
    shift <- c(nim@"qoffset_x", nim@"qoffset_y", nim@"qoffset_z")
    ## X is determined by X <- R %*% S %*% (I - 1) + shift
    ## Now we can reorient the data only if the X axes are aligned with
    ## the I axes i.e. there are only 3 non-zero values in the matrix RS 
    RS <- R %*% S

    ## Now descale RS and work out the permutation
    trans <- sign(round(RS))
    ## We will need to do something with the trans later...
    trans[1,1] <- -1 * trans[1,1]

    if (invert)
      trans <- qr.solve(trans)

    return(performPermutation(trans, real.dimensions, data))
  } 
  if (nim@"sform_code" > 0) {
    ## [x] is given by a general affine transformation from [i]
    ##
    ## [x] <- A%*%(i-1) + shift
    ## 
    ## where A <- nim@srow_[x,y,z][1:3] and shift <- srow_x,y,z[4]
    S <- array(dim=c(3,4))
    S[1,] <- nim@"srow_x"
    S[2,] <- nim@"srow_y"
    S[3,] <- nim@"srow_z"
    shift <- S[,4]
    A <- S[,1:3]

    trans <- sign(round(A))
    trans[1,1] <- -1 * trans[1,1]

    if (invert)
      trans <- qr.solve(trans)

    return(performPermutation(trans, real.dimensions, data))
  }

  ## Finally Method 1.
  scaling <- diag(nim@"pixdim"[2:4])
  trans <- sign(scaling)
  ## Method 1 by default has +x going LEFT so no sign-change for trans[1,1]
  if (invert)
    trans <- qr.solve(trans)

  return(performPermutation(trans, real.dimensions, data))
}

############################################################################
## inverseReorient
############################################################################

inverseReorient <- function(nim, verbose=FALSE)
  reorient(nim, nim@.Data, verbose=verbose, invert=TRUE)

############################################################################
## integerTranslation
############################################################################

integerTranslation <- function(nim, data, verbose=FALSE) {
  ## 3D IMAGE (VOLUME) ORIENTATION AND LOCATION IN SPACE:
  ## There are 3 different methods by which continuous coordinates can
  ## attached to voxels.  The discussion below emphasizes 3D volumes,
  ## and the continuous coordinates are referred to as (x,y,z).  The
  ## voxel index coordinates (i.e., the array indexes) are referred to
  ## as (i,j,k), with valid ranges:
  ##   i = 0 .. dim[1]-1
  ##   j = 0 .. dim[2]-1  (if dim[0] >= 2)
  ##   k = 0 .. dim[3]-1  (if dim[0] >= 3)
  ## The (x,y,z) coordinates refer to the CENTER of a voxel.  In
  ## methods 2 and 3, the (x,y,z) axes refer to a subject-based
  ## coordinate system, with
  ##   +x = Right  +y = Anterior  +z = Superior.
  ## This is a right-handed coordinate system.  However, the exact
  ## direction these axes point with respect to the subject depends on
  ## qform_code (Method 2) and sform_code (Method 3).

  dims <- 2:(1+nim@"dim_"[1])

  if (nim@"qform_code" <= 0 && nim@"sform_code" <= 0 ) {
    if (verbose) cat("  dims =", nim@"dim_"[dims], fill=TRUE)
    return(array(data, nim@"dim_"[dims]))
  } else {
    i <- 0:(nim@"dim_"[2]-1)
    j <- 0:(nim@"dim_"[3]-1)
    k <- 0:(nim@"dim_"[4]-1)
    ijk <- cbind(rep(i, nim@"dim_"[3] * nim@"dim_"[4]),
                 rep(rep(j, each=nim@"dim_"[2]), nim@"dim_"[4]),
                 rep(k, each=nim@"dim_"[2] * nim@"dim_"[3]))
    index.ijk <- (ijk[,1] +
                  ijk[,2] * nim@"dim_"[2] +
                  ijk[,3] * nim@"dim_"[2] * nim@"dim_"[3])
    ## check for qform codes
    if (nim@"qform_code" > 0) {
      if (verbose) cat("  NIfTI-1: qform_code > 0", fill=TRUE)
      qfac <- nim@"pixdim"[1]
      R <- quaternion2rotation(nim@"quatern_b",
                               nim@"quatern_c",
                               nim@"quatern_d")
      ## HACK!!! To ensure matrix is integer-valued
      R <- round(R)
      qoffset <- c(nim@"qoffset_x", nim@"qoffset_y", nim@"qoffset_z")
      if (qfac < 0)
        R[3,3] <- -R[3,3]
      if (all(abs(R) == diag(3))) {
        ## HACK!!! Multiply x-dimension for proper orientation in R
        R[1,] <- -R[1,]
        xyz <-
          t(sweep(R %*% t(sweep(ijk, 2, as.array(nim@"pixdim"[2:4]), "*")),
                  1, as.array(qoffset), "+"))
        index.xyz <- (xyz[,1] +
                      xyz[,2] * nim@"dim_"[2] +
                      xyz[,3] * nim@"dim_"[2] * nim@"dim_"[3])      
        if (verbose) cat("  dims =", nim@"dim_"[dims], fill=TRUE)
        return(array(data[order(index.xyz)], nim@"dim_"[dims]))
      } else {
        stop("-- rotation matrix is NOT (approximately) diagonal with +/- 1s --")
      }
    }
    ## check for sform codes
    if (nim@"sform_code" > 0) {
      if (verbose) cat("  NIfTI-1: sform_code > 0", fill=TRUE)
      xyz <- matrix(0, length(data), 3)
      xyz[,1] <- (nim@"srow_x"[1] * ijk[,1] + nim@"srow_x"[2] * ijk[,2] +
                  nim@"srow_x"[3] * ijk[,3] + nim@"srow_x"[4])
      ## HACK!!! Multiply x-dimension for proper orientation in R
      xyz[,1] <- -xyz[,1]
      xyz[,2] <- (nim@"srow_y"[1] * ijk[,1] + nim@"srow_y"[2] * ijk[,2] +
                  nim@"srow_y"[3] * ijk[,3] + nim@"srow_y"[4])
      xyz[,3] <- (nim@"srow_z"[1] * ijk[,1] + nim@"srow_z"[2] * ijk[,2] +
                  nim@"srow_z"[3] * ijk[,3] + nim@"srow_z"[4])
      index.xyz <- (xyz[,1] +
                    xyz[,2] * nim@"dim_"[2] +
                    xyz[,3] * nim@"dim_"[2] * nim@"dim_"[3])
      if (verbose) cat("  dims =", nim@"dim_"[dims], fill=TRUE)
      return(array(data[order(index.xyz)], nim@"dim_"[dims]))
    }
  }
}

## write commands
##    
##    writeBin(as.vector(nim), fid, size=nim@"bitpix"/8)
##    writeBin(as.vector(nim@.Data[order(index.xyz)]), fid,
##                 size=nim@"bitpix"/8)
##      writeBin(as.vector(nim@.Data[order(index.xyz)]), fid,
##               size=nim@"bitpix"/8)

############################################################################
## invertIntegerTranslation
############################################################################

invertIntegerTranslation <- function(nim, verbose=FALSE) {
  dims <- 2:(1+nim@"dim_"[1])
  if (nim@"qform_code" <= 0 && nim@"sform_code" <= 0) {
    if (verbose)
      cat("  dims =", nim@"dim_"[dims], fill=TRUE)
    return(nim@.Data)
  } else {
    i <- 0:(nim@"dim_"[2]-1)
    j <- 0:(nim@"dim_"[3]-1)
    k <- 0:(nim@"dim_"[4]-1)
    ijk <- cbind(rep(i, nim@"dim_"[3] * nim@"dim_"[4]),
                 rep(rep(j, each=nim@"dim_"[2]), nim@"dim_"[4]),
                 rep(k, each=nim@"dim_"[2] * nim@"dim_"[3]))
    index.ijk <- (ijk[,1] +
                  ijk[,2] * nim@"dim_"[2] +
                  ijk[,3] * nim@"dim_"[2] * nim@"dim_"[3])
    ## check for qform codes
    if (nim@"qform_code" > 0) {
      if (verbose) {
        cat("  NIfTI-1: qform_code > 0", fill=TRUE)
        cat("  dims =", nim@"dim_"[dims], fill=TRUE)
      }
      qfac <- nim@"pixdim"[1]
      R <- quaternion2rotation(nim@"quatern_b",
                               nim@"quatern_c",
                               nim@"quatern_d")
      qoffset <- c(nim@"qoffset_x", nim@"qoffset_y", nim@"qoffset_z")
      if (qfac < 0)
        R[3,3] <- -R[3,3]
      if (all(abs(R) == diag(3))) {
        ## HACK!!! Multiply x-dimension for proper orientation in R
        R[1,] <- -R[1,]
        xyz <-
          t(sweep(R %*% t(sweep(ijk, 2, as.array(nim@"pixdim"[2:4]), "*")),
                  1, as.array(qoffset), "+"))
        index.xyz <- (xyz[,1] +
                      xyz[,2] * nim@"dim_"[2] +
                      xyz[,3] * nim@"dim_"[2] * nim@"dim_"[3])      
        if (verbose) cat("  dims =", nim@"dim_"[dims], fill=TRUE)
	return(nim@.Data[order(index.xyz)])
      } else {
        stop("-- rotation matrix is NOT diagonal with +/- 1s --")
      }
      ## stop("-- qform_code > 0 not implemented --")
    }
    ## check for sform codes
    if (nim@"sform_code" > 0) {
      if (verbose) {
        cat("  NIfTI-1: sform_code > 0", fill=TRUE)
        cat("  dims =", nim@"dim_"[dims], fill=TRUE)
      }
      xyz <- matrix(0, length(nim@.Data), 3)
      xyz[,1] <- (nim@"srow_x"[1] * ijk[,1] + nim@"srow_x"[2] * ijk[,2] +
                  nim@"srow_x"[3] * ijk[,3] + nim@"srow_x"[4])
      ## HACK!!! Multiply x-dimension for proper orientation in R
      xyz[,1] <- -xyz[,1]
      xyz[,2] <- (nim@"srow_y"[1] * ijk[,1] + nim@"srow_y"[2] * ijk[,2] +
                  nim@"srow_y"[3] * ijk[,3] + nim@"srow_y"[4])
      xyz[,3] <- (nim@"srow_z"[1] * ijk[,1] + nim@"srow_z"[2] * ijk[,2] +
                  nim@"srow_z"[3] * ijk[,3] + nim@"srow_z"[4])
      index.xyz <- (xyz[,1] +
                    xyz[,2] * nim@"dim_"[2] +
                    xyz[,3] * nim@"dim_"[2] * nim@"dim_"[3])
      return(nim@.Data[order(index.xyz)])
    }
  }
}

############################################################################
## translateCoordinate
############################################################################

translateCoordinate <- function(i, nim, verbose=FALSE) {
  ## 3D Image orientation and location in space (as per nifti1.h)
  ##
  ## There are three different methods by which xyzt(i,j,k) may be determined
  ## I will henceforth write x for realspace co-ord of i a voxel in the data
  ## array with indices i,j,k i.e. nim@.Data[i,j,k].
  ## 
  ## NB 1: nifti refers to voxel co-ords as:
  ## i = 0:dim[1] etc however, R indices are 1 based
  ##
  ## NB 2: x(i) refers to the centre of the voxel.
  ## 
  ## NB 3: Methods 2 & 3 have subject based co-ords with increasing x,y,z going
  ## Right, Anteriorly, Superiorly respectively
  ## 
  ## This is a right-handed coordinate system. However, the exact direction
  ## these axes point with respect to the subject depends on qform_code
  ## (Method 2) and sform_code (Method 3).
  ##
  ## More NBs:
  ## 1. By default the i index varies most rapidly, etc.
  ## 2. The ANALYZE 7.5 coordinate system is
  ##   +x = Left  +y = Anterior  +z = Superior
  ## (A left-handed co-ordinate system)
  ## 3. The three methods below give the locations of the voxel centres in the 
  ##   x,y,z system. In many cases programs will want to display the data on
  ##   other grids. In which case the program will be required to convert the
  ##   desired (x,y,z) values in to voxel values using the inverse
  ##   transformation.
  ## 4. Method 2 uses a factor qfac which is either -1 or 1. qfac is stored in
  ##   pixdim[0]. if pixdim[0]!= 1 or -1, which should not occur, we assume 1.
  ## 5. The units of the xyzt are set in xyzt_units field
  
  ## Method 2. when qform_code > 0, which should be the "normal" case
  if (nim@"qform_code" > 0) { 
    if (verbose) {
      cat("QForm_code <-", nim@"qform_code", ": Orientation by Method 2.",
	fill=TRUE)
    }
    ## The [x] coordinates are given by the pixdim[] scales:
    scaling <- diag(nim@"pixdim"[2:4])
    ##  a quaternion rotation matrix, R:
    R <- quaternion2rotation(nim@"quatern_b",
	nim@"quatern_c",
	nim@"quatern_d")
    ## A scaling factor, (because quaternions have to have det(R)=1):
    scalingFactor <- nim@"pixdim"[1]
    ## and a shift 
    shift <- c(nim@"qoffset_x", nim@"qoffset_y", nim@"qoffset_z")
    
    ## This method is intended to represent "scanner-anatomical"
    ## coordinates, which are often embedded in the image header, e.g. DICOM
    ## fields (0020,0032), (0020,0037), (0028,0030), and (0018,0050). These
    ## represent the nominal orientation and location of the data. This method
    ## can also be used to represent "aligned" coordinates, which would
    ## typically result from post-acquisition alignment of the volume to a
    ## "standard" orientation e.g.  the same subject on another day, or a rigid
    ## rotation to true anatomical orientation from the tilted position of the
    ## subject in the scanner.
    ## 
    ## [x] =  _R_%*%scaling%*%[i]+shift
    ##
    ## Where scaling[3] = scaling[3]*scalingFactor
    
    ## first enforce scalingFactor = 1 or -1
    if (scalingFactor != 1 && scalingFactor != -1) { 
      if (verbose) {
	cat("ScalingFactor (nim@\"pixdim\"[1]) <-", scalingFactor,
	    "!= -1 or 1. Defaulting to 1", fill=TRUE)
      }
      scalingFactor <- 1 
    }

    scaling[3,3] <- scaling[3,3] * scalingFactor

    ## Now we can reorient the data only if the X axes are aligned with the I
    ## axes i.e. there are only 3 non-zero values in the matrix RS 
    RS <- R %*% S
    return(RS %*% (i - 1) + shift)

  }

  ## Method 3. when sform_code > 0
  if (nim@"sform_code" > 0) {
    if (verbose)
      cat("SForm_code <- ", nim@"sform_code", ": Orientation by Method 3.",
	  fill=TRUE)
    ## [x] is given by a general affine transformation from [i]
    ##
    ## [x] <- A%*%(i-1) + shift
    ## 
    ## where A <- nim@srow_[x,y,z][1:3] and shift <- srow_x,y,z[4]
    S <- array(dim=c(3,4))
    S[1,] <- nim@"srow_x"
    S[2,] <- nim@"srow_y"
    S[3,] <- nim@"srow_z"
    shift <- S[,4]
    A <- S[,1:3]

    return (A %*% (i - 1) + shift)
  }

  ## Method 1. The `old' way used only if "qform_code" is 0
  ## The co-ord mapping from [i] to [x] is the ANALYZE 7.5 way.
  ## A simple scaling relationship applies
  ##
  ## x <- pixdim[1:3] * i[1:3]
  if (verbose)
    cat("QForm_code and SForm_code unset: Orientation by Method 1.", fill=TRUE)

  ## nifti1.h pixdim[1] <- ourpixdim[2]
  scaling <- diag(nim@"pixdim"[2:4])

  ## Remember i. <- i - 1
  return(scaling %*% (i - 1))
  
}
