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

#############################################################################
## setClass("anlz")
#############################################################################

setClass("anlz",
         representation("sizeof_hdr" = "numeric",
                        "db_type" = "character",
                        "db_name" = "character",
                        "extents" = "numeric",
                        "session_error" = "numeric",
                        "regular" = "character",
                        "hkey_un0" = "character",
                        "dim_" = "vector",
                        "vox_units" = "character",
                        "cal_units" = "character",
                        "unused1" = "numeric",
                        "datatype" = "numeric",
                        "bitpix" = "numeric",
                        "dim_un0" = "numeric",
                        "pixdim" = "vector",
                        "vox_offset" = "numeric",
                        "funused1" = "numeric",
                        "funused2" = "numeric",
                        "funused3" = "numeric",
                        "cal_max" = "numeric",
                        "cal_min" = "numeric",
                        "compressed" = "numeric",
                        "verified" = "numeric",
                        "glmax" = "numeric",
                        "glmin" = "numeric",
                        "descrip" = "character",
                        "aux_file" = "character",
                        "orient" = "character",
                        "originator" = "character",
                        "generated" = "character",
                        "scannum" = "character",
                        "patient_id" = "character",
                        "exp_date" = "character",
                        "exp_time" = "character",
                        "hist_un0" = "character",
                        "views" = "numeric",
                        "vols_added" = "numeric",
                        "start_field" = "numeric",
                        "field_skip" = "numeric",
                        "omax" = "numeric",
                        "omin" = "numeric",
                        "smax" = "numeric",
                        "smin" = "numeric"),
         prototype("sizeof_hdr" = 348,
                   "db_type" = "",
                   "db_name" = "",
                   "extents" = numeric(1),
                   "session_error" = numeric(1),
                   "regular" = "r",
                   "hkey_un0" = "",
                   "dim_" = numeric(8),
                   "vox_units" = "mm",
                   "cal_units" = "",
                   "unused1" = numeric(1),
                   "datatype" = 2,
                   "bitpix" = 8,
                   "dim_un0" = numeric(1),
                   "pixdim" = numeric(8),
                   "vox_offset" = numeric(1),
                   "funused1" = numeric(1),
                   "funused2" = numeric(1),
                   "funused3" = numeric(1),
                   "cal_max" = numeric(1),
                   "cal_min" = numeric(1),
                   "compressed" = numeric(1),
                   "verified" = numeric(1),
                   "glmax" = numeric(1),
                   "glmin" = numeric(1),
                   "descrip" = "",
                   "aux_file" = "",
                   "orient" = "0",
                   "originator" = "",
                   "generated" = "",
                   "scannum" = "",
                   "patient_id" = "",
                   "exp_date" = "",
                   "exp_time" = "",
                   "hist_un0" = "",
                   "views" = numeric(1),
                   "vols_added" = numeric(1),
                   "start_field" = numeric(1),
                   "field_skip" = numeric(1),
                   "omax" = numeric(1),
                   "omin" = numeric(1),
                   "smax" = numeric(1),
                   "smin" = numeric(1)),
         contains="array")

#############################################################################
## setMethod("show", "anlz")
#############################################################################

setMethod("show", "anlz", function(object) {
  cat("ANALYZE 7.5 format", fill=TRUE)
  cat("  Type            :", class(object), fill=TRUE)
  cat("  Data Type       : ", object@"datatype", " (",
      convert.datatype.anlz(object@"datatype"), ")", sep="", fill=TRUE)
  cat("  Bits per Pixel  :", object@"bitpix", fill=TRUE)
  cat("  Orient          : ", object@"orient", " (",
      convert.orient.anlz(object@"orient"), ")", sep="", fill=TRUE)
  cat("  Dimension       :",
      paste(object@"dim_"[2:(1+object@"dim_"[1])], collapse=" x "), fill=TRUE)
  cat("  Pixel Dimension :",
      paste(round(object@"pixdim"[2:(1+object@"dim_"[1])], 2),
            collapse=" x "), fill=TRUE)
  cat("  Voxel Units     :", object@"vox_units", fill=TRUE)
})

#############################################################################
## setValidity("anlz")
#############################################################################

setValidity("anlz", function(object) {
  retval <- NULL
  indices <- 2:(1+object@"dim_"[1])
  ## sizeof_hdr must be 348
  if (object@"sizeof_hdr" != 348)
    retval <- c(retval, "sizeof_hdr != 348")
  ## datatype needed to specify type of image data
  if (!object@"datatype" %in% c(0,2^(0:7), 255))
    retval <- c(retval, "datatype not recognized")
  ## bitpix should correspond correctly to datatype
  
  ## dim should be non-zero for dim[1] dimensions
  if (!all(as.logical(object@"dim_"[indices])))
    retval <- c(retval, "dim[1]/dim mismatch")
  ## number of data dimensions should match dim[1]
  if (length(indices) != length(dim(object@.Data)))
    retval <- c(retval, "dim[1]/img mismatch")
  ## pixdim[n] required when dim[n] is required
  if (!all(as.logical(object@"dim_"[indices]) &&
           as.logical(object@"pixdim"[indices])))
    retval <- c(retval, "dim/pixdim mismatch")
  ## data dimensions should match dim 
  if (!all(object@"dim_"[indices] == dim(object@.Data)))
    retval <- c(retval, "dim/img mismatch")
  if (is.null(retval)) return(TRUE)
  else return(retval)
})

#############################################################################
## anlz()
#############################################################################

anlz <- function(img=array(0, dim=rep(1,4)), dim, ...) {
  if (missing(dim)) {
    if (is.array(img))
      dim <- base::dim(img)
    else
      dim <- c(1, length(img))
  }
  ld <- length(dim)
  if (ld < 3)
    stop(sprintf("length(dim) must be at least 3 and is %d.", ld))
  
  x <- c(length(dim), dim[1], dim[2], dim[3],
         ifelse(is.na(dim[4]), 1, dim[4]), rep(1,3))
  y <- c(0.0, rep(1.0,length(dim)), rep(0.0,3))
  cal.max <- quantile(img, probs=0.95, na.rm=TRUE)
  cal.min <- quantile(img, probs=0.05, na.rm=TRUE)
  obj <- new("anlz", .Data=array(img, dim=dim), "dim_"=x, "pixdim"=y,
             "cal_max"=cal.max, "cal_min"=cal.min, ...)
  validObject(obj)
  return(obj)
}

#############################################################################
## anlz()
#############################################################################

is.anlz <- function(x) {
  if (!is(x, "anlz"))
    return(FALSE)
  else
    return (TRUE)
}

#############################################################################
## descrip() and descrip<-()
#############################################################################

if (!isGeneric("descrip")) {
  if (is.function("descrip"))
    setGeneric("descrip", descrip)
  else
    setGeneric("descrip", function(object) { standardGeneric("descrip") })
}
setMethod("descrip", "anlz", function(object) { object@descrip })
setGeneric("descrip<-", function(x, value) { standardGeneric("descrip<-") })
setReplaceMethod("descrip", "anlz",
                 function(x, value) { x@descrip <- value ; x })

#############################################################################
## aux.file() and aux.file<-()
#############################################################################

if (!isGeneric("aux.file")) {
  if (is.function("aux.file"))
    setGeneric("aux.file", aux_file)
  else
    setGeneric("aux.file", function(object) { standardGeneric("aux.file") })
}
setMethod("aux.file", "anlz", function(object) { object@aux_file })
setGeneric("aux.file<-", function(x, value) { standardGeneric("aux.file<-") })
setReplaceMethod("aux.file", "anlz",
                 function(x, value) { x@"aux_file" <- value ; x })
