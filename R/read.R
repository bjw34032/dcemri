read.hdr <- function(fname, verbose=FALSE, warn=-1) {
  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)
  ## Check if any file extensions are present
  ANALYZE <- ifelse(length(grep("hdr", fname)) != 0, TRUE, FALSE)
  NIFTI <- ifelse(length(grep("nii", fname)) != 0, TRUE, FALSE)
  GZ <- ifelse(length(grep("gz", fname)) != 0, TRUE, FALSE)

  if (GZ) {
    if (ANALYZE) {
      ls.hdr.gz <- system(paste("ls", fname), intern=TRUE)
      if (length(ls.hdr.gz) != 0) {
        if (verbose)
          cat(paste("  fname =", fname, "and file =", ls.hdr.gz), fill=TRUE)
        hdr <- read.analyze.hdr(sub(".hdr.gz", "", fname), gzipped=TRUE,
                                verbose=verbose, warn=warn)
        options(warn=oldwarn)
        return(hdr)
      }
    } else {
      if (NIFTI) {
        ls.nii.gz <- system(paste("ls", fname), intern=TRUE)
        if (length(ls.nii.gz) != 0) {
          if (verbose)
            cat(paste("  fname =", fname, "and file =", ls.nii.gz), fill=TRUE)
          hdr <- read.nifti.hdr(sub(".nii.gz", "", fname), gzipped=TRUE,
                                verbose=verbose, warn=warn)
          options(warn=oldwarn)
          return(hdr)
        }
      } else {
        options(warn=oldwarn)
        stop(paste(fname, "is not recognized."))
      }
    }    
  } else {
    if (ANALYZE) {
      ls.hdr <- system(paste("ls", fname), intern=TRUE)
      if (length(ls.hdr) != 0) {
        if (verbose)
          cat(paste("  fname =", fname, "and file =", ls.hdr), fill=TRUE)
        hdr <- read.analyze.hdr(sub(".hdr", "", fname), gzipped=FALSE,
                                verbose=verbose, warn=warn)
        options(warn=oldwarn)
        return(hdr)
      } else {
        fname.hdr.gz <- paste(fname, "gz", sep=".")
        ls.hdr.gz <- system(paste("ls", fname.hdr.gz), intern=TRUE)
        if (length(ls.hdr.gz) != 0) {
          if (verbose)
            cat(paste("  fname =", fname, "and file =", ls.hdr.gz), fill=TRUE)
          hdr <- read.analyze.hdr(sub(".hdr", "", fname), gzipped=TRUE,
                                  verbose=verbose, warn=warn)
          options(warn=oldwarn)
          return(hdr)
        } else {
          options(warn=oldwarn)
          stop(paste(fname, "is not recognized."))
        }
      }
    } else {
      if (NIFTI) {
        ls.nii <- system(paste("ls", fname), intern=TRUE)
        if (length(ls.nii) != 0) {
          if (verbose)
            cat(paste("  fname =", fname, "and file =", ls.nii), fill=TRUE)
          hdr <- read.nifti.hdr(sub(".nii", "", fname), gzipped=FALSE,
                                verbose=verbose, warn=warn)
          options(warn=oldwarn)
          return(hdr)
        } else {
          fname.nii.gz <- paste(fname, "gz", sep=".")
          ls.nii.gz <- system(paste("ls", fname.nii.gz), intern=TRUE)
          if (length(ls.nii.gz) != 0) {
            if (verbose)
              cat(paste("  fname =", fname, "and file =", ls.nii.gz),
                  fill=TRUE)
            hdr <- read.nifti.hdr(sub(".nii", "", fname), gzipped=TRUE,
                                  verbose=verbose, warn=warn)
            options(warn=oldwarn)
            return(hdr)
          } else {
            options(warn=oldwarn)
            stop(paste(fname, "is not recognized."))
          }
        }
      } else {
        fname.hdr <- paste(fname, "hdr", sep=".")
        ls.hdr <- system(paste("ls", fname.hdr), intern=TRUE)
        if (length(ls.hdr) != 0) {
          if (verbose)
            cat(paste("  fname =", fname, "and file =", ls.hdr), fill=TRUE)
          hdr <- read.analyze.hdr(fname, gzipped=FALSE, verbose=verbose,
                                  warn=warn)
          options(warn=oldwarn)
          return(hdr)
        } else {
          fname.hdr.gz <- paste(fname, "hdr.gz", sep=".")
          ls.hdr.gz <- system(paste("ls", fname.hdr.gz), intern=TRUE)
          if (length(ls.hdr.gz) != 0) {
            if (verbose)
              cat(paste("  fname =", fname, "and file =", ls.hdr.gz),
                  fill=TRUE)
            hdr <- read.analyze.hdr(fname, gzipped=TRUE, verbose=verbose,
                                    warn=warn)
            options(warn=oldwarn)
            return(hdr)
          } else {
            fname.nii <- paste(fname, "nii", sep=".")
            ls.nii <- system(paste("ls", fname.nii), intern=TRUE)
            if (length(ls.nii) != 0) {
              if (verbose)
                cat(paste("  fname =", fname, "and file =", ls.nii), fill=TRUE)
              hdr <- read.nifti.hdr(fname, gzipped=FALSE, verbose=verbose,
                                    warn=warn)
              options(warn=oldwarn)
              return(hdr)
            } else {
              fname.nii.gz <- paste(fname, "nii.gz", sep=".")
              ls.nii.gz <- system(paste("ls", fname.nii.gz), intern=TRUE)
              if (length(ls.nii.gz) != 0) {
                if (verbose)
                  cat(paste("  fname =", fname, "and file =", ls.nii.gz),
                      fill=TRUE)
                hdr <- read.nifti.hdr(fname, gzipped=TRUE, verbose=verbose,
                                      warn=warn)
                options(warn=oldwarn)
                return(hdr)
              } else {
                options(warn=oldwarn)
                stop(paste(fname, "is not recognized."))
              }
            }
          }
        }
      }
    }
  }
}

read.analyze.hdr <- function(fname, gzipped=TRUE, verbose=FALSE, warn=-1) {

  if (gzipped)
    fid <- gzfile(paste(fname, ".hdr.gz", sep=""), "rb")
  else
    fid <- file(paste(fname, ".hdr", sep=""), "rb")

  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)
  ## Test for endian properties
  endian <- .Platform$endian
  size.of.hdr <- readBin(fid, integer(), size=4)
  if (size.of.hdr != 348) {
    close(fid)
    endian <- "swap"
    if(gzipped)
      fid <- gzfile(paste(fname, ".hdr.gz", sep=""), "rb")
    else
      fid <- file(paste(fname, ".hdr", sep=""), "rb")
    size.of.hdr <- readBin(fid, integer(), size=4, endian=endian)
    if(verbose)
      cat("  ENDIAN = swap", fill=TRUE)
  }
  
  db.type <- rawToChar(readBin(fid, "raw", n=10))
  db.name <- rawToChar(readBin(fid, "raw", n=18))
  
  extents <- readBin(fid, integer(), size=4, endian=endian)
  session.error <- readBin(fid, integer(), size=2, endian=endian)
  
  regular <- rawToChar(readBin(fid, "raw", n=1))
  hkey.un0 <- rawToChar(readBin(fid, "raw", n=1))
  
  dim <- readBin(fid, integer(), 8, size=2, endian=endian)
  
  vox.units <- rawToChar(readBin(fid, "raw", n=4))
  cal.units <- rawToChar(readBin(fid, "raw", n=8))
  
  skip <- readBin(fid, integer(), size=2, endian=endian)
  datatype <- readBin(fid, integer(), size=2, endian=endian)
  bitpix <- readBin(fid, integer(), size=2, endian=endian)
  dim.un0 <- readBin(fid, integer(), size=2, endian=endian)
  pixdim <- readBin(fid, numeric(), 8, size=4, endian=endian)
  vox.offset <- readBin(fid, numeric(), size=4, endian=endian)

  skip <- readBin(fid, numeric(), n=3, size=4, endian=endian)

  cal.max <- readBin(fid, numeric(), size=4, endian=endian)
  cal.min <- readBin(fid, numeric(), size=4, endian=endian)
  compressed <- readBin(fid, numeric(), size=4, endian=endian)
  verified <- readBin(fid, numeric(), size=4, endian=endian)
  glmax <- readBin(fid, integer(), size=4, endian=endian)
  glmin <- readBin(fid, integer(), size=4, endian=endian)
 
  descrip <- rawToChar(readBin(fid, "raw", n=80))
  aux.file <- rawToChar(readBin(fid, "raw", n=24))
  ## This problem appeared in R-2.2.1 with error message:
  ## Error in readChar(fid, n = 1) : invalid UTF-8 input in readChar()
  ## So... readChar() has been replaced with readBin(..., raw(), ...)
  ## as suggested by Duncan Murdoch.
  ## 27 Apr 06
  ## Must use --disable-mbcs when compiling R to turn off UTF-8
  ## 31 May 06
  ## I am replacing all readChar() calls with
  ##   rawToChar(readBin(..., "raw", ...))
  ## 13 Feb 09 
  orient <- rawToChar(readBin(fid, "raw", n=1))
  originator <- rawToChar(readBin(fid, "raw", n=10))
  generated <- rawToChar(readBin(fid, "raw", n=10))
  scannum <- rawToChar(readBin(fid, "raw", n=10))
  patient.id <- rawToChar(readBin(fid, "raw", n=10))
  exp.date <- rawToChar(readBin(fid, "raw", n=10))
  exp.time <- rawToChar(readBin(fid, "raw", n=10))
  hist.un0 <- rawToChar(readBin(fid, "raw", n=3))
  
  views <- readBin(fid, integer(), size=4, endian=endian)
  vols.added <- readBin(fid, integer(), size=4, endian=endian)
  start.field <- readBin(fid, integer(), size=4, endian=endian)
  field.skip <- readBin(fid, integer(), size=4, endian=endian)
  omax <- readBin(fid, integer(), size=4, endian=endian)
  omin <- readBin(fid, integer(), size=4, endian=endian)
  smax <- readBin(fid, integer(), size=4, endian=endian)
  smin <- readBin(fid, integer(), size=4, endian=endian)
  
  close(fid)
  
  output <- list(size.of.hdr, endian, extents, session.error, dim, datatype,
                 bitpix, dim.un0, pixdim, vox.offset, cal.max, cal.min,
                 compressed, verified, glmax, glmin, views, vols.added,
                 start.field, field.skip, omax, omin, smax, smin, descrip,
                 db.type, db.name, regular, hkey.un0, vox.units,
                 cal.units, aux.file, orient, originator, generated,
                 scannum, patient.id, exp.date, exp.time, hist.un0)

  names(output) <-
    c("size.of.hdr", "endian", "extents", "session.error", "dim", "datatype",
      "bitpix", "dim.un0", "pixdim", "vox.offset", "cal.max", "cal.min",
      "compressed", "verified", "glmax", "glmin", "views", "vols.added",
      "start.field", "field.skip", "omax", "omin", "smax", "smin",
      "descrip",  "db.type", "db.name", "regular", "hkey.un0", "vox.units",
      "cal.units", "aux.file","orient","originator", "generated",
      "scannum", "patient.id", "exp.date", "exp.time", "hist.un0")

  ## Warnings?
  options(warn=oldwarn)

  return(output)
}

make.hdr <- function(X, Y, Z, T, datatype, min=0, max=0) {
  if(X <=0 && Y <= 0 && Z <= 0 && T <= 0)
    stop("")
  if(!(datatype %in% c("UNKNOWN","BINARY","CHAR","SHORT","INT",
                       "FLOAT","COMPLEX","DOUBLE","RGB")))
    stop(paste("Unrecognized Data Type:", datatype))
  
  size.of.hdr <- as.integer(348)
  
  db.type <- rep("", 10) # readChar(fid, n=10)
  db.name <- rep("", 18) # readChar(fid, n=18)
  
  extents <- integer(1) # readBin(fid, integer(), size=4, endian=endian)
  session.error <- integer(1) # readBin(fid, integer(), size=2, endian=endian)
  
  regular <- "r" # readChar(fid, n=1)
  hkey.un0 <- "" # readChar(fid, n=1)

  dim <- integer(8)
  dim[1] = 4 # all Analyze images are taken as 4 dimensional
  dim[2] = as.integer(X) # slice width in pixels
  dim[3] = as.integer(Y) # slice height in pixels
  dim[4] = as.integer(Z) # volume depth in slices
  dim[5] = as.integer(T) # number of volumes per file

  vox.units <- rep("", 4) # readChar(fid, n=4)
  cal.units <- rep("", 8) # readChar(fid, n=8)

  datatype <- switch(datatype,
                     "UNKNOWN" = 0,
                     "BINARY" = 1,
                     "CHAR" = 2,
                     "SHORT" = 4,
                     "INT" = 8,
                     "FLOAT" = 16,
                     "COMPLEX" = 32,
                     "DOUBLE" = 64,
                     "RGB" = 128) # readBin(fid, integer(), size=2, endian=endian)
  bitpix <- switch(datatype,
                     "UNKNOWN" = 0,
                     "BINARY" = 1,
                     "CHAR" = 8,
                     "SHORT" = 16,
                     "INT" = 32,
                     "FLOAT" = 32,
                     "COMPLEX" = 64,
                     "DOUBLE" = 64,
                     "RGB" = 24) # readBin(fid, integer(), size=2, endian=endian)
  dim.un0 <- integer(1) # readBin(fid, integer(), size=2, endian=endian)
  pixdim <- numeric(8) # readBin(fid, numeric(), 8, size=4, endian=endian)
  vox.offset <- 0.0 # readBin(fid, numeric(), size=4, endian=endian)

  cal.max <- 0.0 # readBin(fid, numeric(), size=4, endian=endian)
  cal.min <- 0.0 # readBin(fid, numeric(), size=4, endian=endian)
  compressed <- 0.0 # readBin(fid, numeric(), size=4, endian=endian)
  verified <- 0.0 # readBin(fid, numeric(), size=4, endian=endian)
  glmax <- max # readBin(fid, integer(), size=4, endian=endian)
  glmin <- min # readBin(fid, integer(), size=4, endian=endian)
 
  descrip <- rep("", 80) # readChar(fid, n=80)
  aux.file <- rep("", 24) # readChar(fid, n=24)

  orient <- "" # readChar(fid, n=1)
  originator <- rep("", 10) # readChar(fid, n=10)
  generated <- rep("", 10) # readChar(fid, n=10)
  scannum <- rep("", 10) # readChar(fid, n=10)
  patient.id <- rep("", 10) # readChar(fid, n=10)
  exp.date <- rep("", 10) # readChar(fid, n=10)
  exp.time <- rep("", 10) # readChar(fid, n=10)
  hist.un0 <- rep("", 3) # readChar(fid, n=3)
  
  views <- integer(1) # readBin(fid, integer(), size=4, endian=endian)
  vols.added <- integer(1) # readBin(fid, integer(), size=4, endian=endian)
  start.field <- integer(1) # readBin(fid, integer(), size=4, endian=endian)
  field.skip <- integer(1) # readBin(fid, integer(), size=4, endian=endian)
  omax <- integer(1) # readBin(fid, integer(), size=4, endian=endian)
  omin <- integer(1) # readBin(fid, integer(), size=4, endian=endian)
  smax <- integer(1) # readBin(fid, integer(), size=4, endian=endian)
  smin <- integer(1) # readBin(fid, integer(), size=4, endian=endian)
}

read.img <- function(fname, verbose=FALSE, warn=-1) {
  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)
  ## Check if any file extensions are present
  ANALYZE <- ifelse(length(grep("img", fname)) != 0, TRUE, FALSE)
  NIFTI <- ifelse(length(grep("nii", fname)) != 0, TRUE, FALSE)
  GZ <- ifelse(length(grep("gz", fname)) != 0, TRUE, FALSE)

  if (GZ) {
    if (ANALYZE) {
      ls.img.gz <- system(paste("ls", fname), intern=TRUE)
      if (length(ls.img.gz) != 0) {
        if (verbose)
          cat(paste("  fname =", fname, "and file =", ls.img.gz), fill=TRUE)
        img <- read.analyze.img(sub(".img.gz", "", fname), gzipped=TRUE,
                                verbose=verbose, warn=warn)
        options(warn=oldwarn)
        return(img)
      }
    } else {
      if (NIFTI) {
        ls.nii.gz <- system(paste("ls", fname), intern=TRUE)
        if (length(ls.nii.gz) != 0) {
          if (verbose)
            cat(paste("  fname =", fname, "and file =", ls.nii.gz), fill=TRUE)
          img <- read.nifti.img(sub(".nii.gz", "", fname), gzipped=TRUE,
                                verbose=verbose, warn=warn)
          options(warn=oldwarn)
          return(img)
        }
      } else {
        options(warn=oldwarn)
        stop(paste(fname, "is not recognized."))
      }
    }    
  } else {
    if (ANALYZE) {
      ls.img <- system(paste("ls", fname), intern=TRUE)
      if (length(ls.img) != 0) {
        if (verbose)
          cat(paste("  fname =", fname, "and file =", ls.img), fill=TRUE)
        img <- read.analyze.img(sub(".img", "", fname), gzipped=FALSE,
                                verbose=verbose, warn=warn)
        options(warn=oldwarn)
        return(img)
      } else {
        fname.img.gz <- paste(fname, "gz", sep=".")
        ls.img.gz <- system(paste("ls", fname.img.gz), intern=TRUE)
        if (length(ls.img.gz) != 0) {
          if (verbose)
            cat(paste("  fname =", fname, "and file =", ls.img.gz), fill=TRUE)
          img <- read.analyze.img(sub(".img", "", fname), gzipped=TRUE,
                                  verbose=verbose, warn=warn)
          options(warn=oldwarn)
          return(img)
        } else {
          options(warn=oldwarn)
          stop(paste(fname, "is not recognized."))
        }
      }
    } else {
      if (NIFTI) {
        ls.nii <- system(paste("ls", fname), intern=TRUE)
        if (length(ls.nii) != 0) {
          if (verbose)
            cat(paste("  fname =", fname, "and file =", ls.nii), fill=TRUE)
          img <- read.nifti.img(sub(".nii", "", fname), gzipped=FALSE,
                                verbose=verbose, warn=warn)
          options(warn=oldwarn)
          return(img)
        } else {
          fname.nii.gz <- paste(fname, "gz", sep=".")
          ls.nii.gz <- system(paste("ls", fname.nii.gz), intern=TRUE)
          if (length(ls.nii.gz) != 0) {
            if (verbose)
              cat(paste("  fname =", fname, "and file =", ls.nii.gz),
                  fill=TRUE)
            img <- read.nifti.img(sub(".nii", "", fname), gzipped=TRUE,
                                  verbose=verbose, warn=warn)
            options(warn=oldwarn)
            return(img)
          } else {
            options(warn=oldwarn)
            stop(paste(fname, "is not recognized."))
          }
        }
      } else {
        fname.img <- paste(fname, "img", sep=".")
        ls.img <- system(paste("ls", fname.img), intern=TRUE)
        if (length(ls.img) != 0) {
          if (verbose)
            cat(paste("  fname =", fname, "and file =", ls.img), fill=TRUE)
          img <- read.analyze.img(fname, gzipped=FALSE, verbose=verbose,
                                  warn=warn)
          options(warn=oldwarn)
          return(img)
        } else {
          fname.img.gz <- paste(fname, "img.gz", sep=".")
          ls.img.gz <- system(paste("ls", fname.img.gz), intern=TRUE)
          if (length(ls.img.gz) != 0) {
            if (verbose)
              cat(paste("  fname =", fname, "and file =", ls.img.gz),
                  fill=TRUE)
            img <- read.analyze.img(fname, gzipped=TRUE, verbose=verbose,
                                    warn=warn)
            options(warn=oldwarn)
            return(img)
          } else {
            fname.nii <- paste(fname, "nii", sep=".")
            ls.nii <- system(paste("ls", fname.nii), intern=TRUE)
            if (length(ls.nii) != 0) {
              if (verbose)
                cat(paste("  fname =", fname, "and file =", ls.nii), fill=TRUE)
              img <- read.nifti.img(fname, gzipped=FALSE, verbose=verbose,
                                    warn=warn)
              options(warn=oldwarn)
              return(img)
            } else {
              fname.nii.gz <- paste(fname, "nii.gz", sep=".")
              ls.nii.gz <- system(paste("ls", fname.nii.gz), intern=TRUE)
              if (length(ls.nii.gz) != 0) {
                if (verbose)
                  cat(paste("  fname =", fname, "and file =", ls.nii.gz),
                      fill=TRUE)
                img <- read.nifti.img(fname, gzipped=TRUE, verbose=verbose,
                                      warn=warn)
                options(warn=oldwarn)
                return(img)
              } else {
                options(warn=oldwarn)
                stop(paste(fname, "is not recognized."))
              }
            }
          }
        }
      }
    }
  }
}

read.analyze.img <- function(fname, gzipped=TRUE, signed=FALSE, verbose=FALSE,
                             warn=-1) {
  ##
  ## Datatype Table
  ## -----------------------------------------
  ##   0 = DT_NONE or DT_UNKOWN
  ##   1 = DT_BINARY (1 bit per voxel)
  ##   2 = DT_UNSIGNED_CHAR (8 bits per voxel)
  ##   4 = DT_SIGNED_SHORT (16 bits per voxel)
  ##   8 = DT_SINGED_INT (32 bits per voxel)
  ##  16 = DT_FLOAT (32 bits per voxel)
  ##  32 = DT_COMPLEX (2 x 32 bit single)
  ##  64 = DT_DOUBLE (64 bits per voxel)
  ## 128 = DT_RGB
  ## 255 = DT_ALL
  ##

  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)
  ## Read in data from the header file
  hdr <- read.analyze.hdr(fname, gzipped, verbose, warn)
  n.elements <- prod(hdr$dim[2:5])
  size <- hdr$bitpix / 8
  what <- integer()
  switch(as.character(hdr$datatype),
         "2" = { signed <- FALSE },
         "4" = { signed <- TRUE  },
         "8" = { signed <- TRUE  },
         "16" = { what <- numeric() ; signed <- TRUE },
         "32" = {
           signed <- TRUE
           n.elements <- prod(hdr$dim[2:5]) * 2
           size <- hdr$bitpix / 8 / 2
           what <- numeric()
         },
         "64" = { what <- double() ; signed <- TRUE },
         stop(paste("Data type ", hdr$datatype, "unsupported in ", fname,
                    ".img", sep=""))
         )
  ## Read in data from the image file
  if(gzipped)
    fid <- gzfile(paste(fname, ".img.gz", sep=""), "rb")
  else
    fid <- file(paste(fname, ".img", sep=""), "rb")
  image.data <- readBin(fid, what = what, n = n.elements, size = size,
                        signed = signed, endian = hdr$endian)
  close(fid)
  ## Place vector into four-dimensional array
  if(hdr$datatype != 32)
    image.data <- array(image.data, hdr$dim[2:5])
  else {
    odd <- seq(1, n.elements, by=2)
    even <- seq(2, n.elements, by=2)
    image.vec <- complex(real=image.data[odd], imag=image.data[even])
    image.data <- array(image.vec, hdr$dim[2:5])
  }
  ## Warnings?
  options(warn=oldwarn)

  return(image.data)
}

read.nifti.hdr <- function(fname, onefile=TRUE, gzipped=TRUE,
                           verbose=FALSE, warn=-1) {

  ## Convert codes to names
  convert.datatype <- function(dt) {
    ## defgroup NIFTI1_DATATYPE_ALIASES
    ## brief aliases for the nifti1 datatype codes
    switch(as.character(dt),
           "2" = "UINT8",
           "4" = "INT16",
           "8" = "INT32",
           "16" = "FLOAT32",
           "32" = "COMPLEX64",
           "64" = "FLOAT64",
           "128" = "RGB24",
           "256" = "INT8",
           "512" = "UINT16",
           "768" = "UINT32",
           "1024" = "INT64",
           "1280" = "UINT64",
           "1536" = "FLOAT128",
           "1792" = "COMPLEX128",
           "2048" = "COMPLEX256")
  }
  convert.intent <- function(ic) {
    ## defgroup NIFTI1_INTENT_CODES
    ##-------- These codes are for probability distributions -----------
    ## Most distributions have a number of parameters, below denoted
    ## by p1, p2, and p3, and stored in
    ##  - intent_p1, intent_p2, intent_p3 if dataset doesn't have 5th
    ##    dimension
    ##  - image data array                if dataset does have 5th
    ##    dimension
    ##
    ## Functions to compute with many of the distributions below can
    ## be found in the CDF library from U Texas.
    ##
    ## Formulas for and discussions of these distributions can be
    ## found in the following books:
    ##
    ##  [U] Univariate Discrete Distributions,
    ##      NL Johnson, S Kotz, AW Kemp.
    ##
    ##  [C1] Continuous Univariate Distributions, vol. 1,
    ##       NL Johnson, S Kotz, N Balakrishnan.
    ##
    ##  [C2] Continuous Univariate Distributions, vol. 2,
    ##       NL Johnson, S Kotz, N Balakrishnan.
    ##
    switch(as.character(ic),
           "0" = "None",
           "2" = "Correl",
           "3" = "Ttest",
           "4" = "Ftest",
           "5" = "Zscore",
           "6" = "Chisq",
           "7" = "Beta",
           "8" = "Binom",
           "9" = "Gamma",
           "10" = "Poisson",
           "11" = "Normal",
           "12" = "Ftest_Nonc",
           "13" = "Chisq_Nonc",
           "14" = "Logistic",
           "15" = "Laplace",
           "16" = "Uniform",
           "17" = "Ttest_Nonc",
           "18" = "Weibull",
           "19" = "Chi",
           "20" = "Invgauss",
           "21" = "Extval",
           "22" = "Pval",
           "23" = "Logpval",
           "24" = "Log10pval",
           "1001" = "Estimate",      # estimate of some parameter
           "1002" = "Label",         # index into some set of labels
           "1003" = "Neuroname",     # index into the NeuroNames labels set
           "1004" = "Genmatrix",     # M x N matrix at each voxel
           "1005" = "Symmatrix",     # N x N symmetric matrix at each voxel
           "1006" = "Dispvect",      # a displacement field
           "1007" = "Vector",        # a displacement vector
           "1008" = "Pointset",      # a spatial coordinate
           "1009" = "Triangle",      # triple of indexes 
           "1010" = "Quaternion",    # a quaternion
           "1011" = "Dimless")       # Dimensionless value - no params
  }
  convert.form <- function(fc) {
    ## defgroup NIFTI1_XFORM_CODES
    ## brief nifti1 xform codes to describe the "standard" 
    ## coordinate system
    switch(as.character(fc),
           "0" = "Unkown",        # Arbitrary coordinates (Method 1)
           "1" = "Scanner_Anat",  # Scanner-based anatomical coordinates
           "2" = "Aligned_Anat",  # Coordinates aligned to another file's,
                                  # or to anatomical "truth"
           "3" = "Talairach",     # Coordinates aligned to Talairach-
                                  # Tournoux Atlas; (0,0,0)=AC, etc.
           "4" = "MNI_152")       # MNI 152 normalized coordinates
  }
  convert.units <- function(units) {
    ## defgroup NIFTI1_UNITS
    ## brief nifti1 units codes to describe the unit of measurement for
    ## each dimension of the dataset
    switch(as.character(units),
           "0" = "Unkown",        # unspecified units
           "1" = "meter",         # meters
           "2" = "mm",            # millimeters
           "3" = "micron",        # micrometers
           "8" = "sec",           # seconds
           "16" = "msec",         # milliseconds
           "24" = "usec",         # microseconds
           "32" = "Hz",           # Hertz
           "40" = "ppm",          # parts per million
           "48" = "rads")         # radians per second
  }
  convert.slice <- function(sc) {
    switch(as.character(sc),
           "0" = "Unknown",
           "1" = "Seq_Inc",       # sequential increasing
           "2" = "Seq_Dec",       # sequential decreasing
           "3" = "Alt_Inc",       # alternating increasing
           "4" = "Alt_Dec",       # alternating decreasing
           "5" = "Alt_Inc2",      # alternating increasing #2
           "6" = "Alt_Dec2")      # alternating decreasing #2
  }
  ## Bitwise conversion subroutines
  xyzt2space <- function (xyzt) {
    ## define XYZT_TO_SPACE(xyzt)       ( (xyzt) & 0x07 )
    require("bitops")
    bitAnd(xyzt, 7)
  }
  xyzt2time <- function (xyzt) {
    ## define XYZT_TO_TIME(xyzt)        ( (xyzt) & 0x38 )
    require("bitops")
    bitAnd(xyzt, 56)
  }
  dim2freq <- function(di) {
    ## define DIM_INFO_TO_FREQ_DIM(di)   ( ((di)     ) & 0x03 )
    require("bitops")
    bitAnd(di, 3)
  }
  dim2phase <- function(di) {
    ## define DIM_INFO_TO_PHASE_DIM(di)  ( ((di) >> 2) & 0x03 )
    require("bitops")
    bitAnd(bitShiftR(di, 2), 3)
  }
  dim2slice <- function(di) {
    ## define DIM_INFO_TO_SLICE_DIM(di)  ( ((di) >> 4) & 0x03 )
    require("bitops")
    bitAnd(bitShiftR(di, 4), 3)
  }

  ## Open appropriate file
  if(gzipped) {
    suffix <- ifelse(onefile, "nii.gz", "hdr.gz")
    fname <- paste(fname, suffix, sep=".")
    fid <- gzfile(fname, "rb")
  }
  else {
    suffix <- ifelse(onefile, "nii", "hdr")
    fname <- paste(fname, suffix, sep=".")
    fid <- file(fname, "rb")
  }

  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)
  ## Test for endian properties
  endian <- .Platform$endian
  sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
  if(sizeof.hdr != 348) {
    close(fid)
    endian <- "swap"
    if(gzipped)
      fid <- gzfile(fname, "rb")
    else
      fid <- file(fname, "rb")
    sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
    if(verbose) cat("  ENDIAN = swap", fill=TRUE)
  }

  ## was header_key substruct in Analyze 7.5
  data.type <- rawToChar(readBin(fid, "raw", n=10))
  db.name <- rawToChar(readBin(fid, "raw", n=18))
  extents <- readBin(fid, integer(), size=4, endian=endian)
  session.error <- readBin(fid, integer(), size=2, endian=endian)
  regular <- rawToChar(readBin(fid, "raw", n=1))
  dim.info <- readBin(fid, integer(), size=1, signed=FALSE, endian=endian)

  ## was image_dimension substruct in Analyze 7.5
  dim <- readBin(fid, integer(), 8, size=2, endian=endian)
  intent.p1 <- readBin(fid, numeric(), size=4, endian=endian)
  intent.p2 <- readBin(fid, numeric(), size=4, endian=endian)
  intent.p3 <- readBin(fid, numeric(), size=4, endian=endian)
  intent.code <- readBin(fid, integer(), size=2, endian=endian)
  datatype <- readBin(fid, integer(), size=2, endian=endian)
  bitpix <- readBin(fid, integer(), size=2, endian=endian)
  slice.start <- readBin(fid, integer(), size=2, endian=endian)
  pixdim <- readBin(fid, numeric(), 8, size=4, endian=endian)
  vox.offset <- readBin(fid, numeric(), size=4, endian=endian)
  scl.slope <- readBin(fid, numeric(), size=4, endian=endian)
  scl.inter <- readBin(fid, numeric(), size=4, endian=endian)
  slice.end <- readBin(fid, integer(), size=2, endian=endian)
  slice.code <- readBin(fid, integer(), size=1, signed=FALSE, endian=endian)
  xyzt.units <- readBin(fid, integer(), size=1, signed=FALSE, endian=endian)
  cal.max <- readBin(fid, numeric(), size=4, endian=endian)
  cal.min <- readBin(fid, numeric(), size=4, endian=endian)
  slice.duration <- readBin(fid, numeric(), size=4, endian=endian)
  toffset <- readBin(fid, numeric(), size=4, endian=endian)
  glmax <- readBin(fid, integer(), size=4, endian=endian)
  glmin <- readBin(fid, integer(), size=4, endian=endian)

  ## was data_history substruct in Analyze 7.5
  descrip <- rawToChar(readBin(fid, "raw", n=80))
  aux.file <- rawToChar(readBin(fid, "raw", n=24))

  qform.code <- readBin(fid, integer(), size=2, endian=endian)
  sform.code <- readBin(fid, integer(), size=2, endian=endian)

  quatern.b <- readBin(fid, numeric(), size=4, endian=endian)
  quatern.c <- readBin(fid, numeric(), size=4, endian=endian)
  quatern.d <- readBin(fid, numeric(), size=4, endian=endian)
  qoffset.x <- readBin(fid, numeric(), size=4, endian=endian)
  qoffset.y <- readBin(fid, numeric(), size=4, endian=endian)
  qoffset.z <- readBin(fid, numeric(), size=4, endian=endian)

  srow.x <- readBin(fid, numeric(), 4, size=4, endian=endian)
  srow.y <- readBin(fid, numeric(), 4, size=4, endian=endian)
  srow.z <- readBin(fid, numeric(), 4, size=4, endian=endian)

  intent.name <- rawToChar(readBin(fid, "raw", n=16))

  magic <- rawToChar(readBin(fid, "raw", n=4))
  if(!(magic %in% c("n+1","ni1")))
    stop(" -- Unrecognized \"magic\" field! --")

  ## Additional fields...
  freq.dim <- dim2freq(dim.info)
  phase.dim <- dim2phase(dim.info)
  slice.dim <- dim2slice(dim.info)
  intent <- convert.intent(intent.code)
  if(slice.code != 0 && slice.dim != 0 && slice.duration > 0)
    slice.name <- convert.slice(slice.code)
  else
    slice.name <- "Unknown"
  data.type <- ifelse(data.type == "", convert.datatype(datatype), "")
  vox.units <- convert.units(xyzt2space(xyzt.units))
  time.units <- convert.units(xyzt2time(xyzt.units))
  qform.name <- convert.form(qform.code)
  sform.name <- convert.form(sform.code)

  extender <- readBin(fid, integer(), 4, size=1, signed=FALSE, endian=endian)
  if(extender[1] != 0)
    stop("WARNING: Header extensions exist!")
  
  close(fid)

  nhdr <- list(sizeof.hdr, data.type, db.name, extents, session.error,
               regular, dim.info, dim, intent.p1, intent.p2, intent.p3,
               intent.code, datatype, bitpix, slice.start, pixdim,
               vox.offset, scl.slope, scl.inter, slice.end, slice.code,
               xyzt.units, cal.max, cal.min, slice.duration, toffset,
               glmax, glmin, descrip, aux.file, qform.code,
               sform.code, quatern.b, quatern.c, quatern.d, qoffset.x,
               qoffset.y, qoffset.z, srow.x, srow.y, srow.z, intent.name,
               magic, extender,
               endian, freq.dim, phase.dim, slice.dim, intent, slice.name,
               vox.units, time.units, qform.name, sform.name)

  names(nhdr) <- c("sizeof.hdr", "data.type", "db.name", "extents",
                   "session.error", "regular", "dim.info", "dim",
                   "intent.p1", "intent.p2", "intent.p3", "intent.code",
                   "datatype", "bitpix", "slice.start", "pixdim",
                   "vox.offset", "scl.slope", "scl.inter", "slice.end",
                   "slice.code", "xyzt.units", "cal.max", "cal.min",
                   "slice.duration", "toffset", "glmax", "glmin",
                   "descrip", "aux.file", "qform.code", "sform.code",
                   "quatern.b", "quatern.c", "quatern.d", "qoffset.x",
                   "qoffset.y", "qoffset.z", "srow.x", "srow.y",
                   "srow.z", "intent.name", "magic", "extender",
                   "endian", "freq.dim", "phase.dim", "slice.dim",
                   "intent", "slice.name", "vox.units", "time.units",
                   "qform.name", "sform.name")
  
  ## Warnings?
  options(warn=oldwarn)

  return(nhdr)
}

read.nifti.img <- function(fname, onefile=TRUE, gzipped=TRUE,
                           verbose=FALSE, warn=-1) {
  ## --- the original ANALYZE 7.5 type codes ---
  ## DT_NONE                    0
  ## DT_UNKNOWN                 0     # what it says, dude
  ## DT_BINARY                  1     # binary (1 bit/voxel)
  ## DT_UNSIGNED_CHAR           2     # unsigned char (8 bits/voxel)
  ## DT_SIGNED_SHORT            4     # signed short (16 bits/voxel)
  ## DT_SIGNED_INT              8     # signed int (32 bits/voxel)
  ## DT_FLOAT                  16     # float (32 bits/voxel)
  ## DT_COMPLEX                32     # complex (64 bits/voxel)
  ## DT_DOUBLE                 64     # double (64 bits/voxel)
  ## DT_RGB                   128     # RGB triple (24 bits/voxel)
  ## DT_ALL                   255     # not very useful (?)
  ## ----- another set of names for the same ---
  ## DT_UINT8                   2
  ## DT_INT16                   4
  ## DT_INT32                   8
  ## DT_FLOAT32                16
  ## DT_COMPLEX64              32
  ## DT_FLOAT64                64
  ## DT_RGB24                 128
  ## ------------------- new codes for NIFTI ---
  ## DT_INT8                  256     # signed char (8 bits)
  ## DT_UINT16                512     # unsigned short (16 bits)
  ## DT_UINT32                768     # unsigned int (32 bits)
  ## DT_INT64                1024     # long long (64 bits)
  ## DT_UINT64               1280     # unsigned long long (64 bits)
  ## DT_FLOAT128             1536     # long double (128 bits)
  ## DT_COMPLEX128           1792     # double pair (128 bits)
  ## DT_COMPLEX256           2048     # long double pair (256 bits)
  
  ## Read in data from the header file
  nhdr <- read.nifti.hdr(fname, onefile, gzipped, verbose, warn)

  if(nhdr$datatype == 2) {
    ## 1 byte unsigned integer
    what <- integer()
    size <- nhdr$bitpix / 8
    signed <- FALSE
  }
  else {
    if(nhdr$datatype == 4 | nhdr$datatype == 8) {
      ## 1 or 2 byte short integers
      what <- integer()
      size <- nhdr$bitpix / 8
      signed <- TRUE
    }
    else {
      if(nhdr$datatype == 16) {
        ## 4 byte floats
        what <- numeric()
        size <- nhdr$bitpix / 8
        signed <- TRUE
      }
      else {
        if(nhdr$datatype == 64) {
          ## 8 byte doubles
          what <- double()
          size <- nhdr$bitpix / 8
          signed <- TRUE
        }
        else {
          if(nhdr$datatype == 512) {
            ## 2 byte unsigned short integers
            what <- integer()
            size <- nhdr$bitpix / 8
            signed <- FALSE
          }
          else
            stop(paste("Data type ", nhdr$datatype, "unsupported in ", fname,
                       ".img", sep=""))
        }
      }
    }
  }

  n <- prod(nhdr$dim[2:5])

  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)

  ## Open appropriate file
  if(onefile) {
    if(gzipped)
      fid <- gzfile(paste(fname, "nii.gz", sep="."), "rb")
    else
      fid <- file(paste(fname, "nii", sep="."), "rb")
    skip <- readBin(fid, integer(), nhdr$vox.offset, size=1,
                    endian=nhdr$endian)
    img <-  readBin(fid, what, n, size, signed, endian=nhdr$endian)
  }
  else {
    if(gzipped)
      fid <- gzfile(paste(fname, "img.gz", sep="."), "rb")
    else
      fid <- file(paste(fname, "img", sep="."), "rb")
    img <- readBin(fid, what, n, size, signed, endian=nhdr$endian)
  }
  close(fid)

  ## Warnings?
  options(warn=oldwarn)

  ## return four-dimensional array (depends on nhdr$dim[1])
  array(img, nhdr$dim[1:nhdr$dim[1] + 1])
}
