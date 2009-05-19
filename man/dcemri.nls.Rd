\name{dcemri.nls}
\alias{dcemri.nls}
\title{Routine to Fit Pharmacokinetic Models to DCE-MRI Data}
\description{
  These functions provide kinetic parameter estimates for the solutions
  to compartmental models used in dynamic perfusion MRI.  Non-linear
  least squares is performed using the Levenburg-Marquardt algorithm.
}
\usage{
dcemri.nls(conc, time, img.mask, model="extended", aif="tofts.kermode", nprint=0)
}
\arguments{
  \item{conc}{Contrast agent concentration (4D)}
  \item{time}{Acquisition times in minutes (pre-dose must be negative)}
  \item{img.mask}{Binary mask the same size as the contrast agent
    concentration data}
  \item{model}{Type of compartmental model.  Valid models currently
    include "weinmann", "extended", "orton.exp" and "orton.cos"}
  \item{aif}{Arterial input function.  Only "tofts.kermode" and
    "fritz.hansen" are allowed in combination with the "weinmann" and
    "extended" compartmental models, while "orton.exp" and "orton.cos"
    compartmental models require the "user" option if the default
    "orton.exp" and "orton.cos" AIF parameters are not used.}
  \item{nprint}{An integer that enables controlled printing of iterates
    if it is positive (for Levenburg-Marquardt algorithm)} 
}
\value{
  \item{ktrans}{Estimated $K^{trans}$}
  \item{kep}{Estimated $k_{ep}$}
  \item{ktranserror}{Error associated with $K^{trans}$}
  \item{keperror}{Error associated with $k_{ep}$}
  \item{ve}{Estimated volume of extraction}
  \item{vp}{Estimated volume in the plasma space}
  \item{time}{Acquisition times (in minutes)}
  \item{sse}{Sum of squared errors using the compartmental model}
}
\details{
  See references.
}
\references{
  Parker GJM, Padhani AR.
  {$T_1$-w} {DCE-MRI}: {$T_1$-weighted} dynamic contrast-enhanced {MRI}.
  In P~Tofts, ed., {\em Quantitative MRI of the Brain: Measuring Changes
    Caused by Disease}.
  Chichester, UK: Wiley, 2003; 341--364.

  Weinmann HJ, Laniado M, Mutzel W.
  Pharmacokinetics of {Gd-DTPA}/dimeglumine after intraveneous injection
  into healthy volunteers.
  {\em Physiological Chemistry and Physics and Medical NMR}
  1984;16:167--172.
}
\seealso{
}
%\examples{}
\author{Volker Schmid, Brandon Whitcher}
\keyword{}
