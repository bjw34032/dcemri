\name{dce.bayes}
\alias{dce.bayes}
\title{Bayesian Methods for Pharmacokinetic Modeling of Dynamic
  Contrast-Enhanced MRI Data}
\description{
  Bayesian analysis of contrast agent concentration time curves from
  DCE-MRI.
}
\usage{
dce.bayes(conc, time, ...)
}
\arguments{
  \item{conc}{Matrix or array of concentration time series.
    Last dimension must be time.}
  \item{time}{Time in minutes.}
  \item{mask}{Mask matrix or array. Voxels with mask=0 will be excluded.}
  \item{nriters}{Total number of iterations.}
  \item{thin}{Thining factor.}
  \item{burnin}{Number of iterations for burn-in.}
  \item{tune}{Number for iterations for tuning. Algorithm will be 
    tuned to an acceptance rate of 0.3 to 0.6.}
  \item{tau.ktrans}{Variance parameter for prior on $\log(K^{trans})$.}
  \item{tau.kep}{Variance parameter for prior on $\log(k_{ep})$.}
  \item{ab.vp}{Hyper prior parameters for $v_p$ Beta prior.}
  \item{ab.tauepsilon}{Hyper prior parameters for observation error Gamma Prior.}
  \item{aif.model}{Model for arterial input funktion. So far implemented:
    0 or "tofts-kermode" double exponential model by Tofts,Kermode;
    1 or "orton" exponential model by Orton et al.}
  \item{aif.parameter}{Vector of AIF parameters. 
    For Tofts/Kermode: $a_1$,$m_1$,$a_2$,$m_2$; for Orton et al.: $A_b$,$\mu_b$,$A_g$,$\mu_g$.}
  \item{vp.do}{If TRUE, model includes vascular term.}
  \item{samples}{If TRUE, output includes samples.}
}
\value{
}
\details{
}
\references{
}
\seealso{
}
\examples{
}
\author{Volker Schmid}
\keyword{}
