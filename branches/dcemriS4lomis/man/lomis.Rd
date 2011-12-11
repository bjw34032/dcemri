\name{lomis}
\title{Longitudinal Medical Imaging Studies}
\description{
does, what it does
}
\usage{
lomis(scans, masks, timelines,
scandata, 
fixed=c(1,"post","post*treatment"), fixed.linear=NULL, random=c("patient","patient*post*treatment"), random.linear=NULL,
model="weinmann", aif="tofts-kermode",
prior.control=list("beta"=c(1,1),"gamma"=c(1,0.00001),"theta"=c(1,0.00001),"epsilon"=c(1,0.01),"vp"=c(1,19)),
mcmc.control=list("nriters"=10000,"burnin"=20000,"thinning"=100,"retunecycles"=3,"tunepct"=3),
multicore=(require(multicore)), verbose=FALSE)
}
\arguments{
  \item{scans}{is a list of multidimensional (1D-4D) arrays of contrast agent
    concentrations.  The last dimension is assumed to be temporal, while
    the previous dimensions are assumed to be spatial.}
  \item{masks}{is a list of (logical) multidimensional arrays that identifies
    the voxels to be analyzed. Each array has to have same dimension as \code{scans}
    minus temporal dimension.}
  \item{timelines}{is a list of vectors of acquisition times (in minutes) relative to 
    injection of the contrast agent.}
  \item{scandata}{is a data frame of covariates.}
  \item{fixed}{is a vector of variables from \code{scandata}, which are
    used as fixed categorical effects. Can contain \code{-1} to remove
    intercept. Usually contains post treatment scan or time variable,
    interaction of post treatment scan and group (e.g. treatment), and
    categorical covariates. Interactions are added e.g. by
    \code{"post*treatment"}.}
  \item{fixed.linear}{is a vector of variables from \code{scandata}, which are
    used as fixed linear effects (e.g. continous covariates).}
  \item{random}{is a vector of variables from \code{scandata}, which are
    used as fixed categorical effects. Usually contains patient and
    patient/treatment interaction.}
  \item{random.linear}{is a vector of variables from \code{scandata}, which are
    used as random linear effects.}
  \item{model}{is a character string that identifies the type of
    compartmental model to be used.  Acceptable models include:
    \describe{
      \item{\dQuote{weinmann}}{Tofts & Kermode AIF convolved with single
	compartment model}
      \item{\dQuote{extended}}{Weinmann model extended with additional
	vascular compartment (default)}
      \item{\dQuote{orton.exp}}{Extended model using Orton's exponential
	AIF} 
      \item{\dQuote{kety.orton.exp}}{Kety model using Orton's
	exponential AIF}
    }
  }
  \item{aif}{is a character string that identifies the parameters of the
    type of arterial input function (AIF) used with the above model.
    Acceptable values are:
    \itemize{
      \item{\code{tofts.kermode}}{(default) for the \code{weinmann} and
	\code{extended}	models}
      \item{\code{fritz.hansen}}{for the \code{weinmann} and
	\code{extended} models}
      \item{\code{orton.exp}}{(default) for the \code{orton.exp} and
	\code{kety.orton.exp} model} 
      \item{vector of parameters a1, m1, a2, m2}
      \item{data frame with patient specific paramters (a1, m1, a2, m2)}
    }
  }
  \item{control}{is a list of parameters used by \code{nls.lm.control}
    that are set by default, but may be customized by the user.}
  \item{multicore}{is a logical variable (default = \code{(require(multicore))}) that
    allows parallel processing via \pkg{multicore}.}
  \item{verbose}{is a logical variable (default = \code{FALSE}) that
    allows text-based feedback during execution of the function.}
  \item{...}{Additional parameters to the function.}
}
\details{
See Schmid et al. MRM 2009
}
\value{
So far: nothing...
}
\references{
  Schmid et al. (2009) MRM
}
\seealso{
  \code{\link{dcemri.bayes}}, \code{\link{dcemri.map}},
  \code{\link{dcemri.spline}}
}
\examples{
#dont run
if(0)
{


setwd("~/bioimg/data/DCEMRI/breast")
source("readbreastdata.R")
scans<-masks<-timelines<-list()
k<-0
for (i in 1:6)
for (j in c("pre","post"))
{
k<-k+1
breast<-breastdata(i,j,"./")
scans[[k]]<-breast$conc
masks[[k]]<-breast$mask
timelines[[k]]<-breast$time
}
scandata<-data.frame("post"=rep(c(0,1),6),"patient"=rep(1:6,each=2))
library(dcemriS4)
lom<-lomis(scans,masks,timelines,scandata,fixed=c("post"),random=c("patient","patient*post"),verbose=TRUE,aif="fritz.hansen",model="weinmann",multicore=TRUE)

}
}
\author{Volker Schmid \email{volkerschmid@users.sourceforge.net}}
\keyword{models}
