.do.dcemri.lm<-function(list,model,aif,multicore)
{
dclm<-(dcemri.lm(list$conc,list$time,list$mask,model=model,aif=aif,multicore=multicore))
return(dclm)
}	

.lomis.sample<-function(i,sample,nriters,nrscans,cons,maxvalue,times,I,T,A1,M1,A2,M2,mcmc.control,prior.control,nrcovar,fixedeffects,xcovar,model.b,vpupdate)
{
results<-vector("list",10)

for (i in 1:nriters)
{
sample<-.C("baymem_main",	
   as.integer(nrscans),
   as.integer(concs),
   as.double(maxvalue),
   as.double(times),
   as.integer(I),
   as.integer(T),
   as.double(A1),
   as.double(M1),
   as.double(A2),
   as.double(M2), #10
 
  as.integer(mcmc.control$thin),
  as.integer(thin),
  as.integer(0),
  as.integer(0),
  as.double(0),
  as.integer(0),
  as.double(prior.control$beta),
  as.double(prior.control$gamma),
  as.double(prior.control$theta),
  as.double(prior.control$vp), #20

  as.double(prior.control$epsilon),
  as.integer(nrcovar),
  as.integer(fixedeffects),
  as.double(as.matrix(xcovar)),
  as.integer(model.b),
  as.integer(vpupdate),
  as.integer(0),
  as.double(sample[[28]]),
  as.double(sample[[29]]),
  as.double(sample[[30]]),

  as.double(sample[[31]]),
  as.double(sample[[32]]),
  as.double(sample[[33]]),
  as.double(sample[[34]]),
  as.double(sample[[35]]),
  as.double(sample[[36]]),
  as.double(sample[[37]]),
  as.double(sample[[38]]),
  as.double(sample[[39]]),
  as.double(sample[[40]]),
                    PACKAGE="dcemriS4"   
)

for (i in 1:10)
results[[i]]<-c(results[[i]],as.vector(sample[[27+i]]))
}

return(results)
}

lomis<-function(scans, masks, timelines,
scandata, 
fixed=c(1,"post","post*treatment"), fixed.linear=NULL, random=c("patient","patient*post*treatment"), random.linear=NULL,
model="weinmann", aif="tofts.kermode",
prior.control=list("beta"=c(1,1),"gamma"=c(1,0.00001),"theta"=c(1,0.00001),"epsilon"=c(1,0.01),"vp"=c(1,19)),
mcmc.control=list("nriters"=100000,"burnin"=20000,"thinning"=100,"retunecycles"=3,"tunepct"=10, "tuning"=413, "chains"=20),
multicore=(require(multicore)), verbose=FALSE)
{
if (multicore){require(multicore)}

nrscans <- length(scans)
T<-rep(NA,nrscans)
I<-rep(NA,nrscans)
concs<-c()
maxvalue<-rep(NA,nrscans)
times<-c()
firsttunelist<-list()
for (i in 1:nrscans)
{
I[i]<-sum(masks[[i]])
dims<-dim(scans[[i]])
T[i]<-dims[length(dims)]
conc<-as.vector(scans[[i]])
mask<-rep(as.vector(masks[[i]]),T[i])
conc<-conc[mask==1]
maxvalue[i]<-max(conc)
conc<-as.integer(2^15*conc/maxvalue[i])
concs<-c(concs,conc)
times<-c(times,timelines[[i]])
firsttunelist[[i]]<-list("conc"=scans[[i]],"mask"=masks[[i]],"time"=timelines[[i]])
}

xcovar<-NULL
if(fixed[1]!=(-1))
{
xcovar<-data.frame("Intercept"=rep(1,nrscans))
if(fixed[1]==1)fixed<-fixed[-1]
}

for (i in 1:length(fixed))
{
w<-which(names(scandata)==fixed[i])
if (length(w)==1)
{
x<-sort(unique(scandata[,w]))
for (xx in x[-1])
{
xcovar<-cbind(xcovar,ifelse(scandata[,w]==xx,1,0))
colnames(xcovar)[dim(xcovar)[2]]<-paste(fixed[i],xx,sep="")
}
}
if(length(w)==0)
{
w<-NULL
for (j in 1:length(names(scandata)))
if(grepl(names(scandata)[j],fixed[i]))w<-c(w,j)
x0<-apply(scandata[,w],1,prod)
x<-sort(unique(x0))
for (xx in x[-1])
{
xcovar<-cbind(xcovar,ifelse(x0==xx,1,0))
colnames(xcovar)[dim(xcovar)[2]]<-paste(fixed[i],xx,sep="")
}
}
}

if(!is.null(fixed.linear))
for (i in 1:length(fixed.linear))
{
w<-which(names(scandata)==fixed.linear[i])
if (length(w)==1)
xcovar<-cbind(xcovar,scandata[,w])
if(length(w)==0)
{
w<-NULL
for (j in 1:length(names(scandata)))
if(grepl(names(scandata)[j],fixed.linear[i]))w<-c(w,j)
}
xcovar<-cbind(xcovar,apply(scandata[,w],1,prod))
}

fixedeffects<-dim(xcovar)[2]

for (i in 1:length(random))
{
w<-which(names(scandata)==random[i])
if (length(w)==1)
{
x<-sort(unique(scandata[,w]))
for (xx in x)
{
xcovar<-cbind(xcovar,ifelse(scandata[,w]==xx,1,0))
colnames(xcovar)[dim(xcovar)[2]]<-paste(random[i],xx,sep="")
}
}
if(length(w)==0)
{
w<-NULL
for (j in 1:length(names(scandata)))
if(grepl(names(scandata)[j],random[i]))w<-c(w,j)
x0<-apply(scandata[,w],1,prod)
x<-sort(unique(x0))
for (xx in x[-1])
{
xcovar<-cbind(xcovar,ifelse(x0==xx,1,0))
colnames(xcovar)[dim(xcovar)[2]]<-paste(random[i],xx,sep="")
}
}
}

if(!is.null(random.linear))
for (i in 1:length(random.linear))
{
w<-which(names(scandata)==random.linear[i])
if (length(w)==1)
xcovar<-cbind(xcovar,scandata[,w])
if(length(w)==0)
{
w<-NULL
for (j in 1:length(names(scandata)))
if(grepl(names(scandata)[j],random.linear[i]))w<-c(w,j)
}
xcovar<-cbind(xcovar,apply(scandata[,w],1,prod))
}

switch(model)
{
weinmann={
vpupdate<-0
model.b<-1}
extended={
vpupdate<-1
model.b<-1}
orton.exp={
vpupdate<-1
model.b<-2}
kety.orton.exp={
vpupdate<-0
model.b<-2}
}

if(length(aif)==1)
{
aifp<-aifParameters(aif)
A1=rep(aifp$D*aifp$a1,nrscans)
M1=rep(aifp$m1,nrscans)
A2=rep(aifp$D*aifp$a2,nrscans)
M2=rep(aifp$m2,nrscans)
}
else
  {
    if(length(aif)==4)
      {
        A1<-rep(aif[1],nrscans)
        M1<-rep(aif[2],nrscans)
        A2<-rep(aif[3],nrscans)
        M2<-rep(aif[4],nrscans)
      }
    else
      {
        if(length(aif)==5)
          {
            A1<-rep(aif[1]*aif[2],nrscans)
            M1<-rep(aif[3],nrscans)
            A2<-rep(aif[1]*aif[4],nrscans)
            M2<-rep(aif[5],nrscans)
          }
        else
          {
            if(dim(aif)[2]==4)
              {
                A1<-aif[,1]
                M1<-aif[,2]
                A2<-aif[,3]
                M2<-aif[,4]
              }
            else
              {
                if(dim(aif)[2]==5)
                  {
                    A1<-aif[,1]*aif[,2]
                    M1<-aif[,3]
                    A2<-aif[,1]*aif[,4]
                    M2<-aif[,5]
                  }
              }
          }
      }
  }

if (is.null(mcmc.control$nriters))mcmc.control$nriters<-100000
if (is.null(mcmc.control$burnin))mcmc.control$burnin<-20000
if (is.null(mcmc.control$thinning))mcmc.control$thinning<-100
if (is.null(mcmc.control$tunepct))mcmc.control$tunepct<-10
if (is.null(mcmc.control$retunecycles))mcmc.control$retunecycles<-3
if (is.null(mcmc.control$tuning))mcmc.control$tuning<-3

nrcovar=dim(xcovar)[2]


cat("Data Prepared. (");
cat(proc.time())
cat(")\nEstimating initial values..")
##First Step: Estimate Ktrans values per voxel
##ToDo: changes for user-specific AIF

require(minpack.lm)
if (multicore)
{
firsttunelist<-mclapply(firsttunelist,.do.dcemri.lm,model,aif,multicore=FALSE,mc.preschedule=FALSE)
}
else
{
firsttunelist<-lapply(firsttunelist,.do.dcemri.lm,model,aif,multicore=FALSE)
}

print(firsttunelist[[1]])

cat(".done (")
cat(proc.time())
cat(").\nInitital tuning..")
## Step 2: Tuning
# use estimated kinetic values
ktrans<-kep<-vp<-c()
for (i in 1:nrscans)
{
kt<-as.vector(firsttunelist[[i]]$ktrans)
kp<-as.vector(firsttunelist[[i]]$kep)
if (vpupdate==1)vp0<-as.vector(firsttunelist[[i]]$vp)
mask<-as.vector(masks[[i]])
ktrans<-c(ktrans,kt[mask==1])
kep<-c(kep,kp[mask==1])
if (vpupdate==1)vp<-c(vp,v[mask==1])
}
if (vpupdate==0)vp<-rep(0,sum(I))

ktrans[is.na(ktrans)]<-1
kep[is.na(kep)]<-1
vp[is.na(vp)]<-0.05

initial<-.C("baymem_main",
                                        #int nrscans, int* conclist, double* maxvalue, double* timelist, 

   as.integer(nrscans),
   as.integer(concs),
   as.double(maxvalue),
   as.double(times),
                                        #int* I, int* T,  // vectors of Dimension of datafiles
   as.integer(I),
   as.integer(T),
                                        #double* A1, double* M1, double* A2, double* M2, // vectors with AIF parameters
  as.double(A1),
   as.double(M1),
   as.double(A2),
   as.double(M2), #10
                                        #int nriters, int thinning, int burnin, int tuning, double tunepct, int retunecyles, // MCMC settings
as.integer(mcmc.control$tuning+2),
as.integer(1),
as.integer(mcmc.control$tuning+1),
as.integer(mcmc.control$tuning),
as.double(mcmc.control$tunepct),
as.integer(mcmc.control$retunecycles),
#double* ab_beta, double* ab_gamma,
as.double(prior.control$beta),
as.double(prior.control$gamma),
#double* ab_theta, double ab_vp,
as.double(prior.control$theta),
as.double(prior.control$vp), #20
#double ab_epsilon, //Priors
as.double(prior.control$epsilon),
#int nrcovar, int fixedeff, double* designmatrix,// number of covariates and number of fixedeffects, designmatrix of dim nrscans x nrcovar
as.integer(nrcovar),
as.integer(fixedeffects),
as.double(as.matrix(xcovar)),
#int aifmodel, //1: Tofts-Kermode, 2: Orton.exp
#int vpupdate, int verbose
as.integer(model.b),
as.integer(vpupdate),
   as.integer(verbose),
#double* tau_beta, double* tau_beta_kep, double* tau_gamma, double* tau_theta, double* tau_epsilon, double* beta0, double* beta_kep,
as.double(rep(prior.control$beta[1]/prior.control$beta[2],dim(xcovar)[2])),
as.double(rep(prior.control$beta[1]/prior.control$beta[2],dim(xcovar)[2])),
as.double(rep(prior.control$gamma[1]/prior.control$gamma[2],nrscans)), #30
as.double(rep(prior.control$theta[1]/prior.control$theta[2],nrscans)),
as.double(rep(prior.control$epsilon[1]/prior.control$epsilon[2],nrscans)),
as.double(c(1,rep(0,dim(xcovar)[2]-1))),
as.double(c(1,rep(0,dim(xcovar)[2]-1))),
#double* ktrans0, double* kep0, double* vp0, double* sigmagamma0, double* sigmatheta0, double* sigmaeta0
as.double(ktrans),
as.double(kep),
as.double(vp),
as.double(rep(10,sum(I))),
as.double(rep(1,sum(I))),
as.double(rep(1,sum(I))),
                    PACKAGE="dcemriS4"   
)


cat(".done (")
cat(proc.time())
cat(").\nBurnin..")


sample<-.C("baymem_main",	
   as.integer(nrscans),
   as.integer(concs),
   as.double(maxvalue),
   as.double(times),
   as.integer(I),
   as.integer(T),
   as.double(A1),
   as.double(M1),
   as.double(A2),
   as.double(M2), #10
 
  as.integer(mcmc.control$burnin),
  as.integer(1),
  as.integer(mcmc.control$burnin-1),
  as.integer(0),
  as.double(0),
  as.integer(0),
  as.double(prior.control$beta),
  as.double(prior.control$gamma),
  as.double(prior.control$theta),
  as.double(prior.control$vp), #20

  as.double(prior.control$epsilon),
  as.integer(nrcovar),
  as.integer(fixedeffects),
  as.double(as.matrix(xcovar)),
  as.integer(model.b),
  as.integer(vpupdate),
  as.integer(0),
  as.double(initial[[28]]),
  as.double(initial[[29]]),
  as.double(initial[[30]]),

  as.double(initial[[31]]),
  as.double(initial[[32]]),
  as.double(initial[[33]]),
  as.double(initial[[34]]),
  as.double(initial[[35]]),
  as.double(initial[[36]]),
  as.double(initial[[37]]),
  as.double(initial[[38]]),
  as.double(initial[[39]]),
  as.double(initial[[40]]),
                    PACKAGE="dcemriS4"   
)

cat(".done (")
cat(proc.time())
cat(").\nSampling..")

cores<-1
if (multicore)
{
	if (is.null(mcmc.control$chain))
	{
	cores=multicore:::detectCores()
	}
	else
	{
	cores<-mcmc.control$chains
	}
}
nriters<-round(mcmc.control$nriters/cores)

if(multicore)
{
todo<-1:cores
samples<-mclapply(todo,.lomis.sample,sample,nriters,nrscans,cons,maxvalue,times,I,T,A1,M1,A2,M2,mcmc.control,prior.control,nrcovar,fixedeffects,xcovar,model.b,vpupdate)
results<-vector("list",10)
for (i in todo)
for (j in 1:10)
results[[j]]<-c(results[[j]],samples[[i]][[j]])
}
else
{
results<-.lomis.sample(0,sample,nriters,nrscans,cons,maxvalue,times,I,T,A1,M1,A2,M2,mcmc.control,prior.control,nrcovar,fixedeffects,xcovar,model.b,vpupdate)
}

cat("..done(")
cat(proc.time())
cat(").\nSorry, result reconstruction not implemented yet.")
return(results)
#end lomis()
}

