#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <dirent.h>
#include <string>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <Rinternals.h>

void bcholesky(int* n, double* matrix,  int* lda, int* info) {
  F77_CALL(dpotrf)("L", n, matrix, lda, info);
}

void bloese(double* A, double* z, int* n, int* lda, int* incx) {
  F77_CALL(dtrsv)("L","T","N", n, A, lda, z, incx);
} 

void bloese2(double* A, double* z, int* n, int* lda, int* incx) {
  F77_CALL(dtrsv)("L","N","N", n, A, lda, z, incx);
} 

void baymem_gausssample(double* temp, int* noa)
{

for (int i=0; i< noa[0]; i++)
  {
    temp[i]=rnorm(0.0, 1.0);
  }

return;
}

void multivariate_normal(double* x, double* L, int* n, double* temp) {
  int i, c = 1, d = 0;
  //c = 1;
  //d = 0;
  bcholesky(n, L, n, &d);
  bloese2(L, x, n, n ,&c);
  bloese(L, x, n, n ,&c);
  baymem_gausssample(temp, n);
  bloese(L, temp, n, n , &c);
  for (i=0; i<*n; i++) {
    x[i] = x[i] + temp[i];
  }
}


double power(double x, int n)
{
  if (n==0){return 1;}
  if (n==1){return x;}
  return (x*power(x,n-1));
}


double baymem_aif(double t, double* aif_param, int aifmodel)
{
  double aif = 0.0;
  int a1=aif_param[0];
  int m1=aif_param[0];
  int a2=aif_param[0];
  int m2=aif_param[0];

  switch (aifmodel)
    {
    case 1:
      aif = (a1*exp(-m1*t)+a2*exp(-m2*t));
      break;
    case 2:
      aif = a1*t*exp(-m1*t)+a2*(exp(-m2*t)-exp(-m1*t));
    }
  return(aif);
}

double baymem_extraterm(double vp,double time, double* aif_param, int aifmodel)
{
  return(vp*baymem_aif(time, aif_param, aifmodel));
}

double baymem_convterm(double kep, double time, double* aif_param, int aifmodel)
{
  int a1=aif_param[0];
  int m1=aif_param[0];
  int a2=aif_param[0];
  int m2=aif_param[0];

  double bla=0.0;
  switch (aifmodel)
    {
    case 1:
      bla = (a1*(exp(-m1*(time))-exp(-kep*(time)))/(kep-m1)+a2*(exp(-m2*(time))-exp(-kep*(time)))/(kep-m2));
      break;
    case 2:
      bla = a1*kep*(time*exp(-m1*time)-((exp(-m1*time)-exp(-kep*time))/(kep-m1)))/(kep-m1);
      bla = bla + a2*kep*((exp(-m2*time)-exp(-kep*time))/(kep-m2)-(exp(-m1*time)-exp(-kep*time))/(kep-m1));
    }
  return bla;
}

double baymem_log_fc_gamma(double gamma, double tau_epsilon, double tau_gamma, double* conc, double* time, double kep, double vp, int voxel, int s, int aifmodel, double* aif_param, int nrcovar, double* beta0, double** x_covar, int T, int I)
 {
   double p = 0.0;
 
   for (int i=0;i<nrcovar;i++)
     {
       p += beta0[i]*x_covar[s][i];
     }
   p = -0.5*tau_gamma*pow(gamma-p,2.0);
   for (int t=0; t<T; t++)
     {
       p -= 0.5*tau_epsilon*power((conc[voxel+t*I])-baymem_extraterm(vp,time[t], aif_param, aifmodel)-exp(gamma)*baymem_convterm(kep,time[t], aif_param, aifmodel),2);
     }
   return p;
 }
double baymem_update_gamma2(double gamma, double kep, double vp, double tau_gamma, double tau_epsilon, double* conc, double* time, int voxel, double sigma, int s, int aifmodel, double* aif_param, int nrcovar, double* beta0, double** x_covar, int T, int I){

  double gamma_new = rnorm(gamma,sqrt(sigma));
    
  double logalpha = 0.0;

  logalpha +=  baymem_log_fc_gamma(gamma_new,tau_epsilon,tau_gamma,conc,time,kep,vp,voxel, s, aifmodel, aif_param, nrcovar, beta0, x_covar, T, I);
  logalpha -= baymem_log_fc_gamma(gamma,tau_epsilon,tau_gamma,conc,time,kep,vp,voxel, s, aifmodel, aif_param, nrcovar, beta0, x_covar, T, I);
  
  if (exp(logalpha)>runif(0.0,1.0))
    {
      return gamma_new;
    }
  else
    {
      return gamma;
    }
}

double baymem_log_fc_delta(double delta, double tau_epsilon, double tau_gamma, double* conc, double* time, double kep, double vp, int voxel, int s, int aifmodel, double* aif_param, int nrcovar, double* beta0, double** x_covar, int T, int I)
 {
   double p = 0.0;
   double help=0.0;
   for (int i=0;i<nrcovar;i++)
     {
       help += beta0[i]*x_covar[s][i];
     }
   p = -0.5*tau_gamma*delta*delta;
   for (int t=0; t<T; t++)
     {
       p -= 0.5*tau_epsilon*power((conc[voxel+t*I])-baymem_extraterm(vp,time[t], aif_param, aifmodel)-exp(delta+help)*baymem_convterm(kep,time[t], aif_param, aifmodel),2);
     }
   return p;
 }
double baymem_update_gamma2b(double gamma, double kep, double vp, double tau_gamma, double tau_epsilon, double* conc, double* time, int voxel, double sigma, int s, int aifmodel, double* aif_param, int nrcovar, double* beta0, double** x_covar, int T, int I){

  double gamma_new = rnorm(gamma,sqrt(sigma));
    
  double logalpha = 0.0;

  logalpha +=  baymem_log_fc_delta(gamma_new,tau_epsilon,tau_gamma,conc,time,kep,vp,voxel,s,aifmodel,aif_param,nrcovar,beta0,x_covar,T,I);
  logalpha -= baymem_log_fc_delta(gamma,tau_epsilon,tau_gamma,conc,time,kep,vp,voxel,s,aifmodel,aif_param,nrcovar,beta0,x_covar,T,I);
  
  if (exp(logalpha)>runif(0.0,1.0))
    {
      return gamma_new;
    }
  else
    {
      return gamma;
    }
}

double baymem_log_fc_theta2(double theta, double tau_epsilon, double tau_theta, double* conc, double* time, double ktrans, double vp, int voxel, int s, int aifmodel, double* aif_param, int nrcovar, double* beta_kep, double** x_covar, int T, int I)
 {
   double p = 0.0;
  for (int i=0;i<nrcovar;i++)
     {
       p += beta_kep[i]*x_covar[s][i];
     }
   p = -0.5*tau_theta*pow(theta-p,2.0);
   for (int t=1; t<=T; t++)
     {
       p -= 0.5*tau_epsilon*power((conc[voxel+t*I])-baymem_extraterm(vp,time[t],aif_param,aifmodel)-ktrans*baymem_convterm(exp(theta),time[t], aif_param, aifmodel),2);
     }
   return p;
 }
double baymem_update_theta2(double theta, double ktrans, double vp, double* conc, double* time, double tau_epsilon, double tau_theta, int voxel, double sigma, int s, int aifmodel, double* aif_param, int nrcovar, double* beta_kep, double** x_covar, int T, int I)
{
  double theta_new=rnorm(theta,sqrt(sigma));
  double logalpha = 0.0;

  logalpha += baymem_log_fc_theta2(theta_new,tau_epsilon,tau_theta,conc,time,ktrans,vp,voxel,s,aifmodel,aif_param,nrcovar,beta_kep,x_covar,T, I);
  logalpha -= baymem_log_fc_theta2(theta,tau_epsilon,tau_theta,conc,time,ktrans,vp,voxel,s,aifmodel,aif_param,nrcovar,beta_kep,x_covar,T,I);
    

  if (exp(logalpha)>runif(0.0,1.0))
    {
     return theta_new;

    }
  else
    {
     return theta;
      
    }
}

double baymem_log_fc_eta3(double eta, double tau_epsilon, double* conc, double* time, double kep, double ktrans, int voxel, int s, int aifmodel, double* aif_param, int nrcovar, double** x_covar, int T, int I, double a_vp, double b_vp)
 {
   double p = 0.0;
   double a,b;
   for (int t=1; t<=T; t++)
     {
       a = baymem_aif(time[t], aif_param, aifmodel);
       b=conc[voxel+t*I]-ktrans*baymem_convterm(kep,time[t], aif_param, aifmodel);
       p -= 0.5*tau_epsilon*(a*eta-b)*(a*eta-b);
     }
   p += (a_vp-1)*log(eta)+(b_vp-1)*log(1-eta);
   return p;
 }
double baymem_update_eta3(double vp, double kep, double ktrans, double tau_epsilon, double* conc, double* time, int voxel, double sigma, int s, int aifmodel, double* aif_param, int nrcovar, double** x_covar, int T, int I, double a_vp, double b_vp)
{

  double eta_new = -1;
  while (eta_new>1 || eta_new<0)
    {
      eta_new=rnorm(vp,sqrt(sigma));
    }
  double logalpha = 0.0;
  logalpha +=  baymem_log_fc_eta3(eta_new,tau_epsilon,conc,time,kep,ktrans,voxel,s,aifmodel,aif_param,nrcovar,x_covar,T,I,a_vp,b_vp);
  logalpha -= baymem_log_fc_eta3(vp,tau_epsilon,conc,time,kep,ktrans,voxel,s,aifmodel,aif_param,nrcovar,x_covar,T,I,a_vp,b_vp);
   if (exp(logalpha)>runif(0.0,1.0))
     {
       return eta_new;
     }
   else
     {
       return vp;
     }
}

void baymem_update_tau_epsilon1(double* tau, double a, double b, double** conc, double** vp, double** ktrans, double** kep, double** time, int* T, int* I, int nrscans, double* A1, double* M1, double* A2, double* M2, int aifmodel)
{

  double aa,bb,a1,m1,a2,m2;
  double* aif_param=new double[4];
  for (int s=0; s<nrscans; s++)
    {
      aif_param[0]=A1[s];
      aif_param[1]=M1[s];
      aif_param[2]=A2[s];
      aif_param[3]=M2[s];
      aa=a;
      bb=b;
      for (int voxel=0;voxel<I[s]; voxel++)
	{
	  for (int t=0;t<T[s]; t++)
	    {
	      aa += 0.5;
	      bb += .5*pow(conc[s][voxel+t*I[s]]-baymem_extraterm(vp[s][voxel], time[s][t],aif_param,aifmodel)-ktrans[s][voxel]*baymem_convterm(kep[s][voxel], time[s][t], aif_param, aifmodel),2.0);
	    }
	}
      tau[s]=rgamma(aa,bb);
    }
  return;
}

void baymem_update_beta(double* beta0, double** ktrans, double* tau, double *tau0, double** x_covar, double* M, double* V, double*** Vconst, double* temp1, double* temp2, double* temp3, int nrcovar, int nrscans, int fixedeff, int* I)
{
  int* n=new int[1];
    for (int k=0; k<nrcovar; k++)
    {
      M[k]=0.0;
      for (int s=0; s<nrscans; s++)
	{
	  if (x_covar[s][k]!=0)
	    {
	      for (int voxel=0; voxel<I[s]; voxel++)
		{
			      M[k]=M[k]+tau0[s]*x_covar[s][k]*log(ktrans[s][voxel]);
		}
	    }
	}
    }
 
    for (int i=0; i<(nrcovar*nrcovar); i++)
      {
	V[i]=0.0;
      }
    for (int i=0; i<nrcovar; i++)
      {
	for (int j=i; j<nrcovar; j++)
	  {
	    for (int s=0; s<nrscans; s++)
	      {
		V[i*nrcovar+j]+=(double)I[s]*tau0[s]*Vconst[s][i][j];
	      }
	  }
      }
    
    for (int i=fixedeff; i<nrcovar; i++)
      {
	V[i*nrcovar+i]=V[i*nrcovar+i]+tau[i];
      }
  
  int info=1;
  n[0]=nrcovar;
  multivariate_normal(M,V,n,temp1);

  for (int i=0; i<nrcovar; i++)
    {
      beta0[i]=M[i];
    }  

  return;
}

int baymem_tune(double* what,double acc, int i)
{
  int count=0;
  {
      if (acc<.3){what[i]=what[i]*.85;count=1;}
      if (acc<.2){what[i]=what[i]*.1;}
      if (acc<.1){what[i]=what[i]*.1;}
      if (acc<.02){what[i]=what[i]*.1;}
      if (acc>.6){what[i]=what[i]*1.1;count=1;}
      if (acc>.7){what[i]=what[i]*3;}
      if (acc>.8){what[i]=what[i]*5;}
      if (acc>.9){what[i]=what[i]*5;}
      if (acc>.99){what[i]=what[i]*15;}
  }
  return(count);
}
































extern "C" {
void baymem_main(int* nrscans0, int* conclist, double* maxvalue, double* timelist, 
int* I, int* T,  // vectors of Dimension of datafiles
double* A1, double* M1, double* A2, double* M2, // vectors with AIF parameters
int* nriters0, int* thinning0, int* burnin0, int* tuning0, double* tunepct0, int* retunecycles0, // MCMC settings
double* ab_beta, double* ab_gamma,
double* ab_theta, double* ab_vp,
double* ab_epsilon, //Priors
int* nrcovar0, int* fixedeff0, double* designmatrix,// number of covariates and number of fixedeffects, designmatrix of dim nrscans x nrcovar
int aifmodel, //1: Tofts-Kermode, 2: Orton??
		int* vpupdate0, int* verbose0,
		double* tau_beta, double* tau_beta_kep, double* tau_gamma, double* tau_theta, double* tau_epsilon, double* beta0, double* beta_kep,
double* ktrans0, double* kep0, double* vp0, double* sigmagamma0, double* sigmatheta0, double* sigmaeta0
)
{ 
  GetRNGstate();
  int nrscans=nrscans0[0];
  int nriters=nriters0[0];
  int thinning=thinning0[0];
  int burnin=burnin0[0];
  int tuning=tuning0[0];
  double tunepct=tunepct0[0];
  int retunecycles=retunecycles0[0];
  int nrcovar=nrcovar0[0];
  int fixedeff=fixedeff0[0];
  int vpupdate=vpupdate0[0];
  int verbose=verbose0[0];

  double a_beta=ab_beta[0];
  double a_gamma=ab_gamma[0];
  double a_theta=ab_theta[0];
  double a_vp=ab_vp[0];
  double a_epsilon=ab_epsilon[0];

  double b_beta=ab_beta[1];
  double b_gamma=ab_gamma[1];
  double b_theta=ab_theta[1];
  double b_vp=ab_vp[1];
  double b_epsilon=ab_epsilon[1];

 double** conc=new double*[nrscans];
 double** time=new double*[nrscans];
 int** to_tune=new int*[nrscans];

 double** sigmagamma=new double*[nrscans];
 int** acc_gamma=new int*[nrscans];
 double** sigmatheta=new double*[nrscans];
 int** acc_theta=new int*[nrscans];
 double** sigmaeta=new double*[nrscans];
 int** acc_eta=new int*[nrscans];

 int acccounter=0;

 double** ktrans=new double*[nrscans];
 double** kep=new double*[nrscans];
 double** vp=new double*[nrscans];

 int N1=0;

 int T0=0;
 int IT0=0;
 int I0=0;
 for (int s=0; s<nrscans; s++)
   {
     conc[s]=new double[I[s]*T[s]]; 
     to_tune[s]=new int[I[s]]; 
     time[s]=new double[T[s]]; 
     
     ktrans[s] = new double[I[s]];
     kep[s] = new double[I[s]];
     vp[s] = new double[I[s]];
     
     sigmagamma[s] = new double[I[s]];
     sigmatheta[s] = new double[I[s]];
     sigmaeta[s] = new double[I[s]];
     
     acc_gamma[s] = new int[I[s]];
     acc_theta[s] = new int[I[s]];
     acc_eta[s] = new int[I[s]];
     
     for (int t=0; t<T[s]; t++)
       {
	 time[s][t]=timelist[T0+t];
       }
     T0=T0+T[s];
     for (int i=0; i<(I[s]*T[s]); i++)
       {
	 conc[s][i]=maxvalue[s]*(double(conclist[IT0+i]))/(pow(2,15));
       }
     IT0=IT0+I[s]*T[s];
     for (int i=0; i<I[s]; i++)
       {
	 ktrans[s][i]=ktrans0[I0+i];
	 kep[s][i]=kep0[I0+i];
	 vp[s][i]=vp0[I0+i];

	 sigmagamma[s][i]=sigmagamma0[I0+i];
	 sigmatheta[s][i]=sigmatheta0[I0+i];
	 sigmaeta[s][i]=sigmaeta0[I0+i];
	 
	 acc_gamma[s][i]=0;
	 acc_theta[s][i]=0;
	 acc_eta[s][i]=0;
	 
	 to_tune[s][i]=1;
       }
     I0=I0+I[s];
     N1+=I[s];
   }  

//  if(verbose==1){Rprintf("Concentration prepared\n");}
 

 double* Mx1=new double[2*nrcovar];
 double* Mx2=new double[nrcovar];
 double* temp_vec=new double[1];
 double* M=new double[nrcovar];
 double* V=new double[nrcovar*nrcovar];
 double*** Vconst=new double**[nrscans];
 double* Mx3=new double[nrcovar];
 for (int i=0;i<nrscans;i++)
   {
   Vconst[i]=new double*[nrcovar];
   for (int j=0;j<nrcovar;j++)
     {
       Vconst[i][j]=new double[nrcovar];
     }
   }

 double** x_covar=new double*[nrscans];
 for (int i=0;i<(nrscans);i++)
   {
     x_covar[i]=new double[nrcovar];
     for (int j=0; j<nrcovar; j++)
	{
	  x_covar[i][j] = designmatrix[i+j*nrscans];
	}
   }


 for (int i=0; i<nrscans; i++)
   {
     for (int j=0;j<nrcovar;j++)
       {
	 for (int k=0;k<nrcovar;k++)
	   {
	     Vconst[i][j][k] = x_covar[i][j]*x_covar[i][k];
	   }
       }
   }

 double help;
 double* aif_param=new double[4]; 

















 // if(verbose==1){Rprintf("Starting iterations\n");}

 double temp,delta;
 int retune=0;
 for (int iter=1;iter<=nriters;iter++)
   {
     if (verbose==1)
       {
	 if(fmod(iter,100)==0)
	   {
	     Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b Iteration %i  ",iter);
	   }
       }


     acccounter++;
     for (int s=0; s<nrscans; s++)
       {
       aif_param[0]=A1[s];
       aif_param[1]=M1[s];
       aif_param[2]=A2[s];
       aif_param[3]=M2[s];
       help=0.0;
	 for (int i=0;i<nrcovar;i++)
	   {
	     help += beta0[i]*x_covar[s][i];
	   }
	 for (int voxel=0; voxel<I[s]; voxel++)
	   {
	    if (to_tune[s][voxel]==1)
		       {
			 // Update ktrans
			 temp=baymem_update_gamma2(log(ktrans[s][voxel]), kep[s][voxel], vp[s][voxel], tau_gamma[s], tau_epsilon[s], conc[s], time[s], voxel, sigmagamma[s][voxel], s, aifmodel, aif_param, nrcovar, beta0, x_covar, T[s], I[s]);

					    if (temp!=log(ktrans[s][voxel]))
			   {
			     acc_gamma[s][voxel]=acc_gamma[s][voxel]+1;
			     ktrans[s][voxel]=exp(temp);
			   }

			 // update kep
					    temp=baymem_update_theta2(log(kep[s][voxel]), ktrans[s][voxel], vp[s][voxel], conc[s], time[s], tau_epsilon[s], tau_theta[s], voxel, sigmatheta[s][voxel], s, aifmodel, aif_param, nrcovar, beta_kep, x_covar, T[s], I[s]);
			 if (temp!=log(kep[s][voxel]))
			   {
			     acc_theta[s][voxel]=acc_theta[s][voxel]+1;
			     kep[s][voxel]=exp(temp);
			   }


			 // update vp
			 if (vpupdate==1)
			   {
			     temp=baymem_update_eta3(vp[s][voxel], kep[s][voxel],  ktrans[s][voxel],  tau_epsilon[s], conc[s], time[s], voxel, s, 0, aifmodel, aif_param, nrcovar, x_covar, T[s], I[s], a_vp, b_vp);
			     if (temp!=vp[s][voxel])
			       {
				 acc_eta[s][voxel]=acc_eta[s][voxel]+1;
				 vp[s][voxel]=temp;
			       }  
			   }
	   }
	   } // end voxel loop 
       } // end scan loop

     double aa, bb;


     // sigma^2_s
     baymem_update_tau_epsilon1(tau_epsilon, a_epsilon, b_epsilon, conc, vp,ktrans,kep,time, T, I, nrscans, A1, M1, A2, M2, aifmodel);


     // beta
     
      baymem_update_beta(beta0, ktrans,tau_beta,tau_gamma, x_covar, M, V, Vconst, Mx1,Mx2,Mx3,nrcovar,nrscans,fixedeff, I);
      baymem_update_beta(beta_kep, kep,tau_beta_kep,tau_theta, x_covar, M, V, Vconst, Mx1,Mx2,Mx3,nrcovar,nrscans,fixedeff, I);


    // tau_0
     aa=a_gamma;
     bb=b_gamma;
     for (int s=0; s<nrscans; s++)
       {
	 double help=0.0;
	 for (int i=0;i<nrcovar;i++)
	   {
	     help += beta0[i]*x_covar[s][i];
	   }
   
 	 for (int voxel=0; voxel<I[s]; voxel++)
 	   {
 			 aa+=0.5;
 			 bb+=0.5*pow(log(ktrans[s][voxel])-help,2.0);     
 	   }
       
     tau_gamma[s] = rgamma(aa,bb);   
       }
    // tau_0
     aa=a_theta;
     bb=b_theta;
     for (int s=0; s<nrscans; s++)
       {
	 double help=0.0;
	 for (int i=0;i<nrcovar;i++)
	   {
	     help += beta_kep[i]*x_covar[s][i];
	   }
   
 	 for (int voxel=0; voxel<I[s]; voxel++)
 	   {
 			 aa+=0.5;
 			 bb+=0.5*pow(log(kep[s][voxel])-help,2.0);     
 	   }
       
     tau_theta[s] = rgamma(aa,bb);   
       }

     for (int k=fixedeff; k<nrcovar; k++)
       {
	 //tau
	 aa=a_beta+1/2.0;
	 bb=b_beta;
	 bb += 0.5*beta0[k]*beta0[k];
	 tau_beta[k] = rgamma(aa,1/bb);
	 
	 //tau
	 aa=a_beta+1/2.0;
	 bb=b_beta;
	 bb += 0.5*beta_kep[k]*beta_kep[k];
	 tau_beta_kep[k] = rgamma(aa,bb);
       }

  double* aif_param=new double[4];
       

    if (iter==tuning && retune!=retunecycles )
       {
	 if (verbose==1){Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bTuning..");}

	 int count1 = 0;
	 int count2 = 0;
	 int count3 = 0;
	 int count4 = 0;
	 int count5 = 0;
	 int count6 = 0;
	 int tu=0;
	 for (int s=0; s<nrscans; s++)
	   {
	     for (int i=0; i<I[s]; i++)
	       {
		 if (to_tune[s][i]==1)
		   {
		     count2 += acc_theta[s][i];
		     count3 += acc_gamma[s][i];
		     count6 += acc_eta[s][i];
		     tu=baymem_tune(sigmatheta[s],double(acc_theta[s][i])/double(acccounter),i);
		     tu+=baymem_tune(sigmagamma[s],double(acc_gamma[s][i])/double(acccounter),i);
 		     if (vpupdate==1){tu+=baymem_tune(sigmaeta[s],double(acc_eta[s][i])/double(acccounter),i);}
		     count4 += tuning;
		     if (tu>0)
		       {
			 count5++;
		       }
		     else
		       {
			 to_tune[s][i]=0;
		       }
		   }
	       }
	   }
	 
	 for (int s=0; s<nrscans; s++)
	   {
	     for (int i=0; i<I[s]; i++)
	       {
		 acc_theta[s][i]=0;
		 acc_eta[s][i]=0;
		 acc_gamma[s][i]=0;
		 acccounter=0;
	       }
	   }
	 iter=0;
	 if (verbose==1){Rprintf("%d to do.\n",count5);}

	 if (100*count5<N1*tunepct && retune<=(retunecycles-1))
	   {
	     retune++;
	     tunepct=tunepct/2.0;
	     for (int s=0; s<nrscans; s++)
	       {
		 for (int i=0; i<I[s]; i++)
		   {
		     to_tune[s][i]=1;
		   }
	       }
	   }
       }
   
   } // end iteration loop


 

 return;
}

    } // extern "C"
