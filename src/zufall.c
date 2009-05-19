#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

double nulleins()
{
  return runif(0,1);
}
// -------------------------------------
double RNDGAM(double a, double b){
// -------------------------------------
int accept=0;
double c,d,u,v,w,x,y,z;
if(a>1){                        /*   Algorithmus S.410 Devroye */
  c=a-1;
  d=3*a-3/4;
  while(accept==0){
    u=nulleins();
    v=nulleins();
    w=u*(1-u);
    y=sqrt(d/w)*(u-0.5);
    x=c+y;
    if(x >= 0){
      z=64*w*w*w*v*v;
      if( z<= (1-(2*y*y/x))) accept=1; 
      if (accept==0)
        if (log(z)<=2*((c*log(x/c))-y)) accept=1; 
      }
  }
}
else{                               /* Fall: a<=1; Stuart's theorem */
  x = pow(nulleins(),(1/a))*RNDGAM(a+1,1);
}
return x/b;


}
double reins()
{
	return nulleins() * 2.0 - 1.0;
}



double normal(double m, double s)
{
  return (rnorm(0,1)*sqrt(s)+m);
}


//Erzeugt Normalverteilten Zufallsvektor der Laenge noa
void gausssample(double* temp, int* noa)
{
  int i;
for (i=0; i< *noa; i++)
  {
    temp[i]=rnorm(0,1);
  }
return;
}
