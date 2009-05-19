#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

double nulleins();
double RNDGAM(double a, double b);
double normal(double m, double s);
void gausssample(double* temp, int* noa);
double reins();

