#include <math.h>
#define NRANSI
#include "nrutil.h"
/*#define ITMAX 200*/
#define EPS 1.0e-10
#define FREEALL free_dvector(xi,1,n);free_dvector(h,1,n);free_dvector(g,1,n);

void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	    double (*func)(double []), void (*dfunc)(double [], double []),
	    int ITMAX)
{
  extern void dlinmin(double p[], double xi[], int n, double *fret, 
	       double (*func)(double []),
	       void (*dfunc)(double [], double []));
  int j,its;
  double gg,gam,fp,dgg;
  double *g,*h,*xi;
  
  g=dvector(1,n);
  h=dvector(1,n);
  xi=dvector(1,n);
  fp=(*func)(p);
  (*dfunc)(p,xi);
 
  for (j=1;j<=n;j++) 
    {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j];
    }
  for (its=1;its<=ITMAX;its++) 
    {
      *iter=its;
      dlinmin(p,xi,n,fret,func,dfunc);
   
      if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
	FREEALL
	  return;
      }
      fp=(*func)(p);
     
      (*dfunc)(p,xi);
      dgg=gg=0.0;
      for (j=1;j<=n;j++) {
	gg += g[j]*g[j];
	dgg += (xi[j]+g[j])*xi[j];
      }
      if (gg == 0.0) {
	FREEALL
	  return;
      }
      gam=dgg/gg;
      for (j=1;j<=n;j++) {
	g[j] = -xi[j];
	xi[j]=h[j]=g[j]+gam*h[j];
      }
    }
  
  FREEALL
    return;
  
  /*nrerror("Too many iterations in frprmn");*/
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software ,5#!0JX). */
