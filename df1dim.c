#define NRANSI
#include "nrutil.h"

extern int ncom;
extern double *pcom,*xicom,(*nrfunc)(double []);
extern void (*nrdfun)(double [], double []);

double df1dim(double x)
{
	int j;
	double df1=0.0;
	double *xt,*df;

	xt=dvector(1,ncom);
	df=dvector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	(*nrdfun)(xt,df);
	for (j=1;j<=ncom;j++) df1 += df[j]*xicom[j];
	free_dvector(df,1,ncom);
	free_dvector(xt,1,ncom);
	return df1;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software ,5#!0JX). */
