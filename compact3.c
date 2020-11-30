/* #################################################################### */
/* This function preforms Congugate Gradient steps to compact 
   distant MTs (in the yz plane) */ 

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "global_var.h"
#include "ran1.h"
#define EPS 1.0e-10 
#define  crep 7  /* LJ powers: repulsion */
#define  catt 3  /* LJ powers: attraction */
#define epsilon 0.01 /* LJ energy coef */
#define kbond 1000000.0 /* Spring Const */

extern void dlinmin(double *p, double *xi, int n, double *fret, 
		    double (*potx)(double *xcoor));
extern void frprmn(double *p, int n, double ftol, int *iter, double *fret,
		   double (*potx)(double *xcoor), void (*dpotx)(double *xcoor, 
							     double *grad),
		   int MXIT);

double potx(double *xcoor);
void dpotx(double *xcoor, double *grad);



void compactx()
{
  /* Merge MTs a long x-direction. Clusteres of bounded MTs move together,
     since they are bounded by strong springs. Cohesive interaction between
     MTs is obtained by a LJ potential.
  */  
  
  int iMT,jMT,nMT,nriter;
  double emin,ftol;
  double *xcoor,xi,xj;
  int MXIT;
    
  nMT=MT.number;
  ftol=1.0e-12; /* parameter for minimization convergaence */ 
  MXIT=500;     /* Max number of minimization steps */
  
  /* Initiate xcoor vector */
  xcoor=dvector(1,nMT);
  for (iMT=1;iMT<=nMT;iMT++)
    xcoor[iMT]=(double)MT.cm[iMT].x;
	   
  /* Apply NR Congugate Gradient routine Ch. 10.6; P. 423*/
  frprmn(xcoor,nMT,ftol,&nriter,&emin,potx,dpotx,MXIT);
 
  /* Store the new coordinates */
  if (nriter<MXIT)
    for (iMT=1;iMT<=nMT;iMT++)
      {
	MT.pend[iMT].x += (float)xcoor[iMT]-MT.cm[iMT].x;
	MT.mend[iMT].x += (float)xcoor[iMT]-MT.cm[iMT].x;
	MT.cm[iMT].x=(float)xcoor[iMT];
	
	for(jMT=1;jMT<=nMT;jMT++)
	  if (!overlap(iMT,jMT)>0)
	    {
	      ovlp.type[iMT][jMT]=ZERO;
	      ovlp.type[jMT][iMT]=ZERO;
	    }
	 

	
      }
  else
    printf("xcompact:> compaction not succeeded in %d=MXITERs\n",MXIT);
  
  
  printf("xcompact:> E min %f \n",emin);
  printf("xcompact:> Finished with %d iterations\n\n",nriter);
  
  free_dvector(xcoor,1,nMT);
}


/* ----------------------------------------------------------*/	  


double potx(double *xcoor)
{

  /* Calculates the total energy of interaction along x-axis 
     between all MTs. The energy consists of strong spring 
     interactions between bounded MTs, (essentially this keeps
     all bounded MTs together), and a LJ-energy between non-bounded
     MTs.
  */
  
  int iMT,jMT;
  double xi,xj,yi,yj,zi,zj,r0,xeq,rm1,rLJ,r,energy;

  
  energy=0.0;
  for (iMT=1;iMT<MT.number;iMT++)
    {
      xi=xcoor[iMT];
      yi=(double)MT.cm[iMT].y;
      zi=(double)MT.cm[iMT].z;
      for (jMT=iMT+1;jMT<=MT.number;jMT++)
	{
	  xj=xcoor[jMT];
	  yj=(double)MT.cm[jMT].y;
	  zj=(double)MT.cm[jMT].z;
	  rLJ=sqrt(pow(0.5*(MT.length[iMT]+MT.length[jMT]),2)
		   +4*exclude*exclude); 
	  /* this is the requiered LJ distance between MTs */
	  r0=rLJ*pow(((double)crep/(double)catt), 
		     1./(double)(catt-crep)); /* this is the LJ sigma param */
	  r=sqrt((xj-xi)*(xj-xi)+(yj-yi)*(yj-yi)+(zj-zi)*(zj-zi))+EPS;
	  rm1=r0/r; /* already scaled */
	  xeq=0.5*(MT.length[iMT]+MT.length[jMT])-
	    (double)overlap(iMT,jMT);
	  	  
	  /* Spring Energy */
	  if (overlap(iMT,jMT)>0)
	    energy += 0.5*kbond*(fabs(xj-xi)-xeq)*(fabs(xj-xi)-xeq);
	  else
	  /* Lennard-Johns */
	    energy += epsilon*(/* pow(rm1,crep) */-pow(rm1,catt));

	}
    }	  

  return(energy);
}


/* ------------------------------------------------------*/


void dpotx(double *xcoor, double *grad)
{

 /* Calculates the gradient of Lennard-Jones potential
    given in pot().
  */
  
  int iMT,jMT;
  double xi,xj,yi,yj,zi,zj,rm1,xeq,r0,rLJ,r;
  
  /* Initialize the gradient vector */  
  for (iMT=1;iMT<=MT.number;iMT++)
    grad[iMT]=0.0;
  
  /* Calculate the gradient vector */
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      xi=xcoor[iMT];
      yi=(double)MT.cm[iMT].y;
      zi=(double)MT.cm[iMT].z;
      for (jMT=1;jMT<=MT.number;jMT++)
	    {
	      if (jMT!=iMT)
		{
		  xj=xcoor[jMT];
		  yj=(double)MT.cm[jMT].y;
		  zj=(double)MT.cm[jMT].z;
		  rLJ=0.5*(MT.length[iMT]+MT.length[jMT]); 
		  /* this is the requiered LJ distance between MTs */
		  r0=rLJ*pow(((double)crep/(double)catt), 
			     1./(double)(catt-crep)); 
		  /* this is the LJ sigma param */
		  r=sqrt((xj-xi)*(xj-xi)+(yj-yi)*(yj-yi)+(zj-zi)*(zj-zi))+EPS;
		  rm1=r0/r; /* already scaled */
		  xeq=0.5*(MT.length[iMT]+MT.length[jMT])-
		    (double)overlap(iMT,jMT);
		  
		  
		  /* Spring */
		  if (overlap(iMT,jMT)>0)
		    grad[iMT] += 0.5*kbond*(xi-xj)*(1-xeq/fabs(xi-xj));
		  else
		  /* Lennard-Johns */
		  grad[iMT] += 0.5*epsilon*
		    (catt*pow(rm1,catt)/* -crep*pow(rm1,crep) */)*
		    (xi-xj)/(r*r);
		}
	    }
    }
  
}
