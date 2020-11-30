#include <stdio.h>
#include <math.h>
#include "global_var.h"

void avcm2(struct point *av)
{

/* calculate <xcm^2-xref>, <zcm^2-zref> */

  int iMT;
  float xref;

  av->x=0.0;
  av->z=0.0;

  if (bound.type[0]==ABSORB) xref=lbound0;
  if (bound.type[0]==NON && bound.type[1]==NON) xref=(rbound0-lbound0)/2.;
  
  for (iMT=1;iMT<=MT.number;iMT++)
    av->x += sqrtf((MT.cm[iMT].x-xref)*(MT.cm[iMT].x-xref));
  av->x /= MT.number;

  for (iMT=1;iMT<=MT.number;iMT++)
    av->z += sqrtf(MT.cm[iMT].z*MT.cm[iMT].z);
  av->z /= MT.number;
}

void avvel2(struct point *av)
{
  int iMT;

  /* calculate <vel.x^2>, <vel.z^2> */
  
  av->x=0.0;
  av->z=0.0;

  for (iMT=1;iMT<=MT.number;iMT++)
    av->x += sqrtf(MT.vel[iMT].x*MT.vel[iMT].x);
  av->x /= MT.number;

  for (iMT=1;iMT<=MT.number;iMT++)
    av->z += sqrtf(MT.vel[iMT].z*MT.vel[iMT].z);
  av->z /= MT.number;
}

int iMTxmax()
/* Finds the MT with largest right-x */
{
  int iMT,imax;
  
  imax=1;
  for (iMT=2;iMT<=MT.number;iMT++){
        if ((MT.cm[iMT].x+0.5*MT.length[iMT])>
            (MT.cm[imax].x+0.5*MT.length[imax])) imax=iMT;
        
  }
    
  
  return imax;
}

int iMTxmin()
/* Finds the MT with minimal left-x */
{
  int iMT,imin;
  
  imin=1;
  for (iMT=2;iMT<=MT.number;iMT++){
    if ((MT.cm[iMT].x-0.5*MT.length[iMT])<
	(MT.cm[imin].x-0.5*MT.length[imin])) imin=iMT;
  }
  return imin;
}

int iMTcmmax()
/* Finds the MT with largest xcm */
{
  int iMT,imax;
  
  imax=1;
  for (iMT=2;iMT<=MT.number;iMT++)
    if (MT.cm[iMT].x>MT.cm[imax].x) imax=iMT;
  
    
  return imax;
}

int iMTcmmin()
/* Finds the MT with minimal xcm */
{
  int iMT,imin;
  
  imin=1;
  for (iMT=2;iMT<=MT.number;iMT++)
    if (MT.cm[iMT].x<MT.cm[imin].x) imin=iMT;
  
  return imin;
}

int iMTycmmax()
/* Finds the MT with largest ycm */
{
  int iMT,imax;
  
  imax=1;
  for (iMT=2;iMT<=MT.number;iMT++)
    if (MT.cm[iMT].y>MT.cm[imax].y) imax=iMT;
  
  return imax;
}

int iMTycmmin()
/* Finds the MT with minimal ycm */
{
  int iMT,imin;
  
  imin=1;
  for (iMT=2;iMT<=MT.number;iMT++)
    if (MT.cm[iMT].y<MT.cm[imin].y) imin=iMT;
  
  return imin;
}


int iMTzcmmax()
/* Finds the MT with largest zcm */
{
  int iMT,imax;
  
  imax=1;
  for (iMT=2;iMT<=MT.number;iMT++)
    if (MT.cm[iMT].z>MT.cm[imax].z) imax=iMT;
  
  return imax;
}

int iMTzcmmin()
/* Finds the MT with minimal zcm */
{
  int iMT,imin;
  
  imin=1;
  for (iMT=2;iMT<=MT.number;iMT++)
    if (MT.cm[iMT].z<MT.cm[imin].z) imin=iMT;
  
  return imin;
}


int riMTcmmax()
/* Finds the MT with largest real xcm */
{
  int iMT,imax;
  
  imax=1;
  for (iMT=2;iMT<=MT.number;iMT++)
    if (MT.rcm[iMT].x>MT.rcm[imax].x) imax=iMT;
  
  return imax;
}

int riMTcmmin()
/* Finds the MT with minimal real xcm */
{
  int iMT,imin;
  
  imin=1;
  for (iMT=2;iMT<=MT.number;iMT++)
    if (MT.rcm[iMT].x<MT.rcm[imin].x) imin=iMT;
  
  return imin;
}


int iabsVelmax()
/* Finds the MT with largest absolute velocity */
{
  int iMT,imax;
  
  imax=1;
  for (iMT=2;iMT<=MT.number;iMT++)
    if (fabs(MT.vel[iMT].x/*-velpol*MT.direct[iMT].x*/)>
	fabs(MT.vel[imax].x/*-velpol*MT.direct[imax].x*/)) imax=iMT;
  
  return imax;
}

int iabsVelmin()
/* Finds the MT with minimal absolute velocity */
{
  int iMT,imin;
  
  imin=1;
  for (iMT=2;iMT<=MT.number;iMT++)
    if (fabs(MT.vel[iMT].x/*-velpol*MT.direct[iMT].x*/)<
	fabs(MT.vel[imin].x/*-velpol*MT.direct[imin].x*/)) imin=iMT;
  
  return imin;
}

int iVelmax()
/* Finds the MT with largest velocity */
{
  int iMT,imax;
  
  imax=1;
  for (iMT=2;iMT<=MT.number;iMT++)
    if (MT.vel[iMT].x>MT.vel[imax].x) imax=iMT;
  
  return imax;
}

int iVelmin()
/* Finds the MT with minimal velocity */
{
  int iMT,imin;
  
  imin=1;
  for (iMT=2;iMT<=MT.number;iMT++)
    if (MT.vel[iMT].x<MT.vel[imin].x) imin=iMT;
  
  return imin;
}

int nMTright(float x, float dx)
{
  int iMT,nMT;

  nMT=0;
  for (iMT=1;iMT<=MT.number;iMT++)
    if (MT.direct[iMT].x>0)
      if (MT.cm[iMT].x>x && MT.cm[iMT].x<x+dx)
	nMT++;

  return nMT;
}

int nMTleft(float x, float dx)
{
  int iMT,nMT;

  nMT=0;
  for (iMT=1;iMT<=MT.number;iMT++)
    if (MT.direct[iMT].x<0)
      if (MT.cm[iMT].x>x && MT.cm[iMT].x<x+dx)
	nMT++;

  return nMT;
}


/* ------------------------- Alternative Forms ----------------------------*/
struct point avcm2b()
{
  int iMT;
  struct point av;

  /* calculate <xcm^2-lbound>, <zcm^2-lbound> */
  
  for (iMT=1;iMT<=MT.number;iMT++)
    av.x += sqrtf((MT.cm[iMT].x-lbound0)*(MT.cm[iMT].x-lbound0));
  av.x /= MT.number;

  for (iMT=1;iMT<=MT.number;iMT++)
    av.z += sqrtf(MT.cm[iMT].z*MT.cm[iMT].z);
  av.z /= MT.number;

  
  return av;
}

struct point avvel2b()
{
  int iMT;
  struct point av;

  /* calculate <vel.x^2>, <vel.z^2> */
  
  for (iMT=1;iMT<=MT.number;iMT++)
    av.x += sqrtf(MT.vel[iMT].x*MT.vel[iMT].x);
  av.x /= MT.number;

  for (iMT=1;iMT<=MT.number;iMT++)
    av.z += sqrtf(MT.vel[iMT].z*MT.vel[iMT].z);
  av.z /= MT.number;

  
  return av;
}  

  
  


      
