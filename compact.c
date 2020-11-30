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
#define  crep 8  /* LJ powers: repulsion */
#define  catt 1  /* LJ powers: attraction */
#define  ccen 2 /* power of centering potential */
#define rmin2 0.9*2*exclude 
              /* minimum radius of interaction, 
		 below which there is free sliding */
/* #define  r0 2*exclude*pow((crep/catt),1./(catt-crep))*/
/*#define r0 2*exclude*/
/*#define epsilon 1.0*/ /* LJ energy coef */
/*#define potcen  0.001/pow(exclude,ccen)*/   
#define delta exclude /* step size to move two colliding MTs */

extern void dlinmin(double *p, double *xi, int n, double *fret, 
		    double (*pot)(double *coor));
extern void frprmn(double *p, int n, double ftol, int *iter, double *fret,
		   double (*pot)(double *coor), void (*dpot)(double *coor, 
							     double *grad),
		   int MXIT);

double pot(double *coor);
void dpot(double *coor, double *grad);

/* global variables for use in this file */
int *moveindex,nmove,printgrad;
double r0,potcen,epsilon;


void compact()
{
  /* Scan all MTs. If an MT has NO overlap with any other MT, 
     preform Conjugate Gradient steps tp move it in the yz plane.
  */  
  
  int iy,iz,jy,jz,iMT,jMT,nMT,move;
  int iseg,jseg;
  int nriter,MTcollide;
  double emin,ftol;
  double *coor;
  double *tmpgrad;
  int MXIT,nloop;
  
  double xi1,xi2,yi,zi,xj1,xj2,yj,zj,r,rm1,xov;
  double energy,eij,dy,dz,ycm,zcm,tet;
  
  r0=2*exclude*pow(((double)crep/(double)catt),1./(double)(catt-crep));
  potcen=0.0;  /* Centering potential coef */
  epsilon=1.0/MT.number; /* LJ energy coef */

  printgrad=0;
  nMT=MT.number;
  ftol=1.0e-3; /* parameter for minimization convergaence */
 
  MXIT=500;     /* Max number of minimization steps */
  
  moveindex=ivector(1,nMT);  /* an index to MTs that will be moved */
  
  for (iMT=1;iMT<=nMT;iMT++)
    moveindex[iMT]=0;
  
  /* Find which MTs are not bounded and free to move */
  nmove=0;
  for (iMT=1;iMT<=nMT;iMT++) {
      move=1;
      
      for (jMT=1;jMT<=MT.number;jMT++){
          if (overlap(iMT,jMT)>0) move=0; /* found a neighbor */
      }
      
      if (move==1){
	  /* set nmove to be the n'th MT that moves */ 
          nmove++;
          moveindex[iMT]=nmove;
	  }
  }
  /* nmove is now the number of MTs that move */
  
  
  /* allocate memory for tmp coor vector */
  /* coor == {y1,y2, ..., ynmove, z1,z2, ..., znmove} */
  coor=dvector(1,2*nmove);
  for (iMT=1; iMT<=2*nmove; iMT++) {
      coor[iMT]=0.;
  }
  tmpgrad=dvector(1,2*nmove);
  
  
  /* set the coordinates vector of those MTs that move */
  for (iMT=1;iMT<=nMT;iMT++) {
      if( moveindex[iMT]>0){
          iy=moveindex[iMT];
          iz=moveindex[iMT]+nmove;
	  
          coor[iy]=(double)MT.cm[iMT].y;
          coor[iz]=(double)MT.cm[iMT].z;
	  }
  }

/* calculate center of mass and shift cm to center of axis */
  ycm=0.0;
  zcm=0.0;
  for(iMT=1;iMT<=nMT;iMT++){
      iy=moveindex[iMT];
      iz=moveindex[iMT]+nmove;
      ycm += coor[iy];
      zcm += coor[iz];
  }
  ycm /= (double)nMT;
  zcm /= (double)nMT;
  for(iMT=1;iMT<=nMT;iMT++){
      iy=moveindex[iMT];
      iz=moveindex[iMT]+nmove;
      coor[iy]-=ycm;
      coor[iz]-=zcm;
  }
	  
  
  /* Check that no two MTs overlap. If they do, move each 
     of them randomely by a small Dy, Dz */
  
  nloop=0;
  do
    {
      MTcollide=0;
      nloop += 1;

      /*UN-COMMENT FOR DEBUG */
/*       if (nloop>1) */
/* 	{ */
/* 	  printf("colliding loop: %d\n",nloop); */
/* 	  getchar(); */
/* 	} */

      /* Check if anytwo MT segments collide */
      for (iMT=1;iMT<=nMT;iMT++)
	{
	  if (moveindex[iMT]>0)
	    {
	      iy=moveindex[iMT];
	      iz=moveindex[iMT]+nmove;
	      
	      yi=coor[iy];
	      zi=coor[iz];
	      for (jMT=1;jMT<=nMT;jMT++)
		{
		  jy=moveindex[jMT];
		  jz=moveindex[jMT]+nmove;
		  
		  
		  yj=coor[jy];
		  zj=coor[jz];
		  
		  xov=xoverlap(iMT,jMT);

		  r=sqrt((yj-yi)*(yj-yi)+(zj-zi)*(zj-zi));
		  
		  if (xov>0 && r<rmin2 && iMT!=jMT)
		    {
		      MTcollide = 1;
		      break;
		    }
		}
	      if (MTcollide) break;
	    }
	}
      

      /* Expand the bundle to facilitate separation of colliding 
	 MTs */
      if (MTcollide)
	{
	  potcen=0.01/pow(exclude,ccen); /* initiate centering potential */
	  
	  for(iMT=1;iMT<=nMT;iMT++)
	    {
	      /* don't move MTs that have no x-overlap 
		 with anyother MT, it will drift away */

	      move=0;
	      for (jMT=1;jMT<=MT.number;jMT++)
		{
		  xov=xoverlap(iMT,jMT);
		  if (xov>0) move=1;
		}
	      if (move)
		{
		  iy=moveindex[iMT];
		  iz=moveindex[iMT]+nmove; 
		  dy = 1.0*coor[iy]; 
		  dz = 1.0*coor[iz];
		  coor[iy] += dy;
		  coor[iz] += dz;
		}
	    }

	  /* Shift still colliding MTs in the expanded bundle */
	    for (iMT=1;iMT<=nMT;iMT++)
	    {
	      if (moveindex[iMT]>0)
		{
		  iy=moveindex[iMT];
		  iz=moveindex[iMT]+nmove;
		  
		  yi=coor[iy];
		  zi=coor[iz];
		  for (jMT=1;jMT<=nMT;jMT++)
		    {
		      jy=moveindex[jMT];
		      jz=moveindex[jMT]+nmove;
		      
		      yj=coor[jy];
		      zj=coor[jz];
		      
		      xov=xoverlap(iMT,jMT);
		      r=sqrt((yj-yi)*(yj-yi)+(zj-zi)*(zj-zi));
		      
		      
		      if (xov>0 && r<rmin2 && iMT!=jMT)
			{
			  MTcollide = 2;
			  if (r==0.0)
			    {
			      tet=2*PI*ran1(&idum);
			      dy=2*exclude*cos(tet);
			      dz=2*exclude*sin(tet);
			      coor[jy]=dy;
			      coor[jz]=dz;
			      /*printf("compact:> MTs [%d,%d] collide, r=%f, tet=%f\n",
				MT.index[iMT],MT.index[jMT],r,tet*180/PI);*/
			      /* getchar();  */
			    }
			  else
			    {
			      dy=2.0*exclude*(yj-yi)/r;
			      dz=2.0*exclude*(zj-zi)/r;
			      coor[jy]=coor[iy]+dy;
			      coor[jz]=coor[iz]+dz;
			      /*printf("compact:> MTs [%d,%d] collide, r=%f\n",
				MT.index[iMT],MT.index[jMT],r);*/
			      /* getchar();  */
			    }

		    }
		}
	    }
	}
	}
    }
  while (MTcollide==2);
  
  /* printf("compact:> Finished Colliding iterations: %d\n",nloop);  */ 
  
  /* Apply NR Congugate Gradient routine Ch. 10.6; P. 423*/
  frprmn(coor,2*nmove,ftol,&nriter,&emin,pot,dpot,MXIT);
 


  /* Store the new coordinates */
  if (nriter<MXIT)
    {
     for (iMT=1;iMT<=nMT;iMT++)
	{
	  if (moveindex[iMT]>0)
	    {
	      iy=moveindex[iMT];
	      iz=moveindex[iMT]+nmove;
	      /* printf("compact:>Moved iMT %d zi=%f --> zf=%f\n", */
/* 		     MT.index[iMT],MT.cm[iMT].z,coor[iz]); */
/* 	      printf("compact:>             yi=%f --> yf=%f\n", */
/* 		     MT.cm[iMT].y,coor[iy]); */
/* 	      printf("compact:>iy= %d iz=%d\n",iy,iz); */
	      MT.cm[iMT].y=(float)coor[iy];
	      MT.cm[iMT].z=(float)coor[iz];
	      MT.pend[iMT].y=(float)coor[iy];
	      MT.pend[iMT].z=(float)coor[iz];
	      MT.mend[iMT].y=(float)coor[iy];
	      MT.mend[iMT].z=(float)coor[iz];
	    }
	}
     for (iseg=1;iseg<=2*nMT-1;iseg++)
       for(jseg=iseg;jseg<=2*nMT;jseg++)
	 if (soverlap(iseg,jseg)==0)
	   {
	     ovlp.type[iseg][jseg]=ZERO;
	     ovlp.type[jseg][iseg]=ZERO;
           ovlp.motor_direction[iseg][jseg]=0;
           ovlp.motor_direction[jseg][iseg]=0;
	   }
     
    }
  else
    {
      printf("compact:> compaction not succeeded in %d=MXITERs\n",MXIT);
      printf("compact:> Time iteration step: iter=%d  iav=%d\n",iter,iav);
    }

  /* printf("compact:> Time iteration step: iter=%d  iav=%d\n",iter,iav); */
/*   printf("compact:> Moved %d MTs\n",nmove); */
/*   printf("compact:> Colliding iterations: %d\n",nloop); */
/*   printf("compact:> E min %f \n",emin); */
/*   printf("compact:> Finished with %d iterations\n\n",nriter); */
  /*getchar();*/  

  free_ivector(moveindex,1,nMT);
  free_dvector(coor,1,2*nmove);
  free_dvector(tmpgrad,1,2*nmove);

}


/* ----------------------------------------------------------*/	  


double pot(double *coor)
{

  /* Calculates the total Lennard-Jones energy of interaction
     between all pairs of MTs. The energy depends linearly 
     on the overlap in the x direction.   
  */
  
  int iMT,jMT,nMT,iy,iz,jy,jz;
  double xi,yi,zi,xj,yj,zj,r,rm1;
  double energy,eij,xov;

  

  /* This makes Distance at LJ-energy minimum at 2*exclude */  
  

  energy=0.0;
  for (iMT=1;iMT<MT.number;iMT++)
    {
      if (moveindex[iMT]>0)
	{
	  iy=moveindex[iMT];
	  iz=moveindex[iMT]+nmove;
	  yi=coor[iy];
	  zi=coor[iz];
	}
      else
	{
	  yi=(double)MT.cm[iMT].y;
	  zi=(double)MT.cm[iMT].z;
	}
      for (jMT=iMT+1;jMT<=MT.number;jMT++)
	{
	  if (moveindex[jMT]>0)
	    {
	      jy=moveindex[jMT];
	      jz=moveindex[jMT]+nmove;
	      yj=coor[jy];
	      zj=coor[jz];
	    }
	  else
	    {
	      yj=(double)MT.cm[jMT].y;
	      zj=(double)MT.cm[jMT].z;
	    }
	  
	  
	  r=sqrt(pow((yi-yj),2)+pow((zi-zj),2))+EPS;	  
	  rm1=r0/r; /* already scaled */
	  xov=xoverlap(iMT,jMT);
	  
	  /* Lennard-Johns */
	  eij=epsilon*(pow(rm1,crep)-pow(rm1,catt))*xov;
	  
	  
	  energy += eij;
	  /* if(fabs(eij)>10000) */
/* 	    { */
/* 	      printf("pot:>LJ energy e[%d,%d]=%f\n",MT.index[iMT],MT.index[jMT] */
/* 		     ,eij); */
/* 	      printf("pot:>R[%d,%d]=%f\n",MT.index[iMT],MT.index[jMT] */
/* 		     ,r); */
/* 	      getchar(); */
/* 	    } */

	}
    }

  /*  Add a centering potential to all MTs */  
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      if (moveindex[iMT]>0)
	{
	  iy=moveindex[iMT];
	  iz=moveindex[iMT]+nmove;
	  yi=coor[iy];
	  zi=coor[iz];
	}
      else
	{
	  yi=MT.cm[iMT].y;
	  zi=MT.cm[iMT].z;
	}
      r=sqrt(yi*yi+zi*zi);
      energy += potcen*pow(r,ccen);
    }
  
  /*energy/=(double)MT.number;*/
  return(energy);
}


/* ------------------------------------------------------*/


void dpot(double *coor, double *grad)
{

 /* Calculates the gradient of Lennard-Jones potential
    given in pot().
  */
  
  int iMT,jMT,nMT,iy,iz,jy,jz;
  double xi,yi,zi,xj,yj,zj,r,rm1,xov,gy,gz;
  
  /* Initialize the gradient vector */  
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      if (moveindex[iMT]>0)
	{
	  iy=moveindex[iMT];
	  iz=moveindex[iMT]+nmove;
      	  grad[iy]=0.0;
	  grad[iz]=0.0;
	}
    }
  
  
  /* Calculate the gradient vector */
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      if (moveindex[iMT]>0)
	{
	  iy=moveindex[iMT];
	  iz=moveindex[iMT]+nmove;
	  yi=coor[iy];
	  zi=coor[iz];
	  	  
	  for (jMT=1;jMT<=MT.number;jMT++)
	    {
	      if (jMT!=iMT)
		{
		  if (moveindex[jMT]>0)
		    {
		      jy=moveindex[jMT];
		      jz=moveindex[jMT]+nmove;
		      yj=coor[jy];
		      zj=coor[jz];
		    }
		  else
		    {
		      yj=(double)MT.cm[jMT].y;
		      zj=(double)MT.cm[jMT].z;
		    }
		  
		  r=sqrt(pow((yi-yj),2)+pow((zi-zj),2))+EPS;
		  rm1=r0/r; /* already scaled */
		  xov=xoverlap(iMT,jMT);

		  if (printgrad)
		    {
		      printf("dpot:> r[%d,%d]=%f r0=%f rm1=%f\n",
			     MT.index[iMT],MT.index[jMT],r,r0,rm1);
		    }
		  
		  /* Lennard-Johns */
		      gy=epsilon*xov*
			(catt*pow(rm1,catt)-crep*pow(rm1,crep))
			*(yi-yj)/(r*r);
			/*/(double)MT.number*/
		      gz=epsilon*xov*
			(catt*pow(rm1,catt)-crep*pow(rm1,crep))
			*(zi-zj)/(r*r);
			/*/(double)MT.number*/
		      if (printgrad)
			{
			  printf("dpot:> ov[%d,%d]=%f\n",
				 MT.index[iMT],MT.index[jMT],
				 xoverlap(iMT,jMT));
			  printf("dpot:> gy[%d],%d=%f\n",
				 MT.index[iMT],MT.index[jMT],gy);
			  printf("dpot:> gz[%d],%d=%f\n",
				 MT.index[iMT],MT.index[jMT],gz);
			  
			}
		      
		      grad[iy] += 0.5*gy; /* half because we count each 
					     interaction twice */
		      grad[iz] += 0.5*gz;
		}
	    }
	}
    }
  
  
  if (printgrad)
    {
      printf("dpot:>--------------------------------------\n");
      for (iMT=1;iMT<=MT.number;iMT++) 
	{
	  printf("moveindex[%d]=%d\n",iMT,moveindex[iMT]);
	  iy=moveindex[iMT];
	  iz=moveindex[iMT]+nmove; 
	  printf("dpot:> grad[%d]y=%f\n",MT.index[iMT],grad[iy]);
	  printf("dpot:> grad[%d]z=%f\n",MT.index[iMT],grad[iz]);
	  printf("dpot:> y[%d]=%f  z[%d]=%f\n",MT.index[iMT],
		 coor[iy],MT.index[iMT],coor[iz]);
	}    
      getchar();
    }
  
  /* Add gradient due to centering potential */ 
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      if (moveindex[iMT]>0)
	{
	  iy=moveindex[iMT];
	  iz=moveindex[iMT]+nmove;
	  yi=coor[iy];
	  zi=coor[iz];
  
	  r=sqrt(yi*yi+zi*zi);
	  
	  grad[iy] += potcen*ccen*pow(r,ccen-2)*yi;
	  grad[iz] += potcen*ccen*pow(r,ccen-2)*zi;

	  if (printgrad)
	    {
	      printf("dpot:> grad[%d]y=%f\n",MT.index[iMT],grad[iy]);
	      printf("dpot:> grad[%d]z=%f\n",MT.index[iMT],grad[iz]);
	      printf("dpot:> y[%d]=%f  z[%d]=%f\n",MT.index[iMT],
		     coor[iy],MT.index[iMT],coor[iz]);
	      getchar();
	    }
 	}
    }
}

