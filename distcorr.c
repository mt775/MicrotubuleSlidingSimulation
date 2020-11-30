/* ########################################################### */
/* This function calculates polarity correlation function */
/* i.e., C(x|t)=<n(0)n(x)>_{over all MTs}  

#include <stdio.h>
#include <math.h>
#include "global_var.h"

int getjMT(int i,int it);


void timecorr()
{
  int it1,it2,idt,nMT1,nMTs,jMT,iMT,i,npp,nmm,npm;
  float v1,v2,tmpcor_pp,tmpcor_mm,tmpcor_pm,tmpcor;
  

  for (it1=corr_ti;it1<=corr_tf;it1++)
    {
      nMT1=track.nMT[it1];
      for (it2=it1;it2<=corr_tf;it2++)
	{
	  idt=it2-it1;
	  nMTs=0;
	  npp=0;
	  npm=0;
	  nmm=0;
	  tmpcor=0;
	  tmpcor_pp=0;
	  tmpcor_pm=0;
	  tmpcor_mm=0;
	  
	  for (iMT=1;iMT<=nMT1;iMT++)
	    {
	      i=track.index[it1][iMT];
	      jMT=getjMT(i,it2);
	      if (jMT != 0)
		{
		  v1=track.vel[it1][iMT];
		  v2=track.vel[it2][jMT];
		  tmpcor += (v1*v2);
		  nMTs++;

		  if (MT.direct[iMT].x>0 && MT.direct[jMT].x>0)
		    {
		      tmpcor_pp+=(v1*v2);
		      npp++;
		    }
		  else if ((MT.direct[iMT].x>0 && MT.direct[jMT].x<0) 
			   ||  (MT.direct[iMT].x<0 && MT.direct[jMT].x>0)) 
		    {
		      tmpcor_pm+=(v1*v2);
		      npm++;
		    }
		  else if (MT.direct[iMT].x<0 && MT.direct[jMT].x<0)
		    {
		      tmpcor_mm+=(v1*v2);
		      nmm++;
		    }
		  
		}
	    }
	  if (nMTs>0)
	    tmpcor /= (float)nMTs;
	  if (npp>0)
	    tmpcor_pp /= (float)npp;
	  if (nmm>0)
	    tmpcor_mm /= (float)nmm;
	  if (npm>0)
	    tmpcor_pm /= (float)npm;
	  
	  av.velcort[it1][idt] += tmpcor;
	  av.velcort_pp[it1][idt] += tmpcor_pp;
	  av.velcort_pm[it1][idt] += tmpcor_pm;
	  av.velcort_mm[it1][idt] += tmpcor_mm;
	}
    }
}


/*-----------------------------------------------*/
/* find the jMT whose index at iteration it2, is i */  

int getjMT(int i,int it2)
{
  int jMT,j,nMT2;

  jMT=0;
  nMT2=track.nMT[it2];
  do
    {
      jMT++;
      j=track.index[it2][jMT];
      if (jMT>nMT2)
	{
	  jMT=0;
	  break;
	}
    }
  while (j != i);
  
  /* printf("i=%d j=%d index[%d][%d]=%d",i,j,it2,jMT,track.index[it2][jMT]); */
/*   getchar(); */

  return jMT;
}
