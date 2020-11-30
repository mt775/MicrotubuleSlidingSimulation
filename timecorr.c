/* ########################################################### */
/* This function calculates time velocity correlation function */

#include <stdio.h>
#include <math.h>
#include "global_var.h"

int getjMT(int i,int it);


void timecorr()
{
  int it1,it2,idt,nMT1,nMTs,jMT,jMT2,iMT,iMT2,i,npp,nmm,npm;
  float v1,v2,tmpcor_pp,tmpcor_mm,tmpcor_pm,tmpcor;
  

  for (it1=corr_ti;it1<=corr_tf;it1++)
    {
      nMT1=track.nMT[it1];//number of MTs in it1
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
	      i=track.index[it1][iMT];/* get the index of MT number iMT from it1*/
	      iMT2=getjMT(i,it2);  /* this is the same iMT from it1 at time it2 */
	      if (iMT2 != 0)
		{
		  v1=track.vel[it1][iMT];// velocity of filament at it1
		  v2=track.vel[it2][iMT2];// velocity of filament at it2
		  tmpcor += (v1*v2);//multiply the velocity of filament AT later time
		  nMTs++;//for later normalization

		  if (MT.direct[iMT].x>0 && MT.direct[iMT2].x>0)
		    {
		      tmpcor_pp+=(v1*v2);
		      npp++;
		    }
          /*Check if direction of filaments is conserved so that they are still the 
           same.. Behold that track direct runs over all SEGMENTS => 2*iMT-1*/
		  else if ((track.direct[it1][2*iMT-1]>0 && track.direct[it2][2*iMT2-1]<0)
			   ||  (track.direct[it1][2*iMT-1]<0 && track.direct[it2][2*iMT2-1]>0))
		    {
		      printf("timcorr:> logical error, iMT and iMT2 should be the same MT\n");
                printf("%i %i\n", iMT,iMT2);
                printf("index=%i\n",i);
                printf("it1=%i it2=%i\n", it1,it2);
                printf("direct1=%f direct2=%f\n", track.direct[it1][2*iMT-1],track.direct[it2][2*iMT2-1]);
		      getchar();
		      tmpcor_pm+=(v1*v2);
		      npm++;
		      
		    }
		  else if (MT.direct[iMT].x<0 && MT.direct[iMT2].x<0)
		    {
		      tmpcor_mm+=(v1*v2);
		      nmm++;
		    }
		  
		}

	      /* get cross-correlations */
	      for (jMT=1;jMT<nMT1;jMT++)
		{
		  i=track.index[it1][jMT];
		  jMT2=getjMT(i,it2);  /* this is the same 
					  jMT at different time */
		  if (jMT2 !=0) 
		    {
            /*Check if direction of filaments is conserved so that they are still the
            same.. Behold that track direct runs over all SEGMENTS => 2*iMT-1*/
            if ((track.direct[it1][2*iMT-1]>0 && track.direct[it2][2*iMT2-1]<0)
                    ||  (track.direct[it1][2*iMT-1]<0 && track.direct[it2][2*iMT2-1]>0))
			{
			  v1=track.vel[it1][iMT];
			  v2=track.vel[it2][jMT2];
			  tmpcor_pm+=(v1*v2);
			  npm++;
			}
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

/* int getjMT(int i,int it2) */
/* { */
/*   int jMT,j,nMT2; */

/*   jMT=0; */
/*   nMT2=track.nMT[it2]; */
/*   do */
/*     { */
/*       jMT++; */
/*       j=track.index[it2][jMT]; */
/*       if (jMT>nMT2) */
/* 	{ */
/* 	  jMT=0; */
/* 	  break; */
/* 	} */
/*     } */
/*   while (j != i); */
  
/*   /\* printf("i=%d j=%d index[%d][%d]=%d",i,j,it2,jMT,track.index[it2][jMT]); *\/ */
/* /\*   getchar(); *\/ */

/*   return jMT; */
/* } */
