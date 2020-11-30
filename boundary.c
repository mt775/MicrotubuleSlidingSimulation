/* ################################################################### */
/* boundary -- preformes boundary conditions and fills the struct seg. */
/*             This routine is called twice in each time iteration:    */
/*             1. by mkMT(), to initiate the segment data in accord    */
/*                with old and new MT. data. Segment data is used      */
/*                in equations(), and prior to that, in overlapMap().  */
/*             2. by newpos, to apply absorbing b.c. by setting        */
/*                the absorbed MTs' length to zero.                    */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "global_var.h"
#include "ran1.h"
#include <stdlib.h>

/*This function calculates and prints: (1) ndel, i.e. the number of MTs marked
for deletion (2) adjustedDegprob, i.e. the probability of MT degradation,
adjusted by ndel; (3) n_actually_degraded, i.e. the subset of ndel that is
actually degraded.*/

/*MJ: General comments
Please watch out for tabs and curly brackets. See below how I changed it
*/
void MT_degrade(){
  int iMT;
  int ndel=0;
  int n_actually_degraded=0;
  float adjustedDegprob=0.;
  FILE *fp;

  
  if (DEBUG==1){
    //Using filept here to write to OCHECK or given output directory
    fp = fopen(filept("Debug_MT_degrade.dat"),"a");
    /*Flag up inability to open*/
    if (fp == NULL){
      printf("Error opening file.\n");
     exit(1);
    }
  }
  /********************The following section can be used to fix the outflux. Not necessary ATM********/
  /*Set up counter, ndel, for number of MTs touching left boundary*/
  //for (iMT=1; iMT<=MT.number; iMT++){
  //  if (MT.cm[iMT].x-lbound0 < MT.length[iMT]/2.0){
	//    ndel++;
  //  }
  //}
  /*If >=1 MTs are touching the left boundary, count how many of these
  are actually deleted*/
  // if (ndel>0){
  //   adjustedDegprob=degprob/((float) ndel);
  //   for (iMT=1; iMT<=MT.number; iMT++){
  //     if (ran1(&idum)< adjustedDegprob && MT.cm[iMT].x-lbound0 < MT.length[iMT]/2.0){
	//       MT.length[iMT]=0.;
  //       n_actually_degraded++;
  //       if (DEBUG==1)
	//         printf("Deleted MTcm %f degprob %f touching boud %i actually degraded %i \n",MT.cm[iMT].x,adjustedDegprob,ndel,n_actually_degraded);
  //     }
  //   }
  // }

  for (iMT=1; iMT<=MT.number; iMT++){
    if (ran1(&idum)< degprob && MT.cm[iMT].x-lbound0 < MT.length[iMT]/2.0){
	    MT.length[iMT]=0.;
      n_actually_degraded++;
      if (DEBUG==1)
	      printf("Deleted MTcm %f \n",MT.cm[iMT].x);
    }
  }

  if (DEBUG==1){
    //printf("%i\n",iter);
    fprintf(fp,"%i\t%i\n", iter, n_actually_degraded);
    fclose(fp);
  }
}

void boundary()
{
  int iMT,nMT,iseg1,iseg2,jseg;
  float xcm,lx,rx,l2,direct,rb,lb;

  nMT=MT.number;

  /* For the most generic case that MTs don't
     hit any boundary we first set all first segment data
     identical with MT data; the effectivness of the second
     segment is eliminated by setting its length to zero.
  */

 /*random MT degration*/
  if (degprob>0.){
    MT_degrade();
  }      

  /* Apply Boundary conditions                          */
  /* and change what needs to be changed from the above */

  for (iMT=1;iMT<=nMT;iMT++)
    {
      xcm=MT.cm[iMT].x;
      l2=0.5*MT.length[iMT];
      direct=MT.direct[iMT].x;
      lx=xcm-l2;
      rx=xcm+l2;

      /* Absorbing boundary conditions on either left/right boundaries */
      /*            lb0                   rb0           */
      /* <-----------|                     |            */
      /*             |                     |----------->*/
      if ((bound.type[0]==ABSORB && rx<=lbound0) ||
	  (bound.type[1]==ABSORB && lx>=rbound0))
	{
	  MT.length[iMT]=0.0; /* eliminate iMT next iteration */
	}

        


      //under STICKY and NET bc MTs may fall into cell body also helps preventing "crazy" MTs flying around
        if ( ( bound.type[0]==STICKY || bound.type[0]==NET || bound.type[0]==CONFLX || bound.type[0]==WALL) && rx<=lbound0){
                MT.length[iMT]=0.0;
        }

      /* Popup boundary condition on either left/right boundaries */
      /* is similar to absorbing only the MT leaving the boundary is
	 not eliminated but rather its xcm is randomely chosen
	 inside the domain as if it poped up inside. This is a technical
	 trick to cause prevent new MTs instead of the old ones. */
      /*            lb0                   rb0           */
      /* <-----------|                     |            */
      /*             |                     |----------->*/
      if ((bound.type[0]==POPUP && rx<lbound0) ||
	  (bound.type[1]==POPUP && lx>rbound0))
	{
	  xcm=ran1(&idum)*(rbound0-lbound0)+lbound0;

	  /* store the new position */
	  /* note the real coordinates, MT.rcm[iMT].x remain unchanged! */
	  MT.cm[iMT].x=xcm;
	  MT.mend[iMT].x=xcm+direct*l2;
	  MT.pend[iMT].x=xcm-direct*l2;
	}

      /* Periodic bounadry conditions:                        */
      /* Shift boundary-protruding MT into the doamin         */
      else if (bound.type[0]==PERIODIC && bound.type[1]==PERIODIC
	       && (xcm<lbound0 || xcm>rbound0))
	{
	  if (MT.maxL>=(rbound0-lbound0))
	    {
	      printf("boundary:> using periodic boundary cond.\n");
	      printf("           error MT length must be smaller\n");
              printf("           than domain size rbound0-lbound0\n");
	      printf("           MaxL=%f rbound0=%f, lbound0=%f\n",
		     MT.maxL,rbound0,lbound0);
	      exit(0);
	    }
	  /* shift the MT into the domain */
	  /* case A */
	  while(xcm<lbound0)
	    xcm += (rbound0-lbound0);
	  /* case B */
	  while (xcm>rbound0)
	    xcm -= (rbound0-lbound0);

	  /* now MT.cm should be between lb0 and rb0    */
	  if (xcm <lbound0 || xcm>rbound0)
	    {printf("boundary:> Algoritmic error 2\n"); exit(0);}

	  /* store the new position */
	  MT.cm[iMT].x=xcm;
	  MT.mend[iMT].x=xcm+direct*l2;
	  MT.pend[iMT].x=xcm-direct*l2;
	}

      else if (bound.type[1]==REFLECT || bound.type[1]==RFIXREF)
	{
	  float reflectb;
	  int tmpij,tmpji;
	  if (bound.type[1]==REFLECT)
	    reflectb=(iter>1) ? MT .xmax : rbound0;
	  else if (bound.type[1]==RFIXREF)
	    reflectb=reflectb0;

	  if (xcm>reflectb)
	    {
	      /* if xcm crosses the reflection boundary,
		 invert the MT. This means: shifting the xcm
		 to the mirror location on the other side of
		 the boundary, inverting the MT polarity,
		 switching the overlap segment types. */

	      xcm=reflectb-(xcm-reflectb);
	      MT.cm[iMT].x=xcm;
	      MT.direct[iMT].x = -direct;
	      MT.mend[iMT].x=xcm-direct*l2;
	      MT.pend[iMT].x=xcm+direct*l2;
	      iseg1=2*iMT-1;
	      iseg2=2*iMT;
	      for (jseg=1;jseg<=2*nMT;jseg++)
		{
		  if (jseg!=iseg1 && jseg!=iseg2)
		    {
		      tmpij=ovlp.type[iseg1][jseg];
		      ovlp.type[iseg1][jseg]=ovlp.type[iseg2][jseg];
		      ovlp.type[iseg2][jseg]=tmpij;
		      ovlp.type[jseg][iseg1]=ovlp.type[iseg1][jseg];
		      ovlp.type[jseg][iseg2]=ovlp.type[iseg2][jseg];
		    }
		}
	    }

	}
    }
}


/* ------------------------------------------------------- */

void mkseg(int iMT)
{
  int iseg1,iseg2;
  float xcm,lx,rx,l2,direct,rb,lb;

  /* For the most generic case that MTs don't
     hit any boundary we first set all first segment data
     identical with MT data; the effectivness of the second
     segment is eliminated by setting its length to zero.
  */
  iseg1=2*iMT-1; //uneven numbers
  iseg2=2*iMT;   //even numbers

  seg.length[iseg1]=MT.length[iMT];

  seg.cm[iseg1].x=MT.cm[iMT].x;
  seg.mend[iseg1].x=MT.mend[iMT].x;
  seg.pend[iseg1].x=MT.pend[iMT].x;
  seg.direct[iseg1].x=MT.direct[iMT].x;

  seg.cm[iseg1].y=MT.cm[iMT].y;
  seg.mend[iseg1].y=MT.cm[iMT].y;
  seg.pend[iseg1].y=MT.cm[iMT].y;
  seg.direct[iseg1].y=0.0;

  seg.cm[iseg1].z=MT.cm[iMT].z;
  seg.mend[iseg1].z=MT.cm[iMT].z;
  seg.pend[iseg1].z=MT.cm[iMT].z;
  seg.direct[iseg1].z=0.0;

  seg.length[iseg2]=0.0;

  seg.cm[iseg2].x=MT.cm[iMT].x;
  seg.mend[iseg2].x=MT.mend[iMT].x;
  seg.pend[iseg2].x=MT.pend[iMT].x;
  seg.direct[iseg2].x=MT.direct[iMT].x;

  seg.cm[iseg2].y=MT.cm[iMT].y;
  seg.mend[iseg2].y=MT.cm[iMT].y;
  seg.pend[iseg2].y=MT.cm[iMT].y;
  seg.direct[iseg2].y=0.0;

  seg.cm[iseg2].z=MT.cm[iMT].z;
  seg.mend[iseg2].z=MT.cm[iMT].z;
  seg.pend[iseg2].z=MT.cm[iMT].z;
  seg.direct[iseg2].z=0.0;


  /* Apply Boundary conditions                          */
  /* and change what needs to be changed from the above */

  xcm=MT.cm[iMT].x;
  l2=0.5*MT.length[iMT];
  direct=MT.direct[iMT].x;
  lx=xcm-l2;
  rx=xcm+l2;



  /* Periodic bounadry conditions:                                 */
  /* ||      L-IMAGE    ||    DOMAIN       ||    R-IMAGE      ||   */
  /*   -----------------  -----------------  -----------------     */
  /* Case A.1                                                      */
  /* || <----*--|--     ||         |       ||       |       ||     */
  /* =>                                                            */
  /* ||         |       || <----*--|--     ||       |       ||     */
  /*      A.2                                                      */
  /* ||         |<----*-||---      |       ||       |       ||     */
  /* =>                                                            */
  /* ||         |       ||         |<----*-||---    |       ||     */
  /* Case B.1                                                      */
  /* ||         |       ||         |   <---||-*---- |       ||     */
  /* =>                                                            */
  /* ||         |   <---||-*----   |       ||       |       ||     */

  if (bound.type[0]==PERIODIC && bound.type[1]==PERIODIC)
    {
      /* MT.cm should be between lb0 and rb0    */
      if (xcm <lbound0 || xcm>rbound0) // if center out of bounds => error!
	{printf("boundary:> Algoritmic error 2\n"); exit(0);}

      /* generate equivalent segments      */

      //FIRST left boundary
      if ((xcm-lbound0)<l2) /* eqv. lx<lb0 */

	/*     lb0   xcm                       rb0  */
	/* <-----|-----*---------                |   */

	/* =>                                       */

	/*      |<-----*---------          <-----|   */
	/*      |     iseg1                 iseg2|   */

	{
        //printf("SEGMENTATION\n");
	  seg.length[iseg1]=l2+xcm-lbound0;
	  seg.cm[iseg1].x=lbound0+0.5*seg.length[iseg1];
	  seg.mend[iseg1].x=seg.cm[iseg1].x+0.5*direct*seg.length[iseg1];
	  seg.pend[iseg1].x=seg.cm[iseg1].x-0.5*direct*seg.length[iseg1];
	  seg.direct[iseg1].x=direct;

	  seg.length[iseg2]=l2-(xcm-lbound0);
	  seg.cm[iseg2].x=rbound0-0.5*seg.length[iseg2];
	  seg.mend[iseg2].x=seg.cm[iseg2].x+0.5*direct*seg.length[iseg2];
	  seg.pend[iseg2].x=seg.cm[iseg2].x-0.5*direct*seg.length[iseg2];
	  seg.direct[iseg2].x=direct;


	  /* The next lines perturb the complete periodic b.c */
	  /* The MT environment felt by iseg1 of each MT, is made
	     different from that felf by the iseg2. We make the
             iseg2 environemnt similar to that of some other MT,
             not necessary the one physically close to it.
             Note, by shifting only the z and or y axis
	     keeps the correct length of the segment. Otherwise, at high
             densities we find the MTs to behave differently near the
	     boundaries (they condense). */

	  if (iMT>1)
	    {
	      /* seg.cm[iseg2].y=seg.cm[2*(iMT-1)].y;  */
/* 	      seg.mend[iseg2].y=seg.cm[iseg2].y; */
/* 	      seg.pend[iseg2].y=seg.cm[iseg2].y; */
/* 	      seg.direct[iseg2].y=0.0; */

	      seg.cm[iseg2].z=seg.cm[2*(iMT-1)].z;
	      seg.mend[iseg2].z=seg.cm[iseg2].z;
	      seg.pend[iseg2].z=seg.cm[iseg2].z;
	      seg.direct[iseg2].z=0.0;
	    }
	}

	  //SECOND right boundary
      else if ((rbound0-xcm)<l2) /* eqv. rx>rb0 */

	/*     lb0                          rb0          */
	/*      |                --------*----|----->    */

	/* =>                                            */

	/*      |----->          ------------>|          */
	/*      |iseg1             iseg2      |          */

	{
	  seg.length[iseg2]=l2+rbound0-xcm;
	  seg.cm[iseg2].x=rbound0-0.5*seg.length[iseg2];
	  /* seg.cm[iseg2].x=lbound0+(rbound0-lbound0)*ran1(&idum); */
	  seg.mend[iseg2].x=seg.cm[iseg2].x+0.5*direct*seg.length[iseg2];
	  seg.pend[iseg2].x=seg.cm[iseg2].x-0.5*direct*seg.length[iseg2];
	  seg.direct[iseg2].x=direct;

	  /* make image segment to be unconnected */
	  /* put it somewhere isolated            */

	  if (iMT>1)
	    {
	      /* seg.cm[iseg2].y=seg.cm[2*(iMT-1)].y;  */
/* 	      seg.mend[iseg2].y=seg.cm[iseg2].y; */
/* 	      seg.pend[iseg2].y=seg.cm[iseg2].y; */
/* 	      seg.direct[iseg2].y=0.0; */

	      seg.cm[iseg2].z=seg.cm[2*(iMT-1)].z;
	      seg.mend[iseg2].z=seg.cm[iseg2].z;
	      seg.pend[iseg2].z=seg.cm[iseg2].z;
	      seg.direct[iseg2].z=0.0;
	    }
	  seg.length[iseg1]=l2-(rbound0-xcm);
	  seg.cm[iseg1].x=lbound0+0.5*seg.length[iseg1];
	  seg.mend[iseg1].x=seg.cm[iseg1].x+0.5*direct*seg.length[iseg1];
	  seg.pend[iseg1].x=seg.cm[iseg1].x-0.5*direct*seg.length[iseg1];
	  seg.direct[iseg1].x=direct;

	}
    }

  else if (bound.type[1]==REFLECT || bound.type[1]==RFIXREF)
    /* Apply reflecting boundary condition on right boundary */
    /*     lb0                          rb0   rx  */
    /*      |               --------*----|----->  */

    /* =>                                         */

    /*      |                   iseg1    |        */
    /*      |               -------*---->|        */
    /*      |                      <-----|        */
    /*      |                       iseg2         */
    {
      float reflectb;
      if (bound.type[1]==REFLECT)
	reflectb=(iter>1) ? MT .xmax : rbound0;
      else if (bound.type[1]==RFIXREF)
	reflectb=reflectb0;

      if (rx > reflectb)
	{
	  /* Reflection boundary is at max xcm, not max rx */
	  /* recall, MT.xmax is set to largest xcm in mkMTs(). */
	  seg.length[iseg1]=l2+reflectb-xcm;
	  seg.cm[iseg1].x=reflectb-0.5*seg.length[iseg1];
	  seg.direct[iseg1].x=direct;
	  seg.mend[iseg1].x=seg.cm[iseg1].x+0.5*direct*seg.length[iseg1];
	  seg.pend[iseg1].x=seg.cm[iseg1].x-0.5*direct*seg.length[iseg1];


	  seg.length[iseg2]=l2-(reflectb-xcm);
	  seg.cm[iseg2].x=reflectb-0.5*seg.length[iseg2];
	  seg.direct[iseg2].x=-direct;
	  seg.mend[iseg2].x=seg.cm[iseg2].x-0.5*direct*seg.length[iseg2];
	  seg.pend[iseg2].x=seg.cm[iseg2].x+0.5*direct*seg.length[iseg2];
	}
    }

  if (fabs(seg.length[iseg1]+seg.length[iseg2]-MT.length[iMT]) > 1e-1)
    {
      printf("mkseg:> Error segment length does not sum to MT length\a\n");
      printf("mkseg:> LMT[%d]=%f, lseg[%d.%d]=%f lseg[%d.%d]=%f\n",
	     MT.index[iMT],MT.length[iMT],MT.index[iMT],iseg1-2*iMT+2,
	     seg.length[iseg1],MT.index[iMT],iseg2-2*iMT+2,
	     seg.length[iseg2]);
      getchar();
    }
}
