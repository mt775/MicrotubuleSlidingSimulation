/* #################################################################### */
/* Generate Assembly of MTs                                             */
/* This function initializes and updates the struct MT                  */
/* It also allocates new memory for the vectors used in each iteration  */
/* This is doen in two stages:                                          */
/* 1. Data on the newly positioned array of MTs is organized
      in a new and memory-extended struct                               */
/* 2. New MTs are added                                                 */


#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include "global_var.h"
#include "ran1.h"
#include "nrutil.h"
#include "sampleExpDist.h"
#include "lambert_w.h"
#define MXTRIAL 1000



/* Global Variables for this file */
int nold;
/*Functions*/
void getcoor(int istart, int n); //just a forward definition
void allocate_tmp( struct MTdata *MTtmp, struct OVLPdata *ovlptmp);
void MT_to_tmp( struct MTdata *MTtmp, int *index, int nold );
void ovlp_to_ovlptmp( struct OVLPdata *ovlptmp, int *index, int nold );
void free_old_MTdat();
void allocate_MTdata(int n);
void fill_old_MTdata(int nold, struct MTdata *MTtmp);
void fill_old_OVLPdata( struct OVLPdata *ovlptmp, int nold );
void freetmp( struct MTdata *MTtmp, struct OVLPdata *ovlptmp , int nold , int *index );


void mkMTs(int newMTs)
{
  int i,j,i1,i2,j1,j2,iMT,jMT,n,istart,newgen;
  float ranexp;
  struct OVLPdata ovlptmp;
  struct MTdata MTtmp;
  struct segdata segtmp;
  int *index;
  int nlost;
  int iseg1,iseg2,jseg1,jseg2;



 /* newMTs -- is a new source of MTs to the system */
 /* nlost  -- is the number of MTs that were lost in
	      the previous iteration.              */
 /* nold   -- the number of valid MTs from the previous
	      iteration.                           */
 /* newgen -- is the number of MTs for which the coor
	      need to be regenerated.
	      For fixed total flux rate we have
              newMTs=nlost+nadded;
              for fixed incoming flux:
              newMTs=nadded. */

  nold=0;
  nlost=0; //just setting zero


  /*If it>1 we need to rearrange all data arrays to take
   lost/gained MTs into account*/
  if (iter>1) {
      /* count the number of valid MTs; see boundary.c */
      for (iMT=1;iMT<=MT.number;iMT++)
          if (MT.length[iMT]>0) // maybe boundary sets length zero
              nold++;// only BC surviving MTs are counted
          else{
              nlost++;
          }
      index=ivector(1,nold);
      /*Create tmp array to save all MT data*/
      /*#############################*/
      allocate_tmp( &MTtmp , &ovlptmp);
      /*#############################*/

      /* Generate a pointer to the old indexing. Note that
       only MTs with length>0 "survive"*/
      i=0;
      for (iMT=1;iMT<=MT.number;iMT++){
          if (MT.length[iMT]>0){
              i++;
              index[i]=iMT;
          }
	        else{
		        printf("Deleted MT: %i with cm-lb= %f\n", iMT,MT.cm[iMT].x-lbound0);
	        }
	    }


      /* Store MT array in above defined MTtmp array
       STORING BECAUSE NEW MT ENLARGE THE ARRAY! C DOES NOT SUPPORT DYNAMIC ARRAYS!
       MTs that in their new position cross the
         boundary are eliminated, because nold only counts microtubule with l>0*/

      /*################################################*/
      MT_to_tmp( &MTtmp, index, nold );
      /*################################################*/

      /* Store data on Motor orientation in each overlaping region in temporary ovlp;
         Store overlapping map. Update the overlap map by the new indexing */
      /*###############################################*/
      ovlp_to_ovlptmp( &ovlptmp, index, nold );
      /*###############################################*/

      /* Now, free all the old stored memory */
      /*###############################################*/
      free_old_MTdat();
      /*###############################################*/

      /* Add/Initiate new MTs; and allocate memory
	 for all vectors in struct MT */

      if (TotFluxConst)
          newgen=newMTs+nlost;
      else
          newgen=newMTs;

      MT.number=nold+newgen;
      n=MT.number;

      /*Allocate memory for the new round*/
      /*###############################################*/
      allocate_MTdata(n);
      /*###############################################*/

      /* First, fill in the data on "old" MTs from MTtmp array */
      /*###############################################*/
      fill_old_MTdata(nold, &MTtmp);
      /*###############################################*/

      /* generate segment data in accord with b.c. */
      /*###############################################*/
      for (iMT=1;iMT<=nold;iMT++)  mkseg(iMT);
      /*###############################################*/

      //// transfer the overlaps from the unsegmented
      /*#####################################################################*/
      fill_old_OVLPdata( &ovlptmp, nold );
      /*#####################################################################*/

      /* free the tmporary vectors */
      freetmp( &MTtmp, &ovlptmp, nold, index);
  }

  /* -----------------------------------------------------------------*/
  /* Generate new MTs */

  if (TotFluxConst)
    newgen=nlost+newMTs;
  else
    newgen=newMTs;

  MT.number=nold+newgen;


  /* Store the iteration number in which the MT was generated */
  for (i=1;i<=newgen;i++) {
      iMT=i+nold;
      MT.iter0[iMT]=iter;
  }

  /* Then, first randomly choose the Length of the newly generated MTs.
   And if lboundtype==NET make them not connected to the net, but able to
   become connected (set 0) explaination of the 3 values of MT.is_in_net in
   equations.c
   */
  for (i=1;i<=newgen;i++){
      iMT=nold+i;
      // Two length populations
      /*if (ran1(&idum)<0.5)
	MT.length[iMT]=MT.maxL;
      else
      MT.length[iMT]=MT.minL;*/
      // Linear length distribution
      //MT.length[iMT]=(MT.maxL-MT.minL)*ran1(&idum)+MT.minL;
      //exponential dist
      ranexp=ran1exp(&idum);
      while (ranexp>1) {
          ranexp=ran1exp(&idum);
          printf("rand must be corrected!\n");
      }
      if (ranexp>1) {
          printf("ERROR in expo random!\n");
          exit(0);
      }
      MT.length[iMT]=(MT.maxL-MT.minL)*ranexp+MT.minL;
      MT.is_in_net[iMT]=0;
      ovlp.actin[iMT]=0;
  }
/* Generate random positions in a Cylindrical domain depending
 on b.c and addmodel, and width 2*bundle.rad in the YZ plane.*/
  /*####################################################*/
  if (TotFluxConst)
    {
      istart=nold;
      getcoor(istart,nlost);

      istart=nold+nlost;
      getcoor(istart,newMTs);
    }
  else
    {
      istart=nold;
      getcoor(istart,newMTs);
    }


  /*####################################################*/
  /* Generate an unique index for the new MTs */
  /*####################################################*/


    i=1;
    iMT=nold+1;
    while (iMT<=MT.number) {

        if (i>mxmts*niter) {
            printf("More Mts than expected entered the system. Please adjust the allocation of MT.index_count in initialize_variables.c\n");
            exit(0);
        }

        if (track.index_count[i]==0) {
            MT.index[iMT]=i;
            iMT++;
            track.index_count[i]=1;
        }
        i++;
    }



  /* Set to zero the overlapMap of newly generated MTs */

  /*--------|---------------------------------   */
  /* iMT    |    1   |   2   |   3   |   4   |   */
  /*    ----|--------|-------|-------|-------|   */
  /*   |iseg|  1 | 2 | 1 | 2 | 1 | 2 | 1 | 2 |   */
  /*---|----|----------------|---|---|---|---|   */
  /*   |  1 |                |   |   |   |   |   */
  /* 1 |    |                |---|---|---|---|   */
  /*   |  2 |                |   |   |   |   |   */
  /*---|----|     OLD MTs    |---|---|---|---|   */
  /*   | 1  |                |   |   |   |   |   */
  /* 2 |    |                |---|---|---|---|   */
  /*   | 2  |                |   |   |   |   |   */
  /*---|----|---|---|---|----|---|---|---|---|---*/
  /*   |  1 |   |   |   |    | 0 | 0 |   |   |   */
  /* 3 |    |---|---|---|----|---|---|---|---|   */
  /*   |  2 |   |   |   |    | 0 | 0 |   |   |   */
  /*---|----|---|---|---|----|---|---|---|---| i */
  /*   |  1 |   |   |   |    |   |   | 0 | 0 |   */
  /* 4 |    |---|---|---|----|---|---|---|---|   */
  /*   |  2 |   |   |   |    |   |   | 0 | 0 |   */
  /*-----------------------------------------|---*/
  /*                   j                     |   */
  /*-----------------------------------------|---*/

  for (i=1;i<=newgen;i++)
  {
      iMT=i+nold;

      iseg1=2*iMT-1;
      iseg2=2*iMT;

      ovlp.type[iseg1][iseg1]=ZERO;
      ovlp.type[iseg2][iseg2]=ZERO;
      ovlp.type[iseg1][iseg2]=ZERO;
      ovlp.type[iseg2][iseg1]=ZERO;

      for(jMT=1;jMT<=MT.number;jMT++)
	{
	  jseg1=2*jMT-1;
	  jseg2=2*jMT;

	  ovlp.type[iseg1][jseg1]=ZERO;
	  ovlp.type[jseg1][iseg1]=ZERO;
	  ovlp.type[iseg2][jseg2]=ZERO;
	  ovlp.type[jseg2][iseg2]=ZERO;

	  ovlp.type[iseg1][jseg2]=ZERO;
	  ovlp.type[jseg2][iseg1]=ZERO;
	  ovlp.type[iseg2][jseg1]=ZERO;
	  ovlp.type[jseg1][iseg2]=ZERO;
	}
  }

  /*####################################################*/
  /* Remove MTs that failed to grow for GC add model*/
  /*####################################################*/
  nold=0;
  /* count the number of valid MTs; see boundary.c */
  for (iMT=1;iMT<=MT.number;iMT++)
    if (MT.length[iMT]>0) // maybe boundary sets length zero
      nold++;// only BC surviving MTs are counted

  index=ivector(1,nold);
  /*Create tmp array to save all MT data*/
  /*#############################*/
  allocate_tmp( &MTtmp , &ovlptmp);
  /*#############################*/

  /* Generate a pointer to the old indexing. Note that
       only MTs with length>0 "survive"*/
  i=0;
  for (iMT=1;iMT<=MT.number;iMT++){
    if (MT.length[iMT]>0){
      i++;
      index[i]=iMT;
    }
	  else{
		  printf("MT failed to grow: %i with cm-lb= %f\n", iMT,MT.cm[iMT].x-lbound0);
	  }
	}
  /* Store MT array in above defined MTtmp array
  STORING BECAUSE NEW MT ENLARGE THE ARRAY! C DOES NOT SUPPORT DYNAMIC ARRAYS!
  MTs that in their new position cross the
  boundary are eliminated, because nold only counts microtubule with l>0*/

  /*################################################*/
  MT_to_tmp( &MTtmp, index, nold );
  /*################################################*/

  /* Store data on Motor orientation in each overlaping region in temporary ovlp;
      Store overlapping map. Update the overlap map by the new indexing */
  /*###############################################*/
  ovlp_to_ovlptmp( &ovlptmp, index, nold );
  /*###############################################*/

  /* Now, free all the old stored memory */
  /*###############################################*/
  free_old_MTdat();
  /*###############################################*/

  MT.number=nold;
  n=MT.number;

  /*Allocate memory for the new round*/
  /*###############################################*/
  allocate_MTdata(n);
  /*###############################################*/

  /* First, fill in the data on "old" MTs from MTtmp array */
  /*###############################################*/
  fill_old_MTdata(nold, &MTtmp);
  /*###############################################*/

  /* generate segment data in accord with b.c. */
  /*###############################################*/
  for (iMT=1;iMT<=nold;iMT++)  mkseg(iMT);
  /*###############################################*/

  //// transfer the overlaps from the unsegmented
  /*#####################################################################*/
  fill_old_OVLPdata( &ovlptmp, nold );
  /*#####################################################################*/

  /* free the tmporary vectors */
  freetmp( &MTtmp, &ovlptmp, nold, index);


  /*####################################################*/
  /* Determine the index of the two most extreme MTs
   along x-axis; and their xcm */
  /*####################################################*/

  MT.imin=iMTcmmin();
  MT.imax=iMTcmmax();
  MT.iymin=iMTycmmin();
  MT.iymax=iMTycmmax();
  MT.izmin=iMTzcmmin();
  MT.izmax=iMTzcmmax();

  MT.xmin=MT.cm[MT.imin].x;
  MT.xmax=MT.cm[MT.imax].x;
  MT.ymin=MT.cm[MT.iymin].y;
  MT.ymax=MT.cm[MT.iymax].y;
  MT.zmin=MT.cm[MT.izmin].z;
  MT.zmax=MT.cm[MT.izmax].z;

  MT.rimin=riMTcmmin();
  MT.rimax=riMTcmmax();
  MT.rxmin=MT.rcm[MT.rimin].x;
  MT.rxmax=MT.rcm[MT.rimax].x;

  /*As above only takes the CM movement into account we also need
   other variables to determine the REAL xmax (including different lengths)
   NOTE THAT THIS IS ONLY THE CENTER OF MASS!!!!!!!!!!!!!!!!*/
  MT.imin_real=iMTxmin();
  MT.imax_real=iMTxmax();
  MT.xmin_real=MT.cm[MT.imin_real].x;

  MT.xmax_real=MT.cm[MT.imax_real].x;

  /*####################################################*/
  /* Finaly calculate the number of segments in the system */
  /*####################################################*/
  seg.number=0;
  for (iMT=1;iMT<MT.number;iMT++)
    {
      iseg1=2*iMT-1;
      iseg2=2*iMT;

      if (seg.length[iseg1]>0.0) seg.number++;
      if (seg.length[iseg2]>0.0) seg.number++;
    }





  if (!(iter%wfreq))
    {
      printf("----------------------------------------------\n");
      printf("finished iter %d with %d MTs\n",iter,MT.number);
      printf("Ensemble averaging iteration is: %d\n",iav);
      printf("Writing output to directory: %s\n\n",dirname);
      if (TotFluxConst)
	{
	  printf("Total Flux fixed: %d MTs every %f seconds\n",
		 newMTs0,newfreq*dt);
	  printf("Added %d MTs; Lost %d MTs\n",newgen,nlost);
	}
      else
	{
	 printf("Inward Flux fixed: %d MTs every %f seconds\n",
		 newMTs0,newfreq*dt);
	  printf("Added %d MTs; Lost %d MTs\n",newgen,nlost);

	}
    }


}

/*------------------------------------------------------------*/


void getcoor(int istart, int n)
{

  int i,j,iMT,jMT,k,counter,iseg1,iseg2,r_count,a_count, polarity;
  float xcm,dist2,tmp,lb=0,rb=0;
  float nr,factor,r,criteria,dtip;
  int np;
  double dphi,phi,radi,r1,r2,xov,output;

  for (i=1;i<=n;i++){
    iMT=istart+i;

    /*Get the boundaries for adding MTs. Depends on bc and addmodel*/
    if (istart!=0 || iter==1) {
      add_bounds(iMT,&lb,&rb);
    }
    else{
      rb=0;
      lb=0;
    }

    /* randomly pick x-coor in specified range depending on
    boundary condition and addmodel*/
    ///////////////////////////

    xcm=lb+((rb-lb)*ran1(&idum));
    printf("X-coordinate of new MT centre is: %f\n",xcm);

    ///////////////////////////
    /*Depreceated twopop if statement ... ignore*/
    if(!twopop){
	    MT.cm[iMT].x=xcm;
	    MT.rcm[iMT].x=0.0;
	    /*MT.rcm[iMT].x=xcm;*/
    }
    else{
	    if (iMT%2 != 0)
	      MT.cm[iMT].x=xcm;
	    else
	      MT.cm[iMT].x=-MT.cm[iMT-1].x;
	  }


    /* if first iteration select polarityratio0 as criteria*/
    if (iter==1)
	    criteria=PolarityRatio0;

    else{

      if (polaritymodel==FIXED){
        criteria=PolarityRatioNew;
      }
      /*Local polarity model calculates the local number of MTs and bases new polarity ratio on existing template MTs*/
      if (polaritymodel==LOCAL || polaritymodel==GCLOCAL){
        /*New 16/06/2020 allow for local orientation to affect probability of polarity ratio*/
        a_count=0;
        r_count=0;

	      printf("new MT pos: %f\n",MT.cm[iMT].x);
        for (jMT=1;jMT<=MT.number;jMT++){
          tmp=fabs(MT.cm[jMT].x-MT.cm[iMT].x);
          if(tmp<MT.length[jMT]/2.0){
            if(MT.direct[jMT].x<0.){
              a_count++;
            }
            else{
              r_count++;
            }
          }
        }
        //Total nuber of mts
        tmp=((float)(a_count+r_count));
	      printf("N local MTs: %f\n",tmp);
        //If no mts found just select default polarity ratio
        if (tmp<1.){
            tmp=PolarityRatioNew;
        }
        else{
          tmp=(float)(a_count)/tmp;
        }
        /*NEW END tmp now is the local polarity ratio*/
        criteria=tmp;
      }
      /*If no criteria chosen through local model just chose default*/
      else{
        criteria=PolarityRatioNew;
      }
      if(DEBUG)
        printf("Polarity criteria: %f MTs\n",criteria);
      /* Roll the dice to determine whether new MT is +end out or -end out*/
      r=ran1(&idum);

      ////////////////////////////////////////////////////
      /************** GC POLARITY MODEL *****************/
      ////////////////////////////////////////////////////
      if (polaritymodel==GC || polaritymodel==GCLOCAL){
        /*get distance from tip*/
        dtip = rb - xcm;
        /* Delete MT if growth is predicted to be unbounded*/
        if (r<criteria){

          if ( ran1(&idum) > 0.109647 + (0.664258*exp(-0.0540576*dtip)) ){
            if(DEBUG)
              printf("+end out MT failed to polymerize with threshold %f \n",0.109647 + (0.664258*exp(-0.0540576*dtip)));
            MT.length[iMT]=0.0;
          }
          else{
            if(DEBUG)
              printf("+end out MT polymerized with threshold %f \n",0.109647 + (0.664258*exp(-0.0540576*dtip)));
          }
        }
        else {
          if ( ran1(&idum) > 0.417541*exp(-0.030379*dtip) ){
            if(DEBUG)
              printf("-end out MT failed to polymerize with threshold %f \n",0.417541*exp(-0.030379*dtip));
            MT.length[iMT]=0.0;
          }
          else{
            if(DEBUG)
              printf("-end out MT polymerized with threshold %f \n",0.417541*exp(-0.030379*dtip));
          }
        }
      }
    }

    /* fill MT structure depending on polarity*/
    if (r<criteria){
	    if(!twopop){
	      MT.pend[iMT].x=MT.cm[iMT].x+MT.length[iMT]/2.0;
	      MT.mend[iMT].x=MT.cm[iMT].x-MT.length[iMT]/2.0;
	    }
	    else{
	      if (iMT%2 != 0) /* real MT */{
		      MT.pend[iMT].x=MT.cm[iMT].x+MT.length[iMT]/2.0;
		      MT.mend[iMT].x=MT.cm[iMT].x-MT.length[iMT]/2.0;
		    }
	      else /* image MT */{
		      MT.pend[iMT].x=-MT.pend[iMT-1].x;
		      MT.mend[iMT].x=-MT.mend[iMT-1].x;
		    }
	    }
	  }
    else{
	    if (!twopop){
        MT.pend[iMT].x=MT.cm[iMT].x-MT.length[iMT]/2.0;
        MT.mend[iMT].x=MT.cm[iMT].x+MT.length[iMT]/2.0;
           // MT.pend[iMT].x=MT.cm[iMT].x+MT.length[iMT]/2.0;
            //MT.mend[iMT].x=MT.cm[iMT].x-MT.length[iMT]/2.0;
	    }
	    else{
	      if (iMT%2 != 0){
          MT.pend[iMT].x=MT.cm[iMT].x-MT.length[iMT]/2.0;
          MT.mend[iMT].x=MT.cm[iMT].x+MT.length[iMT]/2.0;
		      //MT.pend[iMT].x=MT.cm[iMT].x+MT.length[iMT]/2.0;
		      //MT.mend[iMT].x=MT.cm[iMT].x-MT.length[iMT]/2.0;
		    }
	      else{
		      MT.pend[iMT].x=-MT.pend[iMT-1].x;
		      MT.mend[iMT].x=-MT.mend[iMT-1].x;
		    }
	    }
	  }
    /* Define the normalized direction of each MT */
    if (MT.length[iMT]>0.){
      MT.direct[iMT].x=(MT.mend[iMT].x-MT.pend[iMT].x)/MT.length[iMT];
    }
    else{
      MT.direct[iMT].x=0.;
    }
        //printf("%i %f %f\n", iMT, MT.direct[iMT].x, MT.length[iMT]);


      /* preform excluded volume interaction in z-y plane
	 namely, no two MTs can get closer than a cutoff
	 distance, '1.8*exclude'. */
      /* We randomely pick new yz-coor until the criteria
	 is met */


      if (!twopop || ((iMT%2) != 0))
	/* no two populations, or the real MTs in a twopop */
	{

      factor=1.0;
      counter=0;
	  do
	    {
	      counter++;
	      k=0;
	      /* choose coordinates randomely from an hexagonal latice */
	      /* pick radius */
	      radi=factor*bundle.rad*ran1(&idum);
	      phi=2*PI*ran1(&idum);
	      if (radi<exclude)
		{
		  MT.cm[iMT].z=0.0;
		  MT.rcm[iMT].z=0.0;
		  MT.cm[iMT].y=0.0;
		  MT.rcm[iMT].y=0.0;
		}
	      else
		{
		  nr=(float)floor(0.5+(radi+2*exclude)/(2*exclude));
		  radi=2*exclude*(nr-1);

		  /* pick angle */
		  if (nr>0) dphi=PI/(3.0*(nr-1));
		  np=(int)floor(0.5+phi/dphi);
		  phi=dphi*(double)floor(0.5+phi/dphi);
		  /* use cosine theorem to get the radius vector to the
		     relevant position on the hexagonal lattice */
		  r1=2*exclude*(np%(int)(nr-1));
		  r2=sqrtf(radi*radi+r1*r1-r1*radi);

            //printf("radius:%f\n",r2);

          MT.cm[iMT].z=r2*sin(phi);
		  MT.rcm[iMT].z=r2*sin(phi);
		  MT.cm[iMT].y=r2*cos(phi);
		  MT.rcm[iMT].y=r2*cos(phi);
       // printf ("ycoord: %f zcoord: %f\n", MT.cm[iMT].y,MT.cm[iMT].z);
		}

	      /* generate MT segments in accord with b.c. */
	      mkseg(iMT);

	      /*if (iter==1)
		{*/
		  for (jMT=1;jMT<iMT;jMT++)
		    {
		      xov=xoverlap(iMT,jMT);
		      dist2=(MT.cm[iMT].z-MT.cm[jMT].z)*
			(MT.cm[iMT].z-MT.cm[jMT].z)+
			(MT.cm[iMT].y-MT.cm[jMT].y)*
			(MT.cm[iMT].y-MT.cm[jMT].y);
		      if  (xov>0 && sqrtf(dist2)<1.8*exclude) k=1;
		    }
		  if (counter>MXTRIAL)
		    {
		      /* Increase the radius for adding MTs */
		      factor+=0.1;
		      /* printf("iter=%d iMT=%d factor=%f\n",iter,iMT,factor); */
		      /* 		      getchar(); */
		      /*printf("\a\a\n\nError in mkMT: Excluded volume criteria cound not be met!\ncounter>MXITER\ncounter=%d MXITER=%d exclude=%f\nNo space to add new MTs!\nIncrease MXITER, or decrease exlude\n",
			     counter,MXTRIAL,exclude);
			     exit(0);*/
		    }
		  /*}
		    else k=0;*/
	    }
	  while(k);

	  MT.pend[iMT].z=MT.cm[iMT].z;
	  MT.mend[iMT].z=MT.cm[iMT].z;
	  MT.direct[iMT].z=0.0;
	  MT.pend[iMT].y=MT.cm[iMT].y;
	  MT.mend[iMT].y=MT.cm[iMT].y;
	  MT.direct[iMT].y=0.0;
	}
      else if (twopop && iMT%2==0)
	/* yz-coor of image MTs in a two population ensemble */
	{
	  MT.cm[iMT].z=MT.cm[iMT-1].z;
	  MT.cm[iMT].y=MT.cm[iMT-1].y;
	  MT.mend[iMT].z=MT.mend[iMT-1].z;
	  MT.mend[iMT].y=MT.mend[iMT-1].y;
	  MT.pend[iMT].z=MT.pend[iMT-1].z;
	  MT.pend[iMT].y=MT.pend[iMT-1].y;
	  MT.direct[iMT].z=MT.direct[iMT-1].z;
	  MT.direct[iMT].y=MT.direct[iMT-1].y;
	  mkseg(iMT);
	}

    }

}

/*---------------------------------------------------------------*/


void allocate_tmp(struct MTdata *MTtmp, struct OVLPdata *ovlptmp){
    /* allocate memory for the index vector and the tmp vectors */

    MTtmp->length=vector(1,nold);
    MTtmp->cm=pvector(1,nold);   /* these are struct point vectors */
    MTtmp->rcm=pvector(1,nold); /* real position of MT */
    MTtmp->mend=pvector(1,nold);
    MTtmp->pend=pvector(1,nold);
    MTtmp->direct=pvector(1,nold);
    MTtmp->index=ivector(1,nold);
    MTtmp->iter0=ivector(1,nold);
    MTtmp->vel=pvector(1,nold);
    MTtmp->cluster=ivector(1,nold);
    MTtmp->is_in_net=ivector(1,nold);
    MTtmp->grow_direct=ivector(1,nold);

    /* ovlp is between MT segments, potentially, by Boundary condition,
	 each MT can include two segments. We keep track of
	 motor orientation in all overlap regions.
     */
    ovlptmp->type=(enum motortype **)imatrix(1,2*nold,1,2*nold);
    ovlptmp->iter0=matrix(1,2*nold,1,2*nold);
    ovlptmp->motor_direction=imatrix(1,2*nold,1,2*nold);
    ovlptmp->actin=ivector(1,nold);
    ovlptmp->actin_it0=ivector(1,nold);


}

void MT_to_tmp( struct MTdata *MTtmp, int *index, int nold ){
    int i,iMT;

    for ( i=1;i<=nold;i++)
	{
        iMT=index[i];

        MTtmp->index[i]=MT.index[iMT];
        MTtmp->iter0[i]=MT.iter0[iMT];
        MTtmp->length[i]=MT.length[iMT];
        MTtmp->cluster[i]=MT.cluster[iMT];
        MTtmp->is_in_net[i]=MT.is_in_net[iMT];
        MTtmp->grow_direct[i]=MT.grow_direct[iMT];

        MTtmp->cm[i].x=MT.cm[iMT].x;
        MTtmp->rcm[i].x=MT.rcm[iMT].x;
        MTtmp->mend[i].x=MT.mend[iMT].x;
        MTtmp->pend[i].x=MT.pend[iMT].x;
        MTtmp->direct[i].x=MT.direct[iMT].x;
        MTtmp->vel[i].x=MT.vel[iMT].x;

        MTtmp->cm[i].z=MT.cm[iMT].z;
        MTtmp->rcm[i].z=MT.rcm[iMT].z;
        MTtmp->mend[i].z=MT.mend[iMT].z;
        MTtmp->pend[i].z=MT.pend[iMT].z;
        MTtmp->direct[i].z=MT.direct[iMT].z;
        MTtmp->vel[i].z=MT.vel[iMT].z;

        MTtmp->cm[i].y=MT.cm[iMT].y;
        MTtmp->rcm[i].y=MT.rcm[iMT].y;
        MTtmp->mend[i].y=MT.mend[iMT].y;
        MTtmp->pend[i].y=MT.pend[iMT].y;
        MTtmp->direct[i].y=MT.direct[iMT].y;
        MTtmp->vel[i].y=MT.vel[iMT].y;
	}
}

void ovlp_to_ovlptmp( struct OVLPdata *ovlptmp, int *index, int nold ){
    int i, i1,i2,iseg1,iseg2,iMT, j,j1, j2, jseg1, jseg2, jMT;

    for (i=1;i<=nold-1;i++)
    {
        iMT=index[i];

        i1=2*i-1;
        i2=2*i;
        iseg1=2*iMT-1;
        iseg2=2*iMT;

        for (j=i;j<=nold;j++)
        {
            jMT=index[j];

            j1=2*j-1;
            j2=2*j;
            jseg1=2*jMT-1;
            jseg2=2*jMT;

            ovlptmp->type[i1][j1]=ovlp.type[iseg1][jseg1];
            ovlptmp->type[j1][i1]=ovlp.type[jseg1][iseg1];
            ovlptmp->type[i2][j2]=ovlp.type[iseg2][jseg2];
            ovlptmp->type[j2][i2]=ovlp.type[jseg2][iseg2];

            ovlptmp->type[i1][j2]=ovlp.type[iseg1][jseg2];
            ovlptmp->type[j2][i1]=ovlp.type[jseg2][iseg1];
            ovlptmp->type[i2][j1]=ovlp.type[iseg2][jseg1];
            ovlptmp->type[j1][i2]=ovlp.type[jseg1][iseg2];

            ovlptmp->iter0[i1][j1]=ovlp.iter0[iseg1][jseg1];
            ovlptmp->iter0[j1][i1]=ovlp.iter0[jseg1][iseg1];
            ovlptmp->iter0[i2][j2]=ovlp.iter0[iseg2][jseg2];
            ovlptmp->iter0[j2][i2]=ovlp.iter0[jseg2][iseg2];

            ovlptmp->iter0[i1][j2]=ovlp.iter0[iseg1][jseg2];
            ovlptmp->iter0[j2][i1]=ovlp.iter0[jseg2][iseg1];
            ovlptmp->iter0[i2][j1]=ovlp.iter0[iseg2][jseg1];
            ovlptmp->iter0[j1][i2]=ovlp.iter0[jseg1][iseg2];

            ovlptmp->motor_direction[i1][j1]=ovlp.motor_direction[iseg1][jseg1];
            ovlptmp->motor_direction[j1][i1]=ovlp.motor_direction[jseg1][iseg1];
            ovlptmp->motor_direction[i2][j2]=ovlp.motor_direction[iseg2][jseg2];
            ovlptmp->motor_direction[j2][i2]=ovlp.motor_direction[jseg2][iseg2];

            ovlptmp->motor_direction[i1][j2]=ovlp.motor_direction[iseg1][jseg2];
            ovlptmp->motor_direction[j2][i1]=ovlp.motor_direction[jseg2][iseg1];
            ovlptmp->motor_direction[i2][j1]=ovlp.motor_direction[iseg2][jseg1];
            ovlptmp->motor_direction[j1][i2]=ovlp.motor_direction[jseg1][iseg2];
        }
        /*Transfer actin ovlp*/
        ovlptmp->actin[i]=ovlp.actin[iMT];
        ovlptmp->actin_it0[i]=ovlp.actin_it0[iMT];
    }

    if (nold>0) {
        ovlptmp->actin[nold]=ovlp.actin[index[nold]];
        ovlptmp->actin_it0[nold]=ovlp.actin_it0[index[nold]];
    }

}


void free_old_MTdat(){


    free_vector(MT.length,1,MT.number);
    free_pvector(MT.cm,1,MT.number);
    free_pvector(MT.rcm,1,MT.number);
    free_pvector(MT.mend,1,MT.number);
    free_pvector(MT.pend,1,MT.number);
    free_pvector(MT.direct,1,MT.number);
    free_pvector(MT.vel,1,MT.number);
    free_ivector(MT.index,1,MT.number);
    free_ivector(MT.iter0,1,MT.number);
    free_ivector(MT.cluster,0,MT.number);

    free_ivector(MT.is_in_net,1,MT.number);
    free_ivector(MT.grow_direct,1,MT.number);
    free_imatrix((int **)ovlp.type,1,2*MT.number,1,2*MT.number);
    free_matrix(ovlp.iter0,1,2*MT.number,1,2*MT.number);
    free_imatrix(ovlp.motor_direction,1,2*MT.number,1,2*MT.number);
    free_ivector(ovlp.actin, 1,MT.number);
    free_ivector(ovlp.actin_it0, 1,MT.number);

    /* free all segment memory */
    free_vector(seg.length,1,2*MT.number);
    free_pvector(seg.cm,1,2*MT.number);
    free_pvector(seg.mend,1,2*MT.number);
    free_pvector(seg.pend,1,2*MT.number);
    free_pvector(seg.direct,1,2*MT.number);
}

void allocate_MTdata( int n){

    MT.length=vector(1,n);
    MT.cm=pvector(1,n);
    MT.rcm=pvector(1,n);
    MT.mend=pvector(1,n);
    MT.pend=pvector(1,n);

    MT.direct=pvector(1,n);
    MT.vel=pvector(1,n);
    MT.index=ivector(1,n);
    MT.iter0=ivector(1,n);

    MT.cluster=ivector(0,n);
    MT.is_in_net=ivector(1,n);
    MT.grow_direct=ivector(1,n);
    ovlp.type=(enum motortype **)imatrix(1,2*n,1,2*n);
    ovlp.iter0=matrix(1,2*n,1,2*n);
    ovlp.motor_direction=imatrix(1,2*n,1,2*n);
    ovlp.actin=ivector(1,n);
    ovlp.actin_it0=ivector(1,n);

    /* Allocate memory for segment data */
    seg.cm=pvector(1,2*n);
    seg.mend=pvector(1,2*n);
    seg.pend=pvector(1,2*n);
    seg.length=vector(1,2*n);
    seg.direct=pvector(1,2*n);
}

void fill_old_MTdata(int nold, struct MTdata *MTtmp){
    int iMT;
    for ( iMT=1;iMT<=nold;iMT++)
    {
        MT.length[iMT]=MTtmp->length[iMT];
        MT.index[iMT]=MTtmp->index[iMT];
        MT.iter0[iMT]=MTtmp->iter0[iMT];
        MT.cluster[iMT]=MTtmp->cluster[iMT];
        MT.is_in_net[iMT]=MTtmp->is_in_net[iMT];
        MT.grow_direct[iMT]=MTtmp->grow_direct[iMT];

        MT.cm[iMT].x=MTtmp->cm[iMT].x;
        MT.rcm[iMT].x=MTtmp->rcm[iMT].x;
        MT.mend[iMT].x=MTtmp->mend[iMT].x;
        MT.pend[iMT].x=MTtmp->pend[iMT].x;
        MT.direct[iMT].x=MTtmp->direct[iMT].x;
        MT.vel[iMT].x=MTtmp->vel[iMT].x;

        MT.cm[iMT].y=MTtmp->cm[iMT].y;
        MT.rcm[iMT].y=MTtmp->rcm[iMT].y;
        MT.mend[iMT].y=MTtmp->mend[iMT].y;
        MT.pend[iMT].y=MTtmp->pend[iMT].y;
        MT.direct[iMT].y=MTtmp->direct[iMT].y;
        MT.vel[iMT].y=MTtmp->vel[iMT].y;

        MT.cm[iMT].z=MTtmp->cm[iMT].z;
        MT.rcm[iMT].z=MTtmp->rcm[iMT].z;
        MT.mend[iMT].z=MTtmp->mend[iMT].z;
        MT.pend[iMT].z=MTtmp->pend[iMT].z;
        MT.direct[iMT].z=MTtmp->direct[iMT].z;
        MT.vel[iMT].z=MTtmp->vel[iMT].z;
    }
}

void fill_old_OVLPdata( struct OVLPdata *ovlptmp, int nold ){
    int iMT, jMT, iseg1,jseg1, iseg2, jseg2;

    for (iMT=1;iMT<=nold-1;iMT++){
        iseg1=2*iMT-1;
        iseg2=2*iMT;

        ovlp.type[iseg1][iseg1]=ovlptmp->type[iseg1][iseg1];
        ovlp.type[iseg1][iseg2]=ovlptmp->type[iseg1][iseg2];
        ovlp.type[iseg2][iseg1]=ovlptmp->type[iseg2][iseg1];
        ovlp.type[iseg2][iseg2]=ovlptmp->type[iseg2][iseg2];

        ovlp.iter0[iseg1][iseg1]=ovlptmp->iter0[iseg1][iseg1];
        ovlp.iter0[iseg1][iseg2]=ovlptmp->iter0[iseg1][iseg2];
        ovlp.iter0[iseg2][iseg1]=ovlptmp->iter0[iseg2][iseg1];
        ovlp.iter0[iseg2][iseg2]=ovlptmp->iter0[iseg2][iseg2];

        ovlp.motor_direction[iseg1][iseg1]=ovlptmp->motor_direction[iseg1][iseg1];
        ovlp.motor_direction[iseg1][iseg2]=ovlptmp->motor_direction[iseg1][iseg2];
        ovlp.motor_direction[iseg2][iseg1]=ovlptmp->motor_direction[iseg2][iseg1];
        ovlp.motor_direction[iseg2][iseg2]=ovlptmp->motor_direction[iseg2][iseg2];

        for (jMT=iMT+1;jMT<=nold;jMT++){
            jseg1=2*jMT-1;
            jseg2=2*jMT;

            ovlp.type[iseg1][jseg1]=ovlptmp->type[iseg1][jseg1];
            ovlp.type[jseg1][iseg1]=ovlptmp->type[jseg1][iseg1];
            ovlp.type[iseg2][jseg2]=ovlptmp->type[iseg2][jseg2];
            ovlp.type[jseg2][iseg2]=ovlptmp->type[jseg2][iseg2];

            ovlp.type[iseg1][jseg2]=ovlptmp->type[iseg1][jseg2];
            ovlp.type[jseg2][iseg1]=ovlptmp->type[jseg2][iseg1];
            ovlp.type[iseg2][jseg1]=ovlptmp->type[iseg2][jseg1];
            ovlp.type[jseg1][iseg2]=ovlptmp->type[jseg1][iseg2];

            ovlp.iter0[iseg1][jseg1]=ovlptmp->iter0[iseg1][jseg1];
            ovlp.iter0[jseg1][iseg1]=ovlptmp->iter0[jseg1][iseg1];
            ovlp.iter0[iseg2][jseg2]=ovlptmp->iter0[iseg2][jseg2];
            ovlp.iter0[jseg2][iseg2]=ovlptmp->iter0[jseg2][iseg2];

            ovlp.iter0[iseg1][jseg2]=ovlptmp->iter0[iseg1][jseg2];
            ovlp.iter0[jseg2][iseg1]=ovlptmp->iter0[jseg2][iseg1];
            ovlp.iter0[iseg2][jseg1]=ovlptmp->iter0[iseg2][jseg1];
            ovlp.iter0[jseg1][iseg2]=ovlptmp->iter0[jseg1][iseg2];

            ovlp.motor_direction[iseg1][jseg1]=ovlptmp->motor_direction[iseg1][jseg1];
            ovlp.motor_direction[jseg1][iseg1]=ovlptmp->motor_direction[jseg1][iseg1];
            ovlp.motor_direction[iseg2][jseg2]=ovlptmp->motor_direction[iseg2][jseg2];
            ovlp.motor_direction[jseg2][iseg2]=ovlptmp->motor_direction[jseg2][iseg2];

            ovlp.motor_direction[iseg1][jseg2]=ovlptmp->motor_direction[iseg1][jseg2];
            ovlp.motor_direction[jseg2][iseg1]=ovlptmp->motor_direction[jseg2][iseg1];
            ovlp.motor_direction[iseg2][jseg1]=ovlptmp->motor_direction[iseg2][jseg1];
            ovlp.motor_direction[jseg1][iseg2]=ovlptmp->motor_direction[jseg1][iseg2];
        }
        /*Transfer actin ovlp*/
        ovlp.actin[iMT]=ovlptmp->actin[iMT];
        ovlp.actin_it0[iMT]=ovlptmp->actin_it0[iMT];
    }
    ovlp.actin[nold]=ovlptmp->actin[nold];
    ovlp.actin_it0[nold]=ovlptmp->actin_it0[nold];

}


void freetmp( struct MTdata *MTtmp, struct OVLPdata *ovlptmp , int nold , int *index ){

    free_ivector(index,1,nold);
    free_pvector(MTtmp->cm,1,nold);
    free_pvector(MTtmp->rcm,1,nold);
    free_pvector(MTtmp->mend,1,nold);
    free_pvector(MTtmp->pend,1,nold);
    free_pvector(MTtmp->direct,1,nold);
    free_vector(MTtmp->length,1,nold);
    free_ivector(MTtmp->index,1,nold);
    free_ivector(MTtmp->iter0,1,nold);
    free_ivector(MTtmp->cluster,1,nold);
    free_ivector(MTtmp->is_in_net,1,nold);
    free_ivector(MTtmp->grow_direct,1,nold);
    free_pvector(MTtmp->vel,1,nold);
    free_imatrix((int **)ovlptmp->type,1,2*nold,1,2*nold);
    free_matrix(ovlptmp->iter0,1,2*nold,1,2*nold);
    free_imatrix(ovlptmp->motor_direction,1,2*nold,1,2*nold);
    free_ivector(ovlptmp->actin,1,nold);
    free_ivector(ovlptmp->actin_it0,1,nold);
}

/*THis function determines the left and right bounds for inserting new MT
 depending on b.c. and addmodel*/
void add_bounds(int iMT,float *lb, float *rb){
    float tmp;

    /* Determine the allowed cm range for for adding the new MTs */
    if (iter>1 && addmodel==ALL){
        /* check left boundary*/
        /*If NON just allow anything in the range of MTs*/
        if (bound.type[0]==NON)
            *lb=MT.xmin;
        /*In case of any other bounds add in a way which d to reach into the bound as it does not matter*/
        else if (bound.type[0]==ABSORB || bound.type[0]==PERIODIC ||
                 bound.type[0]==POPUP || bound.type[0]==STICKY ||
                 bound.type[0]==NET )
            *lb=lbound0;
        /*Obviously for Wall the MTs are not allowed to reach into the bound*/
        else if (bound.type[0]==WALL)
            *lb=lbound0+MT.length[iMT]/2;
        else if (bound.type[0]==FORCE)
            *lb=MT.xmin_real-0.5*MT.length[MT.imin_real]+MT.length[iMT]/2;
        else
        {
            printf("mkMTs:> Left boundary type not properly specified!");
            exit(0);
        }
        /*Check right boundary*/
        /*If NON just allow anything in the range of MTs*/
        if (bound.type[1]==NON  || bound.type[1]==REFLECT)
            *rb=MT.xmax;
        else if (bound.type[1]==PERIODIC || bound.type[1]==ABSORB
                 || bound.type[1]==POPUP )
            *rb=rbound0;
        else if (bound.type[1]==RFIXREF)
            *rb=reflectb0;
        else if(bound.type[1]==FORCE)
            *rb=MT.xmax_real+0.5*MT.length[MT.imax_real]-MT.length[iMT]/2;

        else
        {
            printf("mkMTs:> Right boundary type not properly specified!");
            exit(0);
        }
    }
    /* add MTs to right end */
    else if (iter > 1 && addmodel==RB){
        /* check left boundary*/
        /*If NON just allow anything in the range of MT length*/
        if (bound.type[0]==ABSORB || bound.type[0]==PERIODIC ||
                 bound.type[0]==POPUP || bound.type[0]==STICKY ||
                 bound.type[0]==NET || bound.type[0]==WALL ||
                 bound.type[0]==NON || bound.type[0]==FORCE )
            *lb=MT.xmax_real+0.5*MT.length[MT.imax_real]-MT.length[iMT]/2;
        else
        {
            printf("mkMTs:> Left boundary type not properly specified!");
            exit(0);
        }
        /*Check right boundary*/
        /*If NON just allow anything in the range of MTs*/
        if (bound.type[1]==NON  || bound.type[1]==REFLECT || bound.type[1]==FORCE)
            *rb=MT.xmax_real+0.5*MT.length[MT.imax_real]-MT.length[iMT]/2;
        else if (bound.type[1]==PERIODIC || bound.type[1]==ABSORB
                 || bound.type[1]==POPUP )
            *rb=rbound0;
        else if (bound.type[1]==RFIXREF)
            *rb=reflectb0;
        else
        {
            printf("mkMTs:> Right boundary type not properly specified!");
            exit(0);
        }
    }
    /* add MTs around zero. Only useful for FORCE-FORCE boundary*/
    else if (iter > 1 && addmodel==CENTER){
        if (bound.type[0]==FORCE && bound.type[1]==FORCE) {
            *lb=-1;
            *rb=1;
        }
        else{
            printf("Invalid use of ZERO bound!\n");
            exit(0);
        }
    }

    /* add MTs to left end
     smoothly or POPUP the fixing below is
     wrong for MTs of different lengths...getcoor has a fix implemented*/
    else if (iter > 1 && addmodel==LB ){
        /* check left boundary*/
        /*If NON just take the farthest left MT as reference*/
        if (bound.type[0]==NON || bound.type[0]==FORCE)
            *lb=MT.xmin;
        /*In case of any other bounds position the MTs so that they
         are close to the bound compared to their length*/
        else if (bound.type[0]==ABSORB || bound.type[0]==PERIODIC ||
                 bound.type[0]==POPUP  || bound.type[0]==NET ||
                 bound.type[0]==CONFLX )
            *lb=lbound0;
        else if (bound.type[0]==WALL || bound.type[0]==STICKY){
            printf("You are trying to use a smooth addmodel LB with a WALL!\n See the problem?");
            exit(0);
        }
        else
        {
            printf("mkMTs:> Left boundary type not properly specified!");
            exit(0);
        }
        /*Check right boundary*/

        if (bound.type[1]==NON  || bound.type[1]==REFLECT)
            *rb=MT.xmin+MT.length[iMT]/2;
        else if (bound.type[1]==PERIODIC || bound.type[1]==ABSORB
                 || bound.type[1]==POPUP || bound.type[1]==RFIXREF
                 || bound.type[1]==FORCE)
            *rb=lbound0+MT.length[iMT]/2;
        else
        {
            printf("mkMTs:> Right boundary type not properly specified!");
            exit(0);
        }
    }
    /*LBPOP leads to MTs appearing in the system with full length*/
    else if (iter > 1 && addmodel==LBPOP ){
        /* check left boundary*/
        /*If NON just take the farthest left MT as reference*/
        if (bound.type[0]==NON)
            *lb=MT.xmin;
        /*In case of any other bounds position the MTs so that they all touch the left bound with their left end*/
        else if (bound.type[0]==ABSORB || bound.type[0]==PERIODIC ||
                 bound.type[0]==POPUP  || bound.type[0]==NET ||
                 bound.type[0]==CONFLX || bound.type[0]==WALL ||
                 bound.type[0]==STICKY )
            *lb=lbound0+MT.length[iMT]/2;
        else if (bound.type[0]==FORCE)
            *lb=MT.xmin_real-0.5*MT.length[MT.imin_real]+MT.length[iMT]/2;
        else
        {
            printf("mkMTs:> Left boundary type not properly specified!");
            exit(0);
        }
        /*Check right boundary*/

        if (bound.type[1]==NON  || bound.type[1]==REFLECT)
            *rb=MT.xmin+MT.length[iMT]/2;
        else if (bound.type[1]==PERIODIC || bound.type[1]==ABSORB
                 || bound.type[1]==POPUP || bound.type[1]==RFIXREF )
            *rb=lbound0+MT.length[iMT]/2;
        else if (bound.type[1]==FORCE)
            *rb=*lb;
        else
        {
            printf("mkMTs:> Right boundary type not properly specified!");
            exit(0);
        }

    }
    else if (iter >1 && addmodel==RLB)
    {
        /* add MTs to both ends;
         to each end one every two iterations */
        if ((iter%2) == 0){
            /* check left boundary*/
            /*If NON just take the farthest left MT as reference*/
            if (bound.type[0]==NON)
                *lb=MT.xmin;
            /*In case of any other bounds position the MTs so that they all touch the left bound with their left end*/
            else if (bound.type[0]==ABSORB || bound.type[0]==PERIODIC ||
                     bound.type[0]==POPUP  || bound.type[0]==NET ||
                     bound.type[0]==CONFLX || bound.type[0]==WALL ||
                     bound.type[0]==STICKY )
                *lb=lbound0+MT.length[iMT]/2;
            else if (bound.type[0]==FORCE)
                *lb=MT.xmin_real-0.5*MT.length[MT.imin_real]+MT.length[iMT]/2;
            else
            {
                printf("mkMTs:> Left boundary type not properly specified!");
                exit(0);
            }
            /*Check right boundary*/

            if (bound.type[1]==NON  || bound.type[1]==REFLECT)
                *rb=MT.xmin+MT.length[iMT]/2;
            else if (bound.type[1]==PERIODIC || bound.type[1]==ABSORB
                     || bound.type[1]==POPUP || bound.type[1]==RFIXREF )
                *rb=lbound0+MT.length[iMT]/2;
            else if (bound.type[1]==FORCE)
                *rb=MT.xmin_real-0.5*MT.length[MT.imin_real]+MT.length[iMT]/2;
            else
            {
                printf("mkMTs:> Right boundary type not properly specified!");
                exit(0);
            }

        }
        else{
            /* check left boundary*/
            /*If NON just allow anything in the range of MT length*/
            if (bound.type[0]==NON)
                *lb=MT.xmax-MT.length[iMT]/2;
            /*In case of any other bounds position the MTs so that they
             are close to the bound compared to their length*/
            else if (bound.type[0]==ABSORB || bound.type[0]==PERIODIC ||
                     bound.type[0]==POPUP || bound.type[0]==STICKY ||
                     bound.type[0]==NET || bound.type[0]==WALL )
                *lb=MT.xmax_real+0.5*MT.length[MT.imax_real]-MT.length[iMT]/2;
            else
            {
                printf("mkMTs:> Left boundary type not properly specified!");
                exit(0);
            }
            /*Check right boundary*/
            /*If NON just allow anything in the range of MTs*/
            if (bound.type[1]==NON  || bound.type[1]==REFLECT)
                *rb=MT.xmax;
            else if (bound.type[1]==PERIODIC || bound.type[1]==ABSORB
                     || bound.type[1]==POPUP )
                *rb=rbound0;
            else if (bound.type[1]==RFIXREF)
                *rb=reflectb0;
            else if(bound.type[1]==FORCE)
                *rb=MT.xmax_real+0.5*MT.length[MT.imax_real]-MT.length[iMT]/2;
            else
            {
                printf("mkMTs:> Right boundary type not properly specified!");
                exit(0);
            }

        }
    }
    else
    /* Otherwise, if iter==1 or if no model is selected
     lbound and rbound are given as input to first
     pull of MTs but plus their length so that nobody touches
     any bounds in the beginning
     */
    {
        *lb=lbound0;
        *rb=rbound0;

    }

    if (*rb<*lb)
    {
        printf("mkMTs:> iter=%i\n", iter);
        printf("mkMTs:> Error. rbound>lbound\a\n");
        printf("mkMTs:> AddModel selected is: %d\n",addmodel);
        printf("mkMTs:> lbound=%f rbound=%f\n",*lb,*rb);
        printf("mkMTs:> MT.cm=%f MT.length=%f\n",MT.cm[iMT].x, MT.length[iMT]);
        printf("mkMTs:> nMT=%d\n",MT.number);

        *rb=0.;
        *lb=0.;
        //exit(0);
    }


}
