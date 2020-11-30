/* #################################################################### */
/* update_av  --  Adds items to ensemble averaged quantaties            */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global_var.h"
#include "nrutil.h"
#define REF (bound.type[1]==REFLECT || bound.type[1]==RFIXREF)
#define PERIOD (bound.type[0]==PERIODIC && bound.type[1]==PERIODIC)
#define ABSORBING (bound.type[0]==ABSORB && bound.type[1]==ABSORB)
#define NEIGHBOR_CRITERIA ((neighbors[iMT]>0 && AvNeighboredMTs==TRUE) || AvNeighboredMTs==FALSE)
#define CANCEL TRUE /* allows to locally cancel the effect of 
			NEIGHBORING_CRITERIA; 
			if CANCEL is set to FALSE, it eliminates the 
			cancelling effect.*/  
void update_av()
{
  /* This function updates averaged quantities such as the density 
     profiles: av.nl=<nl(t,x)>, av.nr=<nr(t,x)> and av.ntot=<ntot(t,x)>; 
     where nr and nl are the number of MTs pointing with their 
     minus-end to the right and left. 
     We also calculate av.xmin and av.xmax which are the average 
     right and left extreems of the growing domain.
  */

  int iMT,jMT,iseg,jseg,iseg1,iseg2,jseg1,jseg2,imaxvel,iminvel,ibox,z,iop;
  int zeros,*neighbors; /* zeros - the number of unbounded MTs */
                        /* neighbort[iMT] - the number of neighbors 
                           for each iMT */
  float xrange,zrange,dx,x,dv,rdx;
  float nr,nl,nt,xr,xl,xr2,xl2,xt,xt2,rxr,rxl,rxr2,rxl2,rxt,rxt2;
  float velt,velr,vell,velt2,velr2,vell2;
  float sz,ov,op;
  float xi,xj,vi,vj,s;
  float r,rmin[MXBOX],rmax[MXBOX],boxnr[MXBOX],boxnl[MXBOX],boxn[MXBOX];
  float box_yz[MXBOX2],polarcor_ryz[MXBOX2];
  float vel_nr[MXBOX],vel_nl[MXBOX],ntot,polarcor[MXBOX];
  float vpol; /* MT polimerization velocity */
  float velmin,velmax; /* for histograms, minimum and maximum velocities */
  
  int ncheck; /* just for debugging, 
		 checks that the sum in each iteration gives MT.number */

  int tau,it;
  float *ntau,*nltau,*nrtau,x2;

  /* For use in several calculations below, find which MTs 
     have no neighbors.
  */

  neighbors=ivector(1,MT.number);

  zeros=0; 
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      neighbors[iMT]=0;
      for (jMT=1;jMT<=MT.number;jMT++)
	if (jMT!=iMT && overlap(iMT,jMT)>0) neighbors[iMT]++;
      if (neighbors[iMT]==0) zeros++;
    }  
  

  /* calculate how many MTs with each polarity */
  av.nMT[iter]+=MT.number;
  av.nzeros[iter]+=zeros;
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      if(MT.direct[iMT].x>0)
	  av.nMTr[iter]++;
      else
	av.nMTl[iter]++;
    }

  /* Calculate the average separation between most distant MTs 
     in the shaft; for each time step. We need 2 cases for this:
   The most extreme center of mass and the most extreme end of a microtubule
   because MT lengths may vary so that extremecm != extremeMT*/
  /* Do the same for what we caled real positions of MTs, 
     i.e., positions not subjected to implications of b.c. */
  
  dx=0.00000001+(MT.xmax-MT.xmin)/(nbox-1);
  rdx=0.0000001+(MT.rxmax-MT.rxmin)/(nbox-1);
  
  av.xmin[iter] += MT.xmin;
  av.xmax[iter] += MT.xmax;

  av.xmax_real[iter] += MT.xmax_real+MT.length[MT.imax_real]/2;
  av.xmin_real[iter] += MT.xmin_real-MT.length[MT.imin_real]/2;
    
  av.zmin[iter] += MT.zmin;
  av.zmax[iter] += MT.zmax;

  av.rxmin[iter] += MT.rxmin;
  av.rxmax[iter] += MT.rxmax; 


  /* Calculate an averaged density profile for the current timestep (iter) */
  /* We also calculate at the same time, an order parameter profile, and 
     average velocity profile.                                             */ 
  iop=0;
  op=0.0;
  ncheck=0;
  for (ibox=0;ibox<nbox;ibox++)
    {boxnr[ibox]=0.0; boxnl[ibox]=0.0;}
  
  for (iMT=1;iMT<=MT.number;iMT++){
      ibox=(int)rint((MT.cm[iMT].x-MT.cm[MT.imin].x)/dx);
      /*printf("dx=%f x=%f x/dx=%f ibox=%d\n",dx,MT.cm[iMT].x-MT.cm[MT.imin].x,
	(MT.cm[iMT].x-MT.cm[MT.imin].x)/dx,ibox);getchar();*/
      if (ibox<0 || ibox>=nbox) /* nbox is a global variable */ {
	  printf("Update_av:> nMT=%d  ibox=%d MT.cm[%d].x=%f dx=%f\n",
		 MT.number,ibox,MT.index[iMT],MT.cm[iMT].x,dx);
	  printf("Update_av:> xmax=%f xmin=%f x=xcm-xmin=%f dx=%f x/dx=%f ibox=%d\n",
		 MT.xmax,MT.xmin,MT.cm[iMT].x-MT.cm[MT.imin].x,dx,
		 (MT.cm[iMT].x-MT.cm[MT.imin].x)/dx,ibox);
	  getchar();
	}
   
      
      if (iMT != MT.imin && iMT != MT.imax && (CANCEL || NEIGHBOR_CRITERIA)) {
	  /* The two extreem MTs, are excluded from the histogram calc, 
	     since they are used to calculate the histogram boundaries, 
	     and therefore by their definition they will have a finite 
	     probablity in each configuration.
	     If "|| CANCEL" is used above, it locally cancels the 
	     NEIGHBORING_CRITERIA condition.
	  */
	 
          if (MT.direct[iMT].x>0){
              ncheck++;
              av.nr[iter][ibox]++;
              av.ntot[iter]++;
              av.xvelr[iter][ibox]+=MT.vel[iMT].x;
	      }
          else{
              ncheck++;
              av.nl[iter][ibox]++;
              av.ntot[iter]++;
              av.xvell[iter][ibox]+=MT.vel[iMT].x;
	      }
	  }

      if (MT.direct[iMT].x>0)
          boxnr[ibox]++;
      else
          boxnl[ibox]++;
  }
  
  /*OrderParameter Profile (av.op) at timestep iter*/
  /*op=0;
  iop=0;
  for(ibox=0; ibox<nbox; ibox++){
      if ((boxnr[ibox]+boxnl[ibox])>0){*/
	  /* if ((av.nr[iter][ibox]+av.nl[iter][ibox])>0) */
	  /* 	{ */
	  /* s=(av.nr[iter][ibox]-av.nl[iter][ibox])/ */
	  /* 	    (av.nr[iter][ibox]+av.nl[iter][ibox]); */
         /* s=(boxnr[ibox]-boxnl[ibox])/(boxnr[ibox]+boxnl[ibox]);
          op += s*s;
          iop++;
      }
  }
  av.op[iter]+=(op/(float)iop);*/

  /* ----------------------------------------------------------*/
/*Calculate average tubulin density (also plus end out tubulin and minus end out tubulin) profile*/
    float left_domain_end=MT.xmin_real-0.5*MT.length[MT.imin_real]; //the left end of the domain
    float boxlength= (MT.xmax_real+0.5*MT.length[MT.imax_real]-MT.xmin_real+0.5*MT.length[MT.imin_real])/(float)nbox;
    for (ibox=1; ibox<=nbox; ibox++) {
        /*#################################################################################*/
        for (iseg=1;iseg<=2*MT.number;iseg++) {
            
            // only carry on in this loop if the segment is existing and overlapping with ibox
            if (seg.length[iseg]<0.0001) continue;
            ov = Interval_overlap( (ibox-1)*boxlength+left_domain_end ,seg.cm[iseg].x-seg.length[iseg]/2. , boxlength , seg.length[iseg] );
            if (ov<0.0001) continue;
            
            iMT=segindex(iseg);
            if (iMT != MT.imin_real && iMT != MT.imax_real && (CANCEL || NEIGHBOR_CRITERIA)){
                av.tubulin_dens[iter][ibox]+=ov/t_length;
                if (seg.direct[iseg].x>0)
                    av.tubulin_r[iter][ibox]+=ov/t_length;
                else if (seg.direct[iseg].x<0)
                    av.tubulin_l[iter][ibox]+=ov/t_length;
            }
            /*Get directed tubulin minus end to the right*/
        }
    }

    
    /* ----------------------------------------------------------*/
    /*Get av active,kinesin,dynein and bipolar probability*/
    for (ibox=1; ibox<=nbox; ibox++) {
        av.ovlp_kin[iter][ibox]+= ProbKi_dist[ibox];
        av.ovlp_dyn[iter][ibox]+= ProbDy_dist[ibox];
        av.ovlp_bi[iter][ibox]+= ProbBi_dist[ibox];
        av.ovlp_act[iter][ibox]+= ProbAct_dist[ibox];
    }
    
    
    /* ----------------------------------------------------------*/
    
    
     /* ----------------------------------------------------------*/

 /* Calculate an averaged density profile for the real particle positions 
    current timestep (iter) */
  if (iter>neqsteps+5) /* the +5 is because all "real" coordinates 
			  start as zero, see newpos.c */
    {
  ncheck=0;
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      ibox=(int)rint((MT.rcm[iMT].x-MT.rcm[MT.rimin].x)/rdx);
      /*printf("dx=%f x=%f x/dx=%f ibox=%d\n",dx,MT.cm[iMT].x-MT.cm[MT.imin].x,
	(MT.cm[iMT].x-MT.cm[MT.imin].x)/dx,ibox);getchar();*/
      if (ibox<0 || ibox>=nbox) /* nbox is a global variable */
	{
	  printf("Update_av:> nMT=%d  ibox=%d MT.cm[%d].x=%f dx=%f\n",
		 MT.number,ibox,MT.index[iMT],MT.rcm[iMT].x,dx);
	  printf("Update_av:> xmax=%f xmin=%f x=xcm-xmin=%f dx=%f x/dx=%f ibox=%d\n",
		 MT.rxmax,MT.rxmin,MT.rcm[iMT].x-MT.rcm[MT.rimin].x,dx,
		 (MT.rcm[iMT].x-MT.rcm[MT.imin].x)/rdx,ibox);
	  getchar();
	}
   
      
      if (iMT != MT.rimin && iMT != MT.rimax && (CANCEL || NEIGHBOR_CRITERIA))
	  /* The two extreem MTs, are excluded from the histogram calc, 
	     since they are used to calculate the histogram boundaries, 
	     and therefore by their definition they will have a finite 
	     probablity in each configuration. 
	     If "|| CANCEL" is used above, it locally cancels the 
	     NEIGHBORING_CRITERIA condition.
	  */
	  {
	    if (MT.direct[iMT].x>0)
	      {
		av.rnr[iter][ibox]++;
	      }
	    else
	      {
		av.rnl[iter][ibox]++;
	      }
	  }




    }
    }

  /* if ((AvNeighboredMTs && ncheck != MT.number-2-zeros) ||  */
/*       (!AvNeighboredMTs && ncheck != MT.number-2))  */
    
/*     /\* -2 is due to the extreem MTs that we exclude *\/ */
    
/*     { */
/*       printf("Update_av:> Error in wdensity\nncheck=%d\nMT.number=%d\nMTs with zero neighbors=%d",ncheck,MT.number,zeros); */
/*       getchar(); */
/*     } */

  /*-----------------------------------------------------------------------*/
  /* Calculate Histograms of Absolute Velocity                             */

  for(ibox=0;ibox<nbox;ibox++){
      vel_nr[ibox]=0.0;
      vel_nl[ibox]=0.0;

  }
  ntot=0.0;
  
  iminvel=iabsVelmin();
  imaxvel=iabsVelmax();

  /* below, the term: -velpol*MT.direct[iminvel].x, 
     is the velocity associated with the PLUS-end growth of the MT. 
     The minus sign is because MT.direct is towards the minus end.  
  */
  /*get step size of velocity sum*/
  dv=(fabs(MT.vel[imaxvel].x/*-velpol*MT.direct[imaxvel].x*/)
      -fabs(MT.vel[iminvel].x/*-velpol*MT.direct[iminvel].x*/))/(nbox-1);

  if (dv!=0.0){
      av.absvmin[iter] += fabs(MT.vel[iminvel].x/*-velpol*MT.direct[iminvel].x*/);
      av.absvmax[iter] += fabs(MT.vel[imaxvel].x/*-velpol*MT.direct[imaxvel].x*/);
      
      for (iMT=1;iMT<=MT.number;iMT++){
          ibox=(int)rint((fabs(MT.vel[iMT].x/*-velpol*MT.direct[iMT].x*/)
			  -fabs(MT.vel[iminvel].x/*-velpol*MT.direct[iminvel].x*/))/dv);
	  
          if (ibox<0 || ibox>=nbox) /* nbox is a global variable */{
              printf("Update_av:> Error absVelocity out of range\n");
              printf("Update_av:> iminvel=%d imaxvel=%d\n",iminvel,imaxvel);
              printf("Update_av:> MT.vel[%d].x=%f MT.vel[%d].x=%f dv=%f\n",
		            iminvel,MT.vel[iminvel].x,imaxvel,MT.vel[imaxvel].x,dv);
              printf("Update_av:> ibox=%d MT.vel[%d].x=%f dv=%f\n",
		            ibox,iMT,MT.vel[iMT].x,dv);
	      
	          for (jMT=1;jMT<=MT.number;jMT++)
                  printf("vel[%d]=%f\n",jMT,MT.vel[jMT].x);
	      
	      
              getchar();
	      }

          if (iMT != iminvel && iMT != imaxvel && NEIGHBOR_CRITERIA){
	    /* The two extreem MTs, are excluded from the histogram calc, 
	       since they are used to calculate the histogram boundaries, 
	       and therefore by their definition they will have a finite 
	       probablity in each configuration. */
	    
              if (fabs(MT.vel[iMT].x/*-velpol*MT.direct[iMT].x*/)>0){
                  if (MT.direct[iMT].x>0){
                      ntot++;
                      
                      vel_nr[ibox]++;
                  }
                  else{
                      ntot++;
                      vel_nl[ibox]++;
                  }
              }
          }
      }
  }

  /* normalize average fraction of MTs moving with vel ibox*dv*/
  for (ibox=0;ibox<nbox;ibox++){
      if (ntot>0){

          av.absvel_nl[iter][ibox]+=(vel_nl[ibox]/ntot);
          av.absvel_nr[iter][ibox]+=(vel_nr[ibox]/ntot);
      }
  }
    
 /*-----------------------------------------------------------------------*/

  /*Calculate av cluster velocities*/
    int clust_tot[MT.number+1];
    float vel_clust[MT.number+1][nbox+1];
    for (ibox=0;ibox<=nbox;ibox++){
        for (iMT=0; iMT<=MT.number; iMT++) {
            vel_clust[iMT][ibox]=0.;
            clust_tot[iMT]=0;
        }
    }

    dv=(fabs(MT.vel[imaxvel].x/*-velpol*MT.direct[imaxvel].x*/)
        -fabs(MT.vel[iminvel].x/*-velpol*MT.direct[iminvel].x*/))/nbox;
    if (dv!=0.0) {
        
        for (iMT=1; iMT<=MT.number; iMT++) {
            

            ibox=(int)floor((fabs(MT.vel[iMT].x/*-velpol*MT.direct[iMT].x*/)
                            -fabs(MT.vel[iminvel].x/*-velpol*MT.direct[iminvel].x*/))/dv)+1;
            if (ibox==0)
                ibox=1;
            else if (ibox==nbox)
                ibox=nbox-1;
            if (iMT != iminvel && iMT != imaxvel && NEIGHBOR_CRITERIA){
           
                vel_clust[MT.cluster[iMT]][ibox]++;
                clust_tot[MT.cluster[iMT]]++;

            }
        }
    }
    for (ibox=1;ibox<=nbox;ibox++){
        for (iMT=1; iMT<=MT.number; iMT++) {
            av.nv_cluster[MT.cluster[iMT]][ibox]+=(vel_clust[MT.cluster[iMT]][ibox]/(float)clust_tot[ibox]);
        }
    }
    


  /*-----------------------------------------------------------------------*/
  /* Calculate Histograms of actual Velocity                             */

  for(ibox=0;ibox<nbox;ibox++){
      vel_nr[ibox]=0.0;
      vel_nl[ibox]=0.0;
  }
  ntot=0.0;  
  
  iminvel=iVelmin();
  imaxvel=iVelmax();
  
  
  velmin=MT.vel[iminvel].x;
  velmax=MT.vel[imaxvel].x;
  
  /* Uncomment to use fixed values for the min/max values of Vel */ 
  /* these are defined in rdata.c                                */
  /* also, update wav_veldensity() in write_output.c             */ 

  /*velmin=hist_vel_min;
    velmax=hist_vel_max;
  */

  dv=(velmax-velmin)/(nbox-1);
   

  if (dv!=0.0)
    {
      av.vmin[iter] += velmin;
      av.vmax[iter] += velmax;
      
      ncheck=0;
      for (iMT=1;iMT<=MT.number;iMT++)
	{
	  if (MT.vel[iMT].x>velmin && MT.vel[iMT].x<velmax)
	    {
	      ibox=(int)rint((MT.vel[iMT].x-velmin)/dv);      
	      if (ibox<0 || ibox>=nbox) /* nbox is a global variable */
		{
		  printf("Update_av:> Error velocity out of range\n");
		  printf("Update_av:> iminvel=%d imaxvel=%d\n",iminvel,imaxvel);
		  printf("Update_av:> MT.vel[%d].x=%f MT.vel[%d].x=%f dv=%f\n",
			 iminvel,velmin,imaxvel,velmax,dv);
		  printf("Update_av:> ibox=%d MT.vel[%d].x=%f dv=%f\n",
			 ibox,iMT,MT.vel[iMT].x,dv);
		  
		  for (jMT=1;jMT<=MT.number;jMT++) 
		    printf("vel[%d]=%f\n",jMT,MT.vel[jMT].x);
		  
		  
		  getchar();
		}
	      
	      if (iMT != iminvel && iMT != imaxvel && NEIGHBOR_CRITERIA) 
	      /* The two extreem MTs, are excluded from the histogram calc, 
		 since they are used to calculate the histogram boundaries, 
		 and therefore by their definition they will have a finite 
		 probablity in each configuration. */
	      
	      /*if (NEIGHBOR_CRITERIA)*/
		{
		  if (MT.direct[iMT].x>0)
		    {
		      ntot++;
		      vel_nr[ibox]++;
		    }
		  else
		    {
		      ntot++;
		      vel_nl[ibox]++;
		    }
		  
		}
	      
	    }
	}
    }

 /* normalize */
  for (ibox=0;ibox<nbox;ibox++)
    {
      if (ntot>0)
	{
	  av.vel_nl[iter][ibox]+=(vel_nl[ibox]/ntot);
	  av.vel_nr[iter][ibox]+=(vel_nr[ibox]/ntot);
	}  
    }


  /*---------------------------------------------------------------*/
  /* more average quantaties should go below */


  /* <xr(t)>,<xr^2(t)>,<xl(t)>,<xl^2(t)> <xt(t)>,<xt^2(t)> */

  nr=0.0;
  nl=0.0;

  xr=0.0;
  xr2=0.0;
  xl=0.0;
  xl2=0.0;
  xt=0.0;
  xt2=0.0;
  
  rxr=0.0;
  rxr2=0.0;
  rxl=0.0;
  rxl2=0.0;
  rxt=0.0;
  rxt2=0.0;

  velt=0.0;
  velr=0.0;
  vell=0.0;

  velt2=0.0;
  velr2=0.0;
  vell2=0.0;

  for(iMT=1;iMT<=MT.number;iMT++)
    {      
      if (NEIGHBOR_CRITERIA || CANCEL) 
	/* If "|| CANCEL" is used above, it locally cancels the 
	   NEIGHBORING_CRITERIA condition.*/
	{
	  xt += MT.cm[iMT].x;
	  xt2 += MT.cm[iMT].x*MT.cm[iMT].x;
	  
	  rxt += MT.rcm[iMT].x;
	  rxt2 += MT.rcm[iMT].x*MT.rcm[iMT].x;
	  
	  velt += MT.vel[iMT].x;
	  velt2 += MT.vel[iMT].x*MT.vel[iMT].x;
	  
	  if (MT.direct[iMT].x>0)
	    {
	      nr++;
	      xr += MT.cm[iMT].x;
	      xr2 += (MT.cm[iMT].x*MT.cm[iMT].x);
	      
	      rxr += MT.rcm[iMT].x;
	      rxr2 += (MT.rcm[iMT].x*MT.rcm[iMT].x);
	      
	      velr += MT.vel[iMT].x;
	      velr2 += (MT.vel[iMT].x*MT.vel[iMT].x);
	    }
	  else
	    {
	      nl++;
	      xl += MT.cm[iMT].x;
	      xl2 += (MT.cm[iMT].x*MT.cm[iMT].x);
	      rxl += MT.rcm[iMT].x;
	      rxl2 += (MT.rcm[iMT].x*MT.rcm[iMT].x);
	      vell += MT.vel[iMT].x;
	      vell2 += (MT.vel[iMT].x*MT.vel[iMT].x);
	    }
	  
	  if (bound.type[1]==RFIXREF)
	    {
	      /* add contribution from image MTs */
	      float imxcm;
	      imxcm=2*reflectb0-MT.cm[iMT].x;
	      xt += imxcm;
	      xt2 += imxcm*imxcm; 
	      velt += -MT.vel[iMT].x;
	      velt2 += MT.vel[iMT].x*MT.vel[iMT].x;
	      
	      if (MT.direct[iMT].x>0)
		{
		  nl++;
		  xl += imxcm;
		  xl2 += imxcm*imxcm;
		  vell += -MT.vel[iMT].x;
		  vell2 += (MT.vel[iMT].x*MT.vel[iMT].x);
		}
	      else
		{
		  nr++;
		  xr += imxcm;
		  xr2 += imxcm*imxcm;
		  velr += -MT.vel[iMT].x;
		  velr2 += (MT.vel[iMT].x*MT.vel[iMT].x);
		}
	    }
	  
	}
    }

  if (nr>0) 
    {
      xr /= nr; xr2 /= nr; rxr /= nr; rxr2 /= nr; 
      velr /= nr; velr2 /= nr;
    }

  if (nl>0) 
    {
      xl /= nl; xl2 /= nl; rxl /= nl; rxl2 /= nl; 
      vell /= nl; vell2 /= nl;
    }
  
  
  if (nr+nl>0)
    {
      xt /= (nr+nl);
      xt2 /= (nr+nl);
      
      rxt /= (nr+nl);
      rxt2 /= (nr+nl);
      
      velt /= (nr+nl);
      velt2 /= (nr+nl);
    }

  av.xr[iter] += xr;
  av.xr2[iter] += xr2;
  av.xl[iter] += xl;
  av.xl2[iter] += xl2;
  av.xt[iter] += xt;
  av.xt2[iter] += xt2;

  av.rxr[iter] += rxr;
  av.rxr2[iter] += rxr2;
  av.rxl[iter] += rxl;
  av.rxl2[iter] += rxl2;
  av.rxt[iter] += rxt;
  av.rxt2[iter] += rxt2;
  if (rxt2<0) getchar();

  av.velr[iter] += velr;
  av.velr2[iter] += velr2;
  av.vell[iter] += vell;
  av.vell2[iter] += vell2;
  av.velt[iter] += velt;
  av.velt2[iter] += velt2;

  /* printf("Update_av:> nr=%f nl=%f nr+nl=%f MT.number=%d: vr=%f vl=%f vt=%f", */
/* 	 nr,nl,nr+nl,MT.number,av.velr[iter],av.vell[iter],av.velt[iter]); */
/*   getchar(); */



/* ---------------------------------------------------------------------*/
  /* Calculate av.x2 such that the time MT.iter0[iMT] that iMT enters 
     into the system is properly taken into account  */

  /* determine the number of MTs at iteration 'iter', that are 
     'tau' iterations in the system. -- a normalization factor 
     for the calulation of the averages */

  if (iter>neqsteps)
    {
      for (iMT=1; iMT<=MT.number; iMT++)
	{
	  /* if (MT.index[iMT]>200 && MT.index[iMT]<210){ */
	  /* We start evaluating the RMSD after several equilibration steps */
	  /* if the MT enters after this period, tau is time it is in the 
	     stsyem, i.e., iter-MT.iter0; if the MT enters before this, 
	     tau is iter-; also, rcm starts evaluated only after
	     neqsteps, see newpos.c.
	  */
	  if (MT.iter0[iMT]>neqsteps)
	    tau=iter-MT.iter0[iMT];
	  else
	    tau=iter-neqsteps;
	  

	  if (MT.direct[iMT].x>0) 
	    {
	      sum.nr[tau]++;
	      sum.xr[tau]+=MT.rcm[iMT].x;
	      sum.xr2[tau]+=(MT.rcm[iMT].x*MT.rcm[iMT].x);
	    }
	  else 
	    {
	      sum.nl[tau]++;
	      sum.xl[tau]+=MT.rcm[iMT].x;
	      sum.xl2[tau]+=(MT.rcm[iMT].x*MT.rcm[iMT].x);
	    } 
	  }     
	/* } */ 
    }

  /*--------------------------------------------------------------  */
  /* calculate the average number of neigbores (coordination number)*/
  /* for each MT in the ensemble                                    */
  /* we preform the summation over MT segments because each segment 
     represents a different MT from either the actual or the image 
     system.                                                        */
 
  nt=0;
  ov=0.0;
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      iseg1=2*iMT-1;
      iseg2=2*iMT;
      if (REF) 
	if (seg.length[iseg2]>0.0) {nt++; ov+=2*seg.length[iseg2];}
      for (jMT=1;jMT<=MT.number;jMT++)
	{
	  jseg1=2*jMT-1;
	  jseg2=2*jMT;
	  if (REF)
	    {
	      if (seg.length[iseg2]>0.0 && seg.length[jseg2]>0.0 && iMT!=jMT)
		{
		  /*  ---------iS1----------|----->                  */
		  /*                  <-iS2-|----------------------- */
		  /*                        |                        */
		  /* ----------jS1----------|----->                  */
		  /*                  <-jS2-|-----------------------  */
		
		  
		  if ((sz=soverlap(iseg1,jseg1))>0 
			&& ovlp.type[iseg1][jseg1]!=ZERO) {nt++; ov+=sz;}
		      
		  if ((sz=soverlap(iseg1,jseg2))>0 && 
		      ovlp.type[iseg1][jseg2]!=ZERO) {nt++; ov+=sz;}
		}

	      else if (seg.length[iseg2]>0.0 && seg.length[jseg2]<=0.0 
		       && iMT!=jMT)
		{
		  /*  ---------iS1----------|----->                  */
		  /*                  <-iS2-|----------------------- */
		  /*                        |                        */
		  /* ----------jS1--------> | <-------------------   */
		  
		  if ((sz=soverlap(iseg1,jseg1))>0 && 
		    ovlp.type[iseg1][jseg1]!=ZERO) {nt++;ov+=sz;}
		    
		      
		  if ((sz=soverlap(iseg2,jseg1))>0 && 
			ovlp.type[iseg2][jseg1]!=ZERO) {nt++;ov+=sz;}
		}

	      else if (seg.length[iseg2]<=0.0 && seg.length[jseg2]>0.0)
		{
		  /*  ---------jS1----------|----->                  */
		  /*                  <-jS2-|----------------------- */
		  /*                        |                        */
		  /* ----------iS1--------> | <-------------------   */

		 if ((sz=soverlap(iseg1,jseg1))>0 && 
		       ovlp.type[iseg1][jseg1]!=ZERO) {nt++;ov+=sz;}
		 
		 if ((sz=soverlap(iseg1,jseg2))>0 && 
		       ovlp.type[iseg1][jseg2]!=ZERO) {nt++;ov+=sz;}
		   
		}
	      else if (seg.length[iseg2]<=0.0 && seg.length[jseg2]<=0.0 
		       && iMT!=jMT)
		{
		  
		  if ((sz=soverlap(iseg1,jseg1))>0 
		       && ovlp.type[iseg1][jseg1]!=ZERO) {nt++;ov+=sz;}
		}
	    }
	  else /* regular case */
	    {
	      if ((sz=soverlap(iseg1,jseg1))>0 && ovlp.type[iseg1][jseg1]!=ZERO 
		   && iMT!=jMT)
		{
		  nt++;
		  ov+=sz;
		}

	    }
	}
      
      /* z=0; */
/*       zeros=0; */
/*       for (jMT=1;jMT<=MT.number;jMT++) */
/* 	if(overlap(iMT,jMT)>0) z=1; */
/*       if (z==0) zeros++; */
    }

  if (AvNeighboredMTs) 
    {
      if ((MT.number-zeros)>0)
	{
	  nt/=(MT.number-zeros);
	  ov/=(MT.number-zeros);
	}
      else
	{
	  /* printf("MT.number=%d zeros=%d nt=%f\n",MT.number,zeros,nt);getchar(); */
	  nt=0.0;
	  ov=0.0;
	  
	}
	
    }
  else
    {
      nt/=(MT.number);
      ov/=(MT.number);
    } 
  av.neighbors[iter]+=nt;
  av.ovlp[iter]+=ov;
	  

  /*-----------------------------------------------------------*/
  /* Calculate distance velocity/velocity correlation function */

  dx=fabs(MT.cm[MT.imax].x-MT.cm[MT.imin].x)/(nbox-1);
  
  for(ibox=0;ibox<nbox;ibox++) boxn[ibox]=0.0;
  for(ibox=0;ibox<nbox;ibox++) polarcor[ibox]=0.0;
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      xi=MT.cm[iMT].x;
      for (jMT=1;jMT<=MT.number;jMT++)
	{
	  if (iMT!=jMT)
	    {
	      xj=MT.cm[jMT].x;
	      x=fabs(xj-xi);
	      ibox=(int)rint(x/dx);
	      if (ibox>=0 && ibox<nbox)
		{
		  polarcor[ibox] += (MT.vel[iMT].x*MT.vel[jMT].x);
		  boxn[ibox]++;
		}
	    }
	}
    }

  for (ibox=0;ibox<nbox;ibox++) 
    if (boxn[ibox]>0)
      av.velcor[iter][ibox] += polarcor[ibox]/boxn[ibox];


  /*-----------------------------------------------------*/
  /* Calculate polarity correlation function over distance */
  /* C(x)=<n(0)n(x)> where n is the polarity of the MT, 
     and x is the sum is over all MTs.                     */  
  
    float dL=(rbound0-lbound0)/6.0;
  dx=fabs(MT.cm[MT.imax].x-MT.cm[MT.imin].x)/(nbox-1);
  
  for(ibox=0;ibox<nbox;ibox++) boxn[ibox]=0.0;
  for(ibox=0;ibox<nbox;ibox++) polarcor[ibox]=0.0;
  for (iMT=1;iMT<=MT.number;iMT++){
        xi=MT.cm[iMT].x;
        for (jMT=1;jMT<=MT.number;jMT++){
            if (iMT!=jMT){
                xj=MT.cm[jMT].x;
                x=fabs(xj-xi);
                ibox=(int)rint(x/dx);
                if (ABSORBING){
                    if ((xi>(lbound0+dL) && xi<(rbound0-dL)) && (xj>(lbound0+dL) && xj<(rbound0-dL)))  {
                        boxn[ibox]++;
                        polarcor  [ibox] += (MT.direct[iMT].x*MT.direct[jMT].x);
                    }
                }
                else{
                    boxn[ibox]++;
                    polarcor[ibox] += (MT.direct[iMT].x*MT.direct[jMT].x);
                }
                if (PERIOD){
                    if (nbox-ibox-1>=0)
                        polarcor[nbox-ibox-1] +=
                        (MT.direct[iMT].x*MT.direct[jMT].x);
                    if (nbox-ibox-1>=0)
                        boxn[nbox-ibox-1]++;
                }
            }
        }
  }

  for (ibox=0;ibox<nbox;ibox++) 
    if (boxn[ibox]>0)
      av.polarcor[iter][ibox] += polarcor[ibox]/boxn[ibox];
  
  /*-----------------------------------------------------*/
  /* Calculate radial (in the YZ plane) polarity correlation 
     function over distance */
  /* C(Ryz)=<n(0)n(Ryz)> where n is the polarity of the MT, 
     and x is the sum is over all MTs.                     */  

  {
    float ymin,ymax,zmin,zmax,rmin,rmax,r,dr,zi,zj,yi,yj;
    
  ymin=MT.cm[MT.iymin].y;
  ymax=MT.cm[MT.iymax].y;
  zmin=MT.cm[MT.izmin].z;
  zmax=MT.cm[MT.izmax].z;
  rmin=0.0;
  rmax=sqrt((ymax-ymin)*(ymax-ymin)+(zmax-zmin)*(zmax-zmin));
  rmax=10*2*exclude;
  dr=fabs(rmax-rmin)/(nbox2-1);
  
  for(ibox=0;ibox<nbox2;ibox++) box_yz[ibox]=0.0;
  for(ibox=0;ibox<nbox2;ibox++) polarcor_ryz[ibox]=0.0;
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      yi=MT.cm[iMT].y;
      zi=MT.cm[iMT].z;
      for (jMT=1;jMT<=MT.number;jMT++)
	{
	  if (xoverlap(iMT,jMT)>0 && iMT!=jMT)
	    {
	      yj=MT.cm[jMT].y;
	      zj=MT.cm[jMT].z;
	      r=sqrt((yj-yi)*(yj-yi)+(zj-zi)*(zj-zi));
	      ibox=(int)rint(r/dr);
	      polarcor_ryz[ibox] += (MT.direct[iMT].x*MT.direct[jMT].x);
	      box_yz[ibox]++;
	      /* if (PERIOD)  */
/* 		{ */
/* 		  if (nbox2-ibox-1>=0)  */
/* 		    polarcor_ryz[nbox2-ibox-1] +=  */
/* 		      (MT.direct[iMT].x*MT.direct[jMT].x); */
/* 		  if (nbox2-ibox-1>=0)  */
/* 		    box_yz[nbox2-ibox-1]++; */
/* 		} */
	    }
	}
    }

  for (ibox=0;ibox<nbox2;ibox++) 
    if (box_yz[ibox]>0)
      av.polarcor_ryz[iter][ibox] += polarcor_ryz[ibox]/box_yz[ibox];
  }



  /*-----------------------------------------------------*/
  /* Average Radius of Gyration OF MT bundle in YZ plane */

/*   { */
/*     float r2; */
/*     r2=0.0; */
/*     for (iMT=1;iMT<=MT.number;iMT++) */
/*       { */
/* 	float yi,zi,ri; */
/* 	yi=MT.cm[iMT].y; */
/* 	zi=MT.cm[iMT].z; */
	
/* 	r2 += sqrt(yi*yi+zi*zi); */
/*       } */
/*     r2 /= MT.number; */
/*     av.ryz[iter] += r2; */
/*   } */

  {
    float rav,ycm,zcm,yi,zi;
    rav=0.0;
    ycm=0;
    zcm=0;
    for (iMT=1;iMT<=MT.number;iMT++)
      {
	ycm += MT.cm[iMT].y;
	zcm += MT.cm[iMT].z;
      }
    ycm /= MT.number;
    zcm /= MT.number;
    
    for (iMT=1;iMT<=MT.number;iMT++)
      {
	yi = MT.cm[iMT].y;
	zi = MT.cm[iMT].z;
	rav += sqrt((yi-ycm)*(yi-ycm)+(zi-zcm)*(zi-zcm));
      }
    rav /= MT.number;
    
    av.ryz[iter] += rav;
  }

  /* ------------------------------------------------------------ */
  /* Calculate bundle shape */


  for (ibox=0;ibox<nbox;ibox++)
    {
      rmin[ibox]=bundle.rad;
      rmax[ibox]=-bundle.rad;
    }
  
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      ibox=(int)rint((MT.cm[iMT].x-MT.cm[MT.imin].x)/dx); 
      /*      r=sqrt(MT.cm[iMT].z*MT.cm[iMT].z+MT.cm[iMT].y*MT.cm[iMT].y);*/
      
      r=MT.cm[iMT].z;

      if (r<rmin[ibox]) 
	rmin[ibox]=r;
      if (r>rmax[ibox])
	rmax[ibox]=r;
    }
  
  for (ibox=0;ibox<nbox;ibox++)
    {
      av.rmin[iter][ibox] += rmin[ibox];
      av.rmax[iter][ibox] += rmax[ibox];
    }
  
  /*-----------------------------------------------------------------------*/
  /* GET ALL TRACK DATA */
  //Assaf output
    /*##############################################*/
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      track.vel[iter][iMT]=MT.vel[iMT].x;   /* for time correlations */
      track.index[iter][iMT]=MT.index[iMT]; /* for time correlations */
      track.xcm[iter][iMT]=MT.rcm[iMT].x;
      track.ycm[iter][iMT]=MT.rcm[iMT].y;
      track.zcm[iter][iMT]=MT.rcm[iMT].z;
      
      track.it0[iter][iMT]=MT.iter0[iMT];
    }
        /*##############################################*/
    //MaxOUTPUT
        /*##############################################*/
    
    ////////GET BOUNDARIES FOR VISUALS///////////////////////////////////////
    if (bound.type[0]==PERIODIC || bound.type[0]==WALL
        || bound.type[0]==STICKY || bound.type[0]==NET
        || bound.type[0]==ABSORB)
        track.lbound[iter]=lbound0;
    else track.lbound[iter]=MT.xmin_real-0.5*MT.length[MT.imin_real];
    
    if (bound.type[1]==NON || bound.type[1]==FORCE) track.rbound[iter]=MT.xmax_real+0.5*MT.length[MT.imax_real];
    else                        track.rbound[iter]=rbound0;
     /////////////////////////////////////////////////////////////////////////
    
    /////////////////GET ARROW TRACK DATA (FOR VISUALIZATION)////////////////
    if (D3_output==TRUE) track.arrow_number[iter]=get_mathematica_arrow3D();
    else  track.arrow_number[iter]=get_mathematica_arrow2D();
    /////////////////////////////////////////////////////////////////////////
    /////////////////GET COORDS FOR VISUALS//////////////////////////////////
    track.nMT[iter]=MT.number;

    for (iMT=1; iMT<=2*MT.number; iMT++){
        track.xcm_vis[iter][iMT] = seg.cm[iMT].x;
        track.ycm_vis[iter][iMT] = seg.cm[iMT].y;
        track.zcm_vis[iter][iMT] = seg.cm[iMT].z;
        track.length[iter][iMT] = seg.length[iMT];
        track.direct[iter][iMT]=seg.direct[iMT].x;
    }

  /*-----------------------------------------------------------------------*/
  
  
  free_ivector(neighbors,1,MT.number);
}
  
