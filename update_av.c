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

/*gets thickness of bundle see after update_av*/
float getstrands();

float av_ovsize_perOV();
float av_ovsize_perMT();
float DEBUGGER=0, DEBUGL=0;

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
  float r,rmin[nbox],rmax[nbox],boxnr[nbox],boxnl[nbox],boxn[nbox];
  float box_yz[nbox],polarcor_ryz[nbox];
  float vel_nr[nbox],vel_nl[nbox],ntot,polarcor[nbox];
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
    /*if(iter>1 && ((av.xmax_real[iter]-av.xmin_real[iter])/(float)iav-(av.xmax_real[iter-1]-av.xmin_real[iter-1])/(float)iav)>5){
        printf(" %i peak detected why!?\n",iter);
        printf("xmax=%f xmin=%f\n",av.xmax_real[iter]/(float)iav,av.xmin_real[iter]/(float)iav);
        printf("xmaxold=%f xminold=%f\n",av.xmax_real[iter-1]/(float)iav,av.xmin_real[iter-1]/(float)iav);
        exit(0);
    }*/
    
  av.zmin[iter] += MT.zmin;
  av.zmax[iter] += MT.zmax;

  av.rxmin[iter] += MT.rxmin;
  av.rxmax[iter] += MT.rxmax; 


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
    
    
  /*Calculate av cluster absolute velocities*/
    
    
    float vel_clust[MT.number+1][nbox+1];
    
    for (ibox=0;ibox<=nbox;ibox++){
        for (iMT=0; iMT<=MT.number; iMT++) {
            vel_clust[iMT][ibox]=0.;
            
        }
    }
    
    imaxvel=iabsVelmax();
    iminvel=iabsVelmin();
    //printf("max velocity=%f\n", MT.vel[imaxvel].x);
    
    dv=maxvel/(float)nbox;
    //dv=(fabs(MT.vel[imaxvel].x/*-velpol*MT.direct[imaxvel].x*/)
     //   -fabs(MT.vel[iminvel].x/*-velpol*MT.direct[iminvel].x*/))/nbox;
    if (dv!=0.0) {
        
        for (iMT=1; iMT<=MT.number; iMT++) {
            
            ibox= (int)floor(fabs(MT.vel[iMT].x/dv))+1;
            //ibox=(int)floor((fabs(MT.vel[iMT].x/*-velpol*MT.direct[iMT].x*/)
              //              -fabs(MT.vel[iminvel].x/*-velpol*MT.direct[iminvel].x*/))/dv)+1;
            if (ibox==0)
                ibox=1;
            else if (ibox==nbox+1)
                ibox=nbox;
            else if (ibox > nbox+1 ){
                /*printf("velocity out of range for updating\n");
                printf("MT.vel=%f\n", MT.vel[iMT].x);
                printf("maximum velocity allowed: %f\n", maxvel);
                printf("Hoping that you choose this parameter wisely I will put this velocity in the last velocity box\n");
                printf("If this message is displayed too often the clustvel dist is flawed\n");*/
                ibox=nbox;
            }
            if (iMT != iminvel && iMT != imaxvel && (CANCEL || NEIGHBOR_CRITERIA)){
                vel_clust[MT.cluster[iMT]][ibox]++;
                //clust_tot[MT.cluster[iMT]]++;

            }
        }
    }
    
    
    for (ibox=1;ibox<=nbox;ibox++){
        
        for (iMT=1; iMT<=MT.number; iMT++) {
            if (ibox==1 && 0) {
                
                printf("ve_clust=%f\n",vel_clust[MT.cluster[iMT]][ibox]);
                printf("MT.clust=%i\n",MT.cluster[iMT]);
                printf("mxmts=%i\n",mxmts);
                printf("avnv=%f\n",av.nv_cluster[MT.cluster[iMT]][ibox]);
            }
            if (vel_clust[MT.cluster[iMT]][ibox]>0) {
                av.nv_cluster[MT.cluster[iMT]][ibox] += vel_clust[MT.cluster[iMT]][ibox];
            }
        }
    }
    
    
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

    

    

  /*-----------------------------------------------------------------------*/
  /* Calculate Histograms of actual Velocity                             */

  for(ibox=0;ibox<nbox;ibox++){
      vel_nr[ibox]=0.0;
      vel_nl[ibox]=0.0;
  }
  ntot=0.0;
    nl=0.;
    nr=0.;
  

  /* Uncomment to use fixed values for the min/max values of Vel */ 
  /* these are defined in rdata.c                                */
  /* also, update wav_veldensity() in write_output.c             */ 

  velmin=-maxvel;
  velmax=maxvel;
  

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
	      if (ibox<0 || ibox>=nbox) /* nbox is a global variable */{
              printf("Update_av:> Error velocity out of range\n");
              printf("Update_av:> iminvel=%d imaxvel=%d\n",iminvel,imaxvel);
              printf("Update_av:> MT.vel[%d].x=%f MT.vel[%d].x=%f dv=%f\n",
                 iminvel,velmin,imaxvel,velmax,dv);
              printf("Update_av:> ibox=%d MT.vel[%d].x=%f dv=%f\n",
                 ibox,iMT,MT.vel[iMT].x,dv);
              
              for (jMT=1;jMT<=MT.number;jMT++) 
                printf("vel[%d]=%f\n",jMT,MT.vel[jMT].x);
              
              
              //getchar();
          }
	      
	      if (iMT != iminvel && iMT != imaxvel && (CANCEL || NEIGHBOR_CRITERIA))
	      /* The two extreme MTs, are excluded from the histogram calc,
		 since they are used to calculate the histogram boundaries, 
		 and therefore by their definition they will have a finite 
		 probablity in each configuration. */
	      //printf("vel=%f\n",MT.vel[iMT].x);
	      /*if (NEIGHBOR_CRITERIA)*/
		{
		  if (MT.direct[iMT].x>0)
		    {
		      ntot++;
              nr++;
		      vel_nr[ibox]++;
                
		    }
		  else
		    {
		      ntot++;
		      vel_nl[ibox]++;
              nl++;
            }
		  
		}
	      
	    }
	}
    }

 /* normalize */
  for (ibox=0;ibox<nbox;ibox++)
    {
      if (nl>0)
          av.vel_nl[iter][ibox]+=(vel_nl[ibox]/(float)nl);
	  if (nr>0)
          av.vel_nr[iter][ibox]+=(vel_nr[ibox]/(float)nr);
	  
    }

    DEBUGGER=0;
    DEBUGL=(av.xmax_real[iter]-av.xmin_real[iter])/iav;
   
        /*-----------------------------------------------------*/
        /* Average number of strands */
        for (ibox=1; ibox<=nbox; ibox++){
            av.thick[iter][ibox]+=getstrands(ibox);
             DEBUGGER+=av.thick[iter][ibox]/(iav)*DEBUGL/nbox;
           // DEBUGGER+=getstrands(ibox)*DEBUGL/nbox;

        }
        //printf("length it%i=%f\n",iter,DEBUGL);
        /*if (iter==1000){
                printf("it=%i Dbugstrands=%f dx=%f L=%f\n",iter,DEBUGGER,DEBUGL/nbox,DEBUGL);
            for (ibox=1; ibox<=nbox; ibox++) {
                printf("%f\t",av.thick[iter][ibox]/iav);
            }
            printf("\n");
        }*/
   
    /*-----------------------------------------------------*/
    /* Average ovlp size in bundle */
    av.ovlpsize_perMT[iter]+=av_ovsize_perMT();
    av.ovlpsize_perOV[iter]+=av_ovsize_perOV();
    

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
    else if (bound.type[0]==NON || bound.type[0]==FORCE )
        track.lbound[iter]=MT.xmin_real-0.5*MT.length[MT.imin_real];
    else{
        printf("update_av.c >> lbound undefined!\n");
        exit(0);
    }
    
    if (bound.type[1]==NON || bound.type[1]==FORCE) track.rbound[iter]=MT.xmax_real+0.5*MT.length[MT.imax_real];
    else if (bound.type[1]==PERIODIC)        track.rbound[iter]=rbound0;
    else{
        printf("update_av.c >> lbound undefined!\n");
        exit(0);
    }
    /////////////////////////////////////////////////////////////////////////
    
    /////////////////GET ARROW TRACK DATA (FOR VISUALIZATION)////////////////
    if (D3_output==TRUE /*&& iav==1*/) track.arrow_number[iter]=get_mathematica_arrow3D();
    else if (iav==1)  track.arrow_number[iter]=get_mathematica_arrow2D();
    /////////////////////////////////////////////////////////////////////////
    /////////////////GET COORDS FOR VISUALS//////////////////////////////////
    track.nMT[iter]=MT.number;
    

    for (iMT=1; iMT<=2*MT.number; iMT++){
        track.xcm_vis[iter][iMT] = seg.cm[iMT].x;
        track.ycm_vis[iter][iMT] = seg.cm[iMT].y;
        track.zcm_vis[iter][iMT] = seg.cm[iMT].z;
        track.length[iter][iMT] = seg.length[iMT];
        track.direct[iter][iMT]=seg.direct[iMT].x;
        if (track.direct[iter][iMT] < -1.1 || track.direct[iter][iMT] > 1.1 || (track.direct[iter][iMT] > -0.1 && track.direct[iter][iMT] < 0.1) ){
            printf("iMT=%i nMT=%i seg direct=%f length=%f\n",iMT,MT.number,seg.direct[iMT].x,seg.length[iMT]);
            printf("ERROR!\n direction out of bounds?! \n iMT=%i nMT=%i",iMT,MT.number);
            exit(0);
        }
    }

  /*-----------------------------------------------------------------------*/
  
  
  free_ivector(neighbors,1,MT.number);
    
}

  /*-----------------------------------------------------------------------*/
  /*-----------------------------------------------------------------------*/
  /*-----------------------------------------------------------------------*/

/*function to calculate AvOv length per overlap*/
float av_ovsize_perOV(){
    int iMT,jMT,count=0;
    float avov=0,dummy;
    for (iMT=1; iMT<=MT.number; iMT++) {
        for (jMT=iMT+1; jMT<=MT.number; jMT++) {
            if (iMT!=jMT){
                dummy=overlap(iMT,jMT);
                if (dummy>0){
                    avov+=dummy;
                    count++;
                }
            }
        }
    }
    return(avov/(float)count);
}


/*function to calculate AvOv length per filament*/
float av_ovsize_perMT(){
    int iMT,jMT,count=0;
    float avov=0,dummy;
    for (iMT=1; iMT<=MT.number; iMT++) {
        for (jMT=iMT+1; jMT<=MT.number; jMT++) {
            if (iMT!=jMT){
                dummy=overlap(iMT,jMT);
                if (dummy>0){
                    avov+=dummy;
                }
            }
        }
    }

    return(avov/(float)MT.number);
}

float getstrands( int ibox ){
    int iMT;
    float m=0,ov;
    //float rxmax=av.xmax_real[iter]/iav;
    //float lxmin=av.xmin_real[iter]/iav;
    float rxmax=MT.xmax_real+0.5*MT.length[MT.imax_real];
    float lxmin=MT.xmin_real-0.5*MT.length[MT.imin_real];
    float dx=(rxmax-lxmin)/(float)nbox;
    //printf("r=%f l=%f\n",rxmax,lxmin);
    
    
    //if ((bound.type[0]==WALL || bound.type[0]==NET) && bound.type[1]==FORCE) {
    for (iMT=1; iMT<=MT.number; iMT++){
        ov=Interval_overlap(((ibox-1)*dx)+lxmin, MT.cm[iMT].x-(MT.length[iMT]/2), dx, MT.length[iMT]);
        
        if ( ov >0 ){
            m+=ov/dx;
        }
        
    }
    /*else if (bound.type[1]==FORCE && bound.type[1]==FORCE){
        for (iMT=1; iMT<=MT.number; iMT++){
            if ((MT.cm[iMT].x - MT.length[iMT]/2 -touch_depth <= lxmin) ||
                (MT.cm[iMT].x + MT.length[iMT]/2 +touch_depth >= rxmax)){
                m++;
            }
        }
    }*/

    return (m);
}
  
