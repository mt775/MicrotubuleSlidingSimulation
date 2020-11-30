/* #################################################################### */
/* initialize_variables -- for each av iteration variables get created  */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global_var.h"
#include "nrutil.h"

/*Initializes all averaging variables before first av loop starts*/
void initialize_av_variables(){

    int iMT,ibox;
    long it1,it2,idt,tau;
    /* Initiate average vectors and matrices      */
    /* iter,niter -- time steps                   */
    /* iav, nav -- ensemble averaging steps       */
    /* nbox -- number of grid boxes in x-profiles */
    
    
    /*Initialize the track variables They are for tracking certain MT properties
     mainly for index tracking and later visualization*/
    /*We need to guess the maximum number of MTs to allocate a track struct which is large enough
     to incorporate all MTs at every timepoint*/

    if (newMTs0==0) {
        mxmts=MTfirst;
        /*As we have a hexagonal lattice each MT can have 6 neighbors maximum*/
        mxmotors=mxmts*32;
    }
    else if( newMTs0!=0 ){
        mxmts=(niter/newfreq*newMTs0)+MTfirst+2;
	mxmts=mxmts; /*Account for outflow MJ 14/08/19*/
        mxmotors=mxmts*12;
    }
    else{
        printf("Error in initialize_variables.c\n");
        exit(0);
    }

    
    track.nMT=ivector(1,niter);
    track.index=imatrix(1,niter,1,mxmts);
    track.index_count=ivector(1,niter*mxmts);
    track.vel=matrix(1,niter,1,mxmts);
    track.xcm=matrix(1,niter,1,2*mxmts);
    track.ycm=matrix(1,niter,1,2*mxmts);
    track.zcm=matrix(1,niter,1,2*mxmts);
    track.length=matrix(1,niter,1,2*mxmts);
    track.xcm_vis=matrix(1,niter,1,2*mxmts);
    track.ycm_vis=matrix(1,niter,1,2*mxmts);
    track.zcm_vis=matrix(1,niter,1,2*mxmts);
    track.mathematica_arrows=matrix(1,niter,1,mxmotors*12);
    track.arrow_number=vector(1,niter);
    track.lbound=vector(1,niter);
    track.rbound=vector(1,niter);
    track.direct=matrix(1,niter,1,2*mxmts);
    track.it0=imatrix(1,niter,1,2*mxmts);
    track.cross_type=imatrix(1,niter,1,mxmotors);
    
    
    

    // define matrices and vectors for av. structure
    av.nr =matrix(1,niter,1,nbox);
    av.nl =matrix(1,niter,1,nbox);
    av.rnr =matrix(1,niter,1,nbox);
    av.rnl =matrix(1,niter,1,nbox);
    av.rmin =matrix(1,niter,1,nbox);
    av.rmax =matrix(1,niter,1,nbox);
    av.ntot =vector(1,niter);
    av.xmin =vector(1,niter);
    av.xmax =vector(1,niter);
    av.rxmin =vector(1,niter);
    av.rxmax =vector(1,niter);
    av.zmin =vector(1,niter);
    av.zmax =vector(1,niter);
    av.ymin =vector(1,niter);
    av.ymax =vector(1,niter);
    av.xr =vector(1,niter);
    av.xr2 =vector(1,niter);
    av.xl =vector(1,niter);
    av.xl2 =vector(1,niter);
    av.xt =vector(1,niter);
    av.xt2 =vector(1,niter);
    av.rxr =vector(1,niter);
    av.rxr2 =vector(1,niter);
    av.rxl =vector(1,niter);
    av.rxl2 =vector(1,niter);
    av.rxt =vector(1,niter);
    av.rxt2 =vector(1,niter);
    av.velt =vector(1,niter);
    av.velr =vector(1,niter);
    av.vell =vector(1,niter);
    av.velt2 =vector(1,niter);
    av.vell2 =vector(1,niter);
    av.velr2 =vector(1,niter);
    av.neighbors =vector(1,niter);
    av.ovlp =vector(1,niter);
    av.absvmin =vector(1,niter);
    av.absvmax =vector(1,niter);
    av.vmin =vector(1,niter);
    av.vmax =vector(1,niter);

    av.absvel_ntot =vector(1,niter);
    av.vel_ntot =vector(1,niter);
    av.absvel_nr =matrix(1,niter, 1,nbox);
    av.absvel_nl =matrix(1,niter, 1,nbox);
    av.vel_nr =matrix(1,niter, 1,nbox);
    av.vel_nl =matrix(1,niter, 1,nbox);
    av.xvelr =matrix(1,niter, 1,nbox);
    av.xvell =matrix(1,niter, 1,nbox);
    av.polarcor =matrix(1,niter, 1,nbox);
    av.polarcor_ryz =matrix(1,niter, 1,nbox);
    av.velcor =matrix(1,niter, 1,nbox);
    av.nMT =vector(1,niter);
    av.nMTr =vector(1,niter);
    av.nMTl =vector(1,niter);
    av.nzeros =vector(1,niter);
    av.ryz =vector(1,niter);
    av.op =vector(1,niter);
    
    av.velcort=matrix(0,corr_tf,0,corr_tf-corr_ti);
    av.velcort_pp=matrix(0,corr_tf,0,corr_tf-corr_ti);
    av.velcort_pm=matrix(0,corr_tf,0,corr_tf-corr_ti);
    av.velcort_mm=matrix(0,corr_tf,0,corr_tf-corr_ti);
    av.ovlp_kin=matrix(1,niter, 1,nbox);
    av.ovlp_dyn=matrix(1,niter, 1,nbox);
    av.ovlp_act=matrix(1,niter, 1,nbox);
    av.ovlp_bi=matrix(1,niter, 1,nbox);
    av.tubulin_r=matrix(1,niter, 1,nbox);
    av.tubulin_l=matrix(1,niter, 1,nbox);
    av.tubulin_dens=matrix(0,niter,0,nbox);
    av.nv_cluster=matrix(1,mxmts,1,nbox);
    av.xltau=vector(0,niter);
    av.xl2tau=vector(0,niter);
    av.xrtau=vector(0,niter);
    av.xr2tau=vector(0,niter);
    av.xtau=vector(0,niter);
    av.x2tau=vector(0,niter);
    av.nrtau=vector(0,niter);
    av.nltau=vector(0,niter);
    av.xmax_real=vector(0,niter);
    av.xmin_real=vector(0,niter);
    av.thick=matrix(1,niter,1,nbox);
    av.ovlpsize_perMT=vector(1,niter);
    av.ovlpsize_perOV=vector(1,niter);
    

    
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    // Set all vectors and matrices to zero!!
    for (iter=1;iter<=niter;iter++)
    {
        

        for (ibox=0;ibox<nbox;ibox++)
        {

            av.nr[iter][ibox]=0;
            av.nl[iter][ibox]=0;
            av.rnr[iter][ibox]=0;
            av.rnl[iter][ibox]=0;
            av.xvelr[iter][ibox]=0;
            av.xvell[iter][ibox]=0;
            av.tubulin_dens[iter][ibox+1]=0;
            av.tubulin_l[iter][ibox]=0;
            av.tubulin_r[iter][ibox]=0;
            av.ovlp_kin[iter][ibox]=0;
            av.ovlp_dyn[iter][ibox]=0;
            av.ovlp_act[iter][ibox]=0;
            av.ovlp_bi[iter][ibox]=0;
            av.absvel_nr[iter][ibox]=0;
            av.absvel_nl[iter][ibox]=0;
            av.vel_nr[iter][ibox]=0;
            av.vel_nl[iter][ibox]=0;
            av.velcor[iter][ibox]=0;
            av.polarcor[iter][ibox]=0;
            av.rmin[iter][ibox]=0;
            av.rmax[iter][ibox]=0;
            av.thick[iter][ibox+1]=0.0;
            
            
        }
        
        for (ibox=1;ibox<nbox;ibox++)
            av.polarcor_ryz[iter][ibox]=0;
        
        
        for (iMT=1;iMT<=mxmts;iMT++)
        {
            track.vel[iter][iMT]=0.0;   /* for time correlations */
            track.index[iter][iMT]=0.0; /* for time correlations */

            
        }

        for (iMT=1; iMT<= mxmotors*12; iMT++) {
            track.mathematica_arrows[iter][iMT]=0.0;
        }

        for (iMT=1; iMT<= mxmts*2; iMT++) {
            track.direct[iter][iMT]=0.0;
            track.length[iter][iMT]=0.0;
            
            track.xcm[iter][iMT]=0.0;
            track.ycm[iter][iMT]=0.0;    /* for calculating trajectories */
            track.zcm[iter][iMT]=0.0;
            // Max edit mathematica tracks
            track.xcm_vis[iter][iMT] = 0.0;
            track.ycm_vis[iter][iMT] = 0.0;
            track.zcm_vis[iter][iMT] = 0.0;
        }
        for (iMT=1; iMT<= mxmotors; iMT++) {
            track.cross_type[iter][iMT]=0.0;
        }
    
        track.nMT[iter]=0.0;
        track.arrow_number[iter]=0.0;
        
        
        av.op[iter]=0.0;
        av.ntot[iter]=0.0;
        av.absvel_ntot[iter]=0.0;
        av.vel_ntot[iter]=0.0;
        av.xmin[iter]=0.0;
        av.xmax[iter]=0.0;
        av.ymin[iter]=0.0;
        av.ymax[iter]=0.0;
        av.zmin[iter]=0.0;
        av.zmax[iter]=0.0;
        av.rxmin[iter]=0.0;
        av.rxmax[iter]=0.0;
        av.absvmin[iter]=0.0;
        av.absvmax[iter]=0.0;
        av.vmin[iter]=0.0;
        av.vmax[iter]=0.0;
        av.xmax_real[iter]=0.0;
        av.xmin_real[iter]=0.0;
        av.xr[iter]=0.0;
        av.xr2[iter]=0.0;
        av.xl[iter]=0.0;
        av.xl2[iter]=0.0;
        av.xt[iter]=0.0;
        av.xt2[iter]=0.0;
        av.rxr[iter]=0.0;
        av.rxr2[iter]=0.0;
        av.rxl[iter]=0.0;
        av.rxl2[iter]=0.0;
        av.rxt[iter]=0.0;
        av.rxt2[iter]=0.0;
        av.velr[iter]=0.0;
        av.vell[iter]=0.0;
        av.velt[iter]=0.0;
        av.velr2[iter]=0.0;
        av.vell2[iter]=0.0;
        av.velt2[iter]=0.0;
        av.neighbors[iter]=0.0;
        av.ovlp[iter]=0.0;
        av.nMT[iter]=0.0;
        av.nMTr[iter]=0.0;
        av.nMTl[iter]=0.0;
        av.nzeros[iter]=0.0;
        av.ryz[iter]=0.0;
        
        av.ovlpsize_perOV[iter]=0.0;
        av.ovlpsize_perMT[iter]=0.0;
    }
    
    for (iMT=1; iMT<=mxmts; iMT++) {
        for (ibox=1; ibox<=nbox; ibox++) {
            av.nv_cluster[iMT][ibox]=0;
        }
    }
    
    
    for(it1=corr_ti;it1<=corr_tf;it1++)
    {
        for(it2=it1;it2<=corr_tf;it2++)
        {
            idt=it2-it1;
            av.velcort[it1][idt]=0.0;
            av.velcort_pp[it1][idt]=0.0;
            av.velcort_mm[it1][idt]=0.0;
            av.velcort_pm[it1][idt]=0.0;
        }
        
    }
    
    for (tau=0; tau<=niter; tau++)
    {
        av.xrtau[tau] =0.0;
        av.xr2tau[tau] =0.0;
        av.xltau[tau] =0.0;
        av.xl2tau[tau] =0.0;
        av.xtau[tau] =0.0;
        av.x2tau[tau] =0.0;
    }

}

/*INitializes variables for each av iteration*/
void initialize_it_variables(int newMTs){
    int iMT, ibox;
    long tau;
    
    perc_prob=0.;
    printf("init>>MT.number =%i\n",newMTs);
    MT.cm=pvector(1,newMTs);
    MT.rcm=pvector(1,newMTs);
    MT.mend=pvector(1,newMTs);
    MT.pend=pvector(1,newMTs);
    MT.length=vector(1,newMTs);
    MT.direct=pvector(1,newMTs);
    MT.vel=pvector(1,newMTs);
    MT.index=ivector(1,newMTs);
    emptyint(MT.index,1,newMTs);
    MT.iter0=ivector(1,newMTs);
    emptyint(MT.iter0,1,newMTs);
    MT.cluster=ivector(0,newMTs);
    emptyint(MT.cluster,0,newMTs);
    
    MT.is_in_net=ivector(1,newMTs);
    emptyint(MT.is_in_net,1,newMTs);
    MT.grow_direct=ivector(1,newMTs);
    emptyint(MT.grow_direct,1,newMTs);
    
    ovlp.type=(enum motortype **)imatrix(1,2*newMTs,1,2*newMTs); //creates an overlap matrix
    ovlp.iter0=matrix(1,2*newMTs,1,2*newMTs);                     // each filament has two segments
    ovlp.motor_direction=imatrix(1,2*newMTs,1,2*newMTs);
    
    /*Allocate memory for actin ov data*/
    ovlp.actin=ivector(1,newMTs);
    ovlp.actin_it0=ivector(1,newMTs);
    emptyint(ovlp.actin, 1,newMTs);
    
    /* Allocate memory for segment data */
    seg.cm=pvector(1,2*newMTs);
    seg.mend=pvector(1,2*newMTs);
    seg.pend=pvector(1,2*newMTs);
    seg.length=vector(1,2*newMTs);
    seg.direct=pvector(1,2*newMTs);
    
    
    /* Allocate memory for trajectory averages data */
    sum.xl=vector(0,niter);
    sum.xl2=vector(0,niter);
    sum.xr=vector(0,niter);
    sum.xr2=vector(0,niter);
    sum.nl=vector(0,niter);
    sum.nr=vector(0,niter);
    /* Zero single trajectory averages */
    for (tau=0;tau<=niter;tau++)
	{
        sum.xl[tau]=0.0;
        sum.xl2[tau]=0.0;
        
        sum.xr[tau]=0.0;
        sum.xr2[tau]=0.0;
        
        sum.nr[tau]=0.0;
        sum.nl[tau]=0.0;
	}

    
}
