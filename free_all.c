/* #################################################################### */
/* free_all -- frees allocated memory after each av it and in the end   */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global_var.h"
#include "nrutil.h"

void free_MTdat(){
    free_vector(MT.length,1,MT.number);
    free_pvector(MT.cm,1,MT.number);
    free_pvector(MT.rcm,1,MT.number);
    free_pvector(MT.mend,1,MT.number);
    free_pvector(MT.pend,1,MT.number);
    free_pvector(MT.direct,1,MT.number);
    free_pvector(MT.vel,1,MT.number);
    free_ivector(MT.index,1,MT.number);
    free_ivector(MT.cluster,0,MT.number);
    free_ivector(MT.is_in_net,1,MT.number);
    free_ivector(MT.grow_direct,1,MT.number);
    free_ivector(MT.iter0,1,MT.number);
    free_imatrix((int **)ovlp.type,1,2*MT.number,1,2*MT.number);
    free_imatrix(ovlp.motor_direction,1,2*MT.number,1,2*MT.number);
    free_matrix(ovlp.iter0,1,2*MT.number,1,2*MT.number);
    free_vector( dist_dyn,1,nbox);
    free_vector( dist_kin,1,nbox);
    free_vector( ProbBi_dist,1,nbox);
    free_vector( ProbKi_dist,1,nbox);
    free_vector( ProbDy_dist,1,nbox);
    free_vector( ProbAct_dist,1,nbox);
    free_ivector( ovlp.actin, 1, MT.number);
    free_ivector( ovlp.actin_it0, 1, MT.number);
    
    /* free all segment memory */
    free_vector(seg.length,1,2*MT.number);
    free_pvector(seg.cm,1,2*MT.number);
    free_pvector(seg.mend,1,2*MT.number);
    free_pvector(seg.pend,1,2*MT.number);
    free_pvector(seg.direct,1,2*MT.number);
}

void free_all(){
    /*free av structure*/
    free_matrix( av.nr ,1,niter,1,nbox);
    free_matrix( av.nl ,1,niter,1,nbox);
    free_matrix( av.rnr ,1,niter,1,nbox);
    free_matrix( av.rnl ,1,niter,1,nbox);
    free_matrix( av.rmin ,1,niter,1,nbox);
    free_matrix( av.rmax ,1,niter,1,nbox);
    free_vector( av.ntot ,1,niter);
    free_vector( av.xmin ,1,niter);
    free_vector( av.xmax ,1,niter);
    free_vector( av.rxmin ,1,niter);
    free_vector( av.rxmax ,1,niter);
    free_vector( av.zmin ,1,niter);
    free_vector( av.zmax ,1,niter);
    free_vector( av.ymin ,1,niter);
    free_vector( av.ymax ,1,niter);
    free_vector( av.xr ,1,niter);
    free_vector( av.xr2 ,1,niter);
    free_vector( av.xl ,1,niter);
    free_vector( av.xl2 ,1,niter);
    free_vector( av.xt, 1,niter);
    free_vector( av.xt2 ,1,niter);
    free_vector( av.rxr ,1,niter);
    free_vector( av.rxr2 ,1,niter);
    free_vector( av.rxl   ,1,niter);
    free_vector( av.rxl2   ,1,niter);
    free_vector( av.rxt   ,1,niter);
    free_vector( av.rxt2   ,1,niter);
    free_vector( av.velt   ,1,niter);
    free_vector( av.vell ,1,niter);
    free_vector( av.velr ,1,niter);
    free_vector( av.velt2 ,1,niter);
    free_vector( av.vell2 ,1,niter);
    free_vector( av.velr2 ,1,niter);
    free_vector( av.neighbors ,1,niter);
    free_vector( av.ovlp ,1,niter);
    free_vector( av.absvmin ,1,niter);
    free_vector( av.absvmax ,1,niter);
    free_vector( av.vmin ,1,niter);
    free_vector( av.vmax ,1,niter);
    free_vector( av.absvel_ntot ,1,niter);
    free_vector( av.vel_ntot ,1,niter);
    free_matrix( av.absvel_nr ,1,niter, 1,nbox);
    free_matrix( av.absvel_nl ,1,niter, 1,nbox);
    free_matrix( av.vel_nr ,1,niter, 1,nbox);
    free_matrix( av.vel_nl ,1,niter, 1,nbox);
    free_matrix( av.xvelr ,1,niter, 1,nbox);
    free_matrix( av.xvell ,1,niter, 1,nbox);
    free_matrix( av.polarcor ,1,niter, 1,nbox);
    free_matrix( av.polarcor_ryz ,1,niter, 1,nbox);
    free_matrix( av.velcor ,1,niter, 1,nbox);
    free_vector( av.nMT ,1,niter);
    free_vector( av.nMTr ,1,niter);
    free_vector( av.nMTl ,1,niter);
    free_vector( av.nzeros ,1,niter);
    free_vector( av.ryz ,1,niter);
    free_vector( av.op ,1,niter);
    
    free_matrix(av.velcort,0,corr_tf+1,0,corr_tf-corr_ti+1);
    free_matrix(av.velcort_pp,0,corr_tf+1,0,corr_tf-corr_ti+1);
    free_matrix(av.velcort_pm,0,corr_tf+1,0,corr_tf-corr_ti+1);
    free_matrix(av.velcort_mm,0,corr_tf+1,0,corr_tf-corr_ti+1);
    free_matrix(av.tubulin_dens,0,niter,0,nbox);
    free_matrix(av.ovlp_kin,1,niter,1,nbox);
    free_matrix(av.ovlp_dyn,1,niter,1,nbox);
    free_matrix(av.ovlp_act,1,niter,1,nbox);
    free_matrix(av.ovlp_bi,1,niter,1,nbox);
    free_matrix(av.nv_cluster,1,mxmts,1,nbox);
    free_matrix(av.tubulin_r,1,niter,1,nbox);
    free_matrix(av.tubulin_l,1,niter,1,nbox);
    free_vector(av.xltau,0,niter);
    free_vector(av.xl2tau,0,niter);
    free_vector(av.xrtau,0,niter);
    free_vector(av.xr2tau,0,niter);
    free_vector(av.xtau,0,niter);
    free_vector(av.x2tau,0,niter);
    free_vector(av.nltau,0,niter);
    free_vector(av.nrtau,0,niter);
    free_vector(av.xmax_real,0,niter);
    free_vector(av.xmin_real,0,niter);
    free_matrix(av.thick,1,niter,1,nbox);
    free_vector(av.ovlpsize_perOV,1,niter);
    free_vector(av.ovlpsize_perMT,1,niter);
    
    /*free the track structure*/
    free_ivector(track.nMT,1,niter);
    free_imatrix( track.index,1,niter,1,mxmts);
    free_ivector( track.index_count,1,niter*mxmts);
    free_matrix( track.vel,1,niter,1,mxmts);
    free_matrix( track.xcm,1,niter,1,2*mxmts);
    free_matrix( track.ycm,1,niter,1,2*mxmts);
    free_matrix( track.zcm,1,niter,1,2*mxmts);
    free_matrix( track.length,1,niter,1,2*mxmts);
    free_matrix( track.xcm_vis,1,niter,1,2*mxmts);
    free_matrix( track.ycm_vis,1,niter,1,2*mxmts);
    free_matrix( track.zcm_vis,1,niter,1,2*mxmts);
    free_matrix( track.mathematica_arrows,1,niter,1,mxmotors*12);
    free_vector( track.arrow_number,1,niter);
    free_vector( track.lbound,1,niter);
    free_vector( track.rbound,1,niter);
    free_matrix( track.direct,1,niter,1,2*mxmts);
    free_imatrix( track.it0,1,niter,1,2*mxmts);
    free_imatrix( track.cross_type,1,niter,1,mxmotors);
    
}
