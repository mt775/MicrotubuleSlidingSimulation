/* #################################################################### */
/* normalize  --  Normalizes average quantities                         */

#include <stdio.h>
#include "global_var.h"

float vel_nr_tot,vel_nl_tot,absvel_nr_tot,absvel_nl_tot; 

void normalize()
{
  long it,it1,it2,idt,ibox,iMT,tau;
  

  /* calculate time averages for the entire trajectory */
  for (tau=0;tau<=niter;tau++)
    {
      if (sum.nr[tau]>0) av.xrtau[tau] += sum.xr[tau]/sum.nr[tau];
      if (sum.nr[tau]>0) av.xr2tau[tau] += sum.xr2[tau]/sum.nr[tau];
      
      if (sum.nl[tau]>0) av.xltau[tau] += sum.xl[tau]/sum.nl[tau];
      if (sum.nl[tau]>0) av.xl2tau[tau] += sum.xl2[tau]/sum.nl[tau];
      
      if ((sum.nr[tau]+sum.nl[tau])>0) 
	av.xtau[tau] += (sum.xr[tau]+sum.xl[tau])/(sum.nl[tau]+sum.nr[tau]);
      if ((sum.nr[tau]+sum.nl[tau])>0)
	av.x2tau[tau] += (sum.xr2[tau]+sum.xl2[tau])/(sum.nl[tau]+sum.nr[tau]);
    }
    
    for (iMT=1; iMT<=mxmts; iMT++) {
        for (ibox=1; ibox<=nbox; ibox++) {
           
            av.nv_cluster[iMT][ibox]/= (float)iav;
        }
        
    }
    
  /* normalize over the ensemble of trajectories */
  for (it=1;it<=niter;it++)
    {
      av.nMT[it]/=(float)iav;
      av.nMTr[it]/=(float)iav;
      av.nMTl[it]/=(float)iav;
      av.nzeros[it]/=(float)iav;

      av.xmin[it] /= (float)iav;
      av.xmax[it] /= (float)iav;
      av.xmin_real[it] /= (float)iav;
      av.xmax_real[it] /= (float)iav;
      av.ymin[it] /= (float)iav;
      av.ymax[it] /= (float)iav;
      av.zmin[it] /= (float)iav;
      av.zmax[it] /= (float)iav;
      av.rxmin[it] /= (float)iav;
      av.rxmax[it] /= (float)iav;

      av.absvmin[it] /= (float)iav;
      av.absvmax[it] /= (float)iav;
      av.vmin[it] /= (float)iav;
      av.vmax[it] /= (float)iav;

      av.xr[it] /= (float)iav;
      av.xr2[it] /= (float)iav;
      av.xl[it] /= (float)iav;
      av.xl2[it] /= (float)iav;
      av.xt[it] /= (float)iav;
      av.xt2[it] /= (float)iav;

      av.rxr[it] /= (float)iav;
      av.rxr2[it] /= (float)iav;
      av.rxl[it] /= (float)iav;
      av.rxl2[it] /= (float)iav;
      av.rxt[it] /= (float)iav;
      av.rxt2[it] /= (float)iav;

      av.velr[it] /= (float)iav;
      av.velr2[it] /= (float)iav;
      av.vell[it] /= (float)iav;
      av.vell2[it] /= (float)iav;
      av.velt[it] /= (float)iav;
      av.velt2[it] /= (float)iav;

      av.neighbors[it] /= (float)iav;
      av.ovlp[it] /= (float)iav;

      av.ryz[it] /= (float)iav;
      av.op[it] /= (float)iav;
        

        av.ovlpsize_perMT[it] /= (float)iav;
        av.ovlpsize_perOV[it] /= (float)iav;

/*       for (ibox=0;ibox<nbox;ibox++) */
/* 	{ */
/* 	  av.nr[it][ibox] /= av.ntot[it];  */
/* 	  av.nl[it][ibox] /= av.ntot[it]; */
/* 	} */
      for (ibox=0;ibox<nbox;ibox++)
	{
	  av.nr[it][ibox] /= (float)iav;
	  av.nl[it][ibox] /= (float)iav;

      
  
      /* for (ibox=0;ibox<nbox;ibox++) */
/* 	{ */
/* 	  if (av.absvel_ntot[it]>0) */
/* 	    { */
/* 	      av.absvel_nr[it][ibox] /= av.absvel_ntot[it]; */
/* 	      av.absvel_nl[it][ibox] /= av.absvel_ntot[it]; */
/* 	    } */
/* 	} */
/*       for (ibox=0;ibox<nbox;ibox++) */
/* 	{ */
/* 	  if (av.vel_ntot[it]>0) */
/* 	    { */
/* 	      av.vel_nr[it][ibox] /= av.vel_ntot[it]; */
/* 	      av.vel_nl[it][ibox] /= av.vel_ntot[it]; */
/* 	    } */
/* 	} */
      av.thick[it][ibox+1] /= (float)iav;
	  av.absvel_nr[it][ibox] /= (float)iav;
	  av.absvel_nl[it][ibox] /= (float)iav;
	  av.vel_nr[it][ibox] /= (float)iav;
	  av.vel_nl[it][ibox] /= (float)iav;
	  av.xvelr[it][ibox] /= (float)iav;
	  av.xvell[it][ibox] /= (float)iav;
	  av.polarcor[it][ibox] /= (float)iav;
      av.rmin[it][ibox] /= (float)iav;
      av.rmax[it][ibox] /= (float)iav;
      av.velcor[it][ibox] /= (float)iav;
      av.tubulin_dens[it][ibox+1]/= (float)iav;
      av.tubulin_r[it][ibox+1]/= (float)iav;
      av.tubulin_l[it][ibox+1]/= (float)iav;
        
     
      av.ovlp_act[it][ibox+1]/= (float)iav;
      av.ovlp_kin[it][ibox+1]/= (float)iav;
      av.ovlp_dyn[it][ibox+1]/= (float)iav;
      av.ovlp_bi[it][ibox+1]/= (float)iav;
  
    }
        
      for (ibox=0;ibox<nbox2;ibox++)
	av.polarcor_ryz[it][ibox] /= (float)iav;
	
    }
    
  for(it1=corr_ti;it1<=corr_tf;it1++)
    for (it2=it1;it2<=corr_tf;it2++)
      {
	idt=it2-it1;
	av.velcort[it1][idt] /= (float)iav;
	av.velcort_pp[it1][idt] /= (float)iav;
	av.velcort_mm[it1][idt] /= (float)iav;
	av.velcort_pm[it1][idt] /= (float)iav;
      }

  for (tau=0;tau<=niter;tau++)
    {
      av.xrtau[tau] /= (float)iav;
      av.xr2tau[tau] /= (float)iav;
      av.xltau[tau] /= (float)iav;
      av.xl2tau[tau] /= (float)iav;
      av.xtau[tau] /= (float)iav;
      av.x2tau[tau] /= (float)iav;
    }
  
}


/* #################################################################### */
/* unnormalize  --  Un Normalize average quantities                     */
/* to recover the original values before the normalization              */



void unnormalize()
{
  long it,it1,it2,idt,ibox,tau, iMT;
  
  
  for (it=1;it<=niter;it++)
    {
      av.nMT[it] *= (float)iav;
      av.nMTr[it] *= (float)iav;
      av.nMTl[it] *= (float)iav;
      av.nzeros[it]*=(float)iav;

      av.xmin[it] *= (float)iav;
      av.xmax[it] *= (float)iav;
      av.xmin_real[it] *= (float)iav;
      av.xmax_real[it] *= (float)iav;
      av.ymin[it] *= (float)iav;
      av.ymax[it] *= (float)iav;
      av.zmin[it] *= (float)iav;
      av.zmax[it] *= (float)iav;

      av.rxmin[it] *= (float)iav;
      av.rxmax[it] *= (float)iav;

      av.absvmin[it] *= (float)iav;
      av.absvmax[it] *= (float)iav;
      av.vmin[it] *= (float)iav;
      av.vmax[it] *= (float)iav;

      av.xr[it] *= (float)iav;
      av.xr2[it] *= (float)iav;
      av.xl[it] *= (float)iav;
      av.xl2[it] *= (float)iav;
      av.xt[it] *= (float)iav;
      av.xt2[it] *= (float)iav;
      
      av.rxr[it] *= (float)iav;
      av.rxr2[it] *= (float)iav;
      av.rxl[it] *= (float)iav;
      av.rxl2[it] *= (float)iav;
      av.rxt[it] *= (float)iav;
      av.rxt2[it] *= (float)iav;

      av.velr[it] *= (float)iav;
      av.velr2[it] *= (float)iav;
      av.vell[it] *= (float)iav;
      av.vell2[it] *= (float)iav;
      av.velt[it] *= (float)iav;
      av.velt2[it] *= (float)iav;

      av.neighbors[it] *= (float)iav;
      av.ovlp[it] *= (float)iav;

      av.ryz[it] *= (float)iav;
      av.op[it] *= (float)iav;
        

      av.ovlpsize_perMT[it] *= (float)iav;
        av.ovlpsize_perOV[it] *= (float)iav;
      
 /*      for (ibox=0;ibox<nbox;ibox++) */
/* 	{ */
/* 	  av.nr[it][ibox] *= av.ntot[it];  */
/* 	  av.nl[it][ibox] *= av.ntot[it]; */
/* 	} */
      for (ibox=0;ibox<nbox;ibox++)
	{
	  av.nr[it][ibox] *= (float)iav;
	  av.nl[it][ibox] *= (float)iav;

     
 /*      for (ibox=0;ibox<nbox;ibox++) */
/* 	{ */
/* 	  if (av.absvel_ntot[it]>0) */
/* 	    { */
/* 	      av.absvel_nr[it][ibox] *= av.absvel_ntot[it]; */
/* 	      av.absvel_nl[it][ibox] *= av.absvel_ntot[it]; */
/* 	    } */
/* 	} */
/*       for (ibox=0;ibox<nbox;ibox++) */
/* 	{ */
/* 	  if (av.vel_ntot[it]>0) */
/* 	    { */
/* 	      av.vel_nr[it][ibox] *= av.vel_ntot[it]; */
/* 	      av.vel_nl[it][ibox] *= av.vel_ntot[it]; */
/* 	    } */
/* 	} */

	  av.absvel_nr[it][ibox] *= (float)iav;
	  av.absvel_nl[it][ibox] *= (float)iav;
        av.vel_nr[it][ibox] *= (float)iav;
        av.vel_nl[it][ibox] *= (float)iav;
        av.xvelr[it][ibox] *= (float)iav;
        av.xvell[it][ibox] *= (float)iav;
        av.polarcor[it][ibox] *= (float)iav;
        av.velcor[it][ibox] *= (float)iav;
        av.rmin[it][ibox] *= (float)iav;
        av.rmax[it][ibox] *= (float)iav;
        av.tubulin_dens[it][ibox+1]*= (float)iav;
        av.tubulin_r[it][ibox+1]*= (float)iav;
        av.tubulin_l[it][ibox+1]*= (float)iav;
        av.ovlp_act[it][ibox+1]*= (float)iav;
        av.ovlp_kin[it][ibox+1]*= (float)iav;
        av.ovlp_dyn[it][ibox+1]*= (float)iav;
        av.ovlp_bi[it][ibox+1]*= (float)iav;
        av.thick[it][ibox+1] *= (float)iav;
	}
	

      for (ibox=0;ibox<nbox2;ibox++)
	av.polarcor_ryz[it][ibox] *= (float)iav;


    }
    
    for (iMT=1; iMT<=mxmts; iMT++) {
        for (ibox=1; ibox<=nbox; ibox++) {
            av.nv_cluster[iMT][ibox]*= (float)iav;
        }
        
    }
  
  for(it1=corr_ti;it1<=corr_tf;it1++)
    for (it2=it1;it2<=corr_tf;it2++)
      {
	idt=it2-it1;
	av.velcort[it1][idt] *= (float)iav;
	av.velcort_pp[it1][idt] *= (float)iav;
	av.velcort_mm[it1][idt] *= (float)iav;
	av.velcort_pm[it1][idt] *= (float)iav;
      }

  for (tau=0;tau<=niter;tau++)
    {
      av.xrtau[tau] *= (float)iav;
      av.xr2tau[tau] *= (float)iav;
      av.xltau[tau] *= (float)iav;
      av.xl2tau[tau] *= (float)iav;
      av.xtau[tau] *= (float)iav;
      av.x2tau[tau] *= (float)iav;
    }
}

      

  
