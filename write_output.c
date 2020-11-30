/* #################################################################### */
/* write_output  --  Collection of output files                         */

/* color 1 black */
/* color 2 red */
/* color 3 green */
/* color 4 blue */
/* color 5 yellow */
/* color 6 brown */
/* color 7 grey */
/* color 8 violet */
/* color 9 cyan */
/* color 10 magenta */
/* color 11 orange */
/* color 12 indigo */
/* color 13 maroon */
/* color 14 turquoise */
/* color 15 green4 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "global_var.h"
#include "nrutil.h"
#define PERIOD (bound.type[0]==PERIODIC && bound.type[1]==PERIODIC)
#define SI12 (iseg-2*iMT+2)   /* seg = {1,2} */


int gn=0; /* graph number */


/* ----------------------------------------------------------------*/
/*This function calculates the av tubulin density above neqsteps*/
void wav_tubulin_dens(){
    int ibox,it;
    char filename[100]="Av_Tubulin_dens.dat";
    float av_min=0.,av_max=0.,boxlength,av_tub=0.;
    FILE *fp;
    
    if (remove (filept_max(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    fprintf(fp,"# %i averaging iterations\n# %i equilibrium iterations\n ",nav, neqsteps);
    fprintf(fp,"# BLA \n");
    fprintf(fp,"#BoxCenter[microM]\tAv number Tubulin dimers\n");
    
    /*First determine av length of ROI*/
    for (it=neqsteps+1; it<=niter; it++) {
        av_min += av.xmin_real[it];
        av_max += av.xmax_real[it];
    }
    av_max/=(niter-neqsteps);
    av_min/=(niter-neqsteps);
    boxlength=((av_max-av_min)/(float)nbox);
    for (ibox=1; ibox<=nbox; ibox++) {
        av_tub=0.;
        for (it=neqsteps+1; it<=niter; it++){
            av_tub += av.tubulin_dens[it][ibox];
        }
        av_tub /= (niter-neqsteps);
        fprintf(fp,"%f\t%f\n", av_min+(ibox*boxlength/2.), av_tub);
    }
    fclose(fp);
}
/* ----------------------------------------------------------------*/
/* ----------------------------------------------------------------*/
/*This function calculates the av tubulin density at given it*/
void wav_tubulin_dens_it(int filenum){
    int ibox,it;
    char filename[100]="Av_Tubulin_dens.dat";
    float av_min=0.,av_max=0.,boxlength,av_tub=0.;
    FILE *fp;
    
    if (remove (filept_max(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    fprintf(fp,"# %i averaging iterations\n# %i equilibrium iterations\n ",nav, neqsteps);
    fprintf(fp,"# BLA \n");
    fprintf(fp,"#BoxCenter[microM]\tAv number Tubulin dimers\n");
    
    /*First determine av length of ROI*/
    av_min = av.xmin_real[filenum];
    av_max = av.xmax_real[filenum];
    
    boxlength=((av_max-av_min)/(float)nbox);
    for (ibox=1; ibox<=nbox; ibox++) {
        av_tub = av.tubulin_dens[filenum][ibox];
        fprintf(fp,"%f\t%f\n", av_min+(ibox*boxlength/2.), av_tub);
    }
    fclose(fp);
}
/* ----------------------------------------------------------------*/

/*This function calculates the av vel of a cluster of size*/
void wav_clustervel(){
    int ibox,it,i;
    char filename[100]="Av_clustervel.dat";
    float absvmin,absvmax,dv, *norm;
    int maxclust;
    FILE *fp;
    
    /*largest n-cluster is perdefault 30 and only different if mxmts<30*/
    if (mxmts >= 30)
        maxclust=30;
    else
        maxclust=mxmts;
    
    norm=vector(1,30);
    
    if (remove (filept_max(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    fprintf(fp,"# %i averaging iterations\n ",nav);
    fprintf(fp,"# BLA \n");
    fprintf(fp,"#vel[micron/s]\t");
    for (i=2; i<=maxclust;i++){
        fprintf(fp,"Avfractsize%i\t",i);
        norm[i]=0.;
    }
    fprintf(fp, "\n");

    dv=maxvel/(float)nbox;

    /*normalize*/
    for (i=2; i<=maxclust; i++) {
        for (ibox=1; ibox<=nbox; ibox++) {
            norm[i]+= av.nv_cluster[i][ibox];
        }
    }
    
    for (ibox=1; ibox<=nbox; ibox++) {
        
        fprintf(fp,"%f\t", (ibox-0.5)*dv);
        for (i=2; i<=maxclust; i++) {
            fprintf(fp,"%f\t", av.nv_cluster[i][ibox]/norm[i]);
        }
        fprintf(fp,"\n");
    }
    
    fclose(fp);
}
/* ----------------------------------------------------------------*/

/* ----------------------------------------------------------------*/



void wav_percprob ()
{
  /* This function calculates the percolation probability of a give chi
  */

  int iMT,jMT,it, ibox;
  float dens=0., bundle_rad=0., av_max=0.,av_min=0.,av_strands=0;
  char filename[100]="percprob_vs_chi.dat";
  FILE *fp;

  if(iav==1){
      if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
      if ((fp=fopen(filept(filename),"a"))==NULL)
          printf("COULDNT OPEN FILE: %s",filename);
  }
  else{
      if ((fp=fopen(filept(filename),"a"))==NULL)
          printf("COULDNT OPEN FILE: %s",filename);
  }
  
  fprintf(fp,"%f\t%f\n", ProbActive, perc_prob/(float)niter);
    
    
  fclose(fp);
}

void wav_percprob_bi ()
{
    /* This function calculates the percolation probability of a give chi
     */
    
    int iMT,jMT,it, ibox;
    float dens=0., bundle_rad=0., av_max=0.,av_min=0.,av_strands=0;
    char filename[100]="percprob_vs_bi.dat";
    FILE *fp;
    
    if(iav==1){
        if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
        if ((fp=fopen(filept(filename),"a"))==NULL)
            printf("COULDNT OPEN FILE: %s",filename);
    }
    else{
        if ((fp=fopen(filept(filename),"a"))==NULL)
            printf("COULDNT OPEN FILE: %s",filename);
    }
    
    fprintf(fp,"%f\t%f\n", ProbBipolar, perc_prob/(float)niter);
    
    
    fclose(fp);
}
/* ----------------------------------------------------------------*/
void wav_stalllength_extf ()
{
    /* This function calculates the stall force of a given fext
     */
    int it;
    float f=0.;
    char filename[100]="Stalllength_vs_fext.dat";
    FILE *fp;
    
    if(iav==1)
        if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);

    if ((fp=fopen(filept(filename),"a"))==NULL)
    printf("COULDNT OPEN FILE: %s",filename);
    
    for (it=neqsteps+1; it<=niter; it++) {
            f += track.rbound[it];
    }
    
    f /= (niter-neqsteps);
    
    fprintf(fp,"%f\t%f\n", ext_fr, f);
    
    
    fclose(fp);
}

/* ----------------------------------------------------------------*/
void wav_stallforce ()
{
    /* This function calculates the stall length depending on Probactive
     */
    int it;
    float f_r=0., f_l=0.;
    char filename[100]="Stallforce_vs_Chi.dat";
    FILE *fp;
    
    if(iav==1)
        if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
    
    if ((fp=fopen(filept(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    for (it=neqsteps+1; it<=niter; it++) {
        /*only count positive force*/
        if (track.rbound[it]>rbound0)
            f_r += track.rbound[it]-rbound0;
            f_l += lbound0-track.lbound[it];
    }
    
    f_r /= (niter-neqsteps)/spring_r;
    f_l /= (niter-neqsteps)/spring_l;
    fprintf(fp,"%f\t%f\n", ProbActive, f_l+f_r);
    
    
    fclose(fp);
}



/* ----------------------------------------------------------------*/

void wav_xvelprof (int filenum)
{
    /* This function calculates an averaged profile <nl(t,x)>,
     <nr(t,x)> and <ntot(t,x)>; where nr and nl are the number
     of MTs pointing with their minus-end to the right and left.
     */
    
    int iMT,jMT;
    char filename[100]="AvXvelProf_";
    char title[100]="Average X-Velocity PROFILE FOR ITERATION No_";
    float xmin,zmin,xmax,zmax,xrange,zrange,dx,x,vel;
    FILE *fp;
    int ibox;
    
    /* These files are denoted graph number 1 == G1 */
    /*gn=1;*/
    
    
    xmin=av.xmin[iter];
    xmax=av.xmax[iter];
    
    zmin=av.vmin[iter];
    zmax=av.vmax[iter];
    dx=(xmax-xmin)/(float)(nbox-1);
    
    /* Open a seperate file for each filenum */
    strcat(filename,(const char *)itoa(filenum));
    strcat(filename,".dat");
    strcat(title,(const char *)itoa(filenum));
    
    if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    /* Write Standard configuration data For XMGR */
    agr_comment(fp,"\n");
    agr_comment(fp,"-----------------------------------------");
    agr_comment(fp,"DENSITY PROFILE");
    
    /* void agr_std(FILE *fp,
     int gn,
     const char *Title,
     const char *Xaxis,
     const char *Yaxis,
     float xmin,
     float ymin,
     float xmax,
     float ymax) */
    agr_std(fp,gn,title,"x/\\f{12}m\\1m","Local averaged Velocity <v(x)>/v0"
            ,xmin,zmin,xmax,zmax);
    
    /* void agr_begset(FILE *fp,
     int gn,
     iset k,
     int lncolor,
     int lnstyle,
     float lnwidth,
     int symbol,
     float symbsz,
     int symbcl)
     */
    
    /*for(ibox=0;ibox<nbox;ibox++)
     printf("av.nr[%d][%d]=%f\n",iter,ibox,av.nr[iter][ibox]);getchar();*/
    
    agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);
    for(ibox=0;ibox<nbox;ibox++)
    {
     
        if (av.nr[iter][ibox]>0)
        {

            x=ibox*dx+av.xmin[iter];
            vel=av.xvelr[iter][ibox]/av.nr[iter][ibox];
            agr_wcoor(fp,x,vel/vel0_d);
        }
    }
    
    agr_begset(fp,gn,1,4,1,2.0,0,0.0,2);
    for(ibox=0;ibox<nbox;ibox++)
    {
     
        if (av.nl[iter][ibox]>0)
        {
            x=ibox*dx+av.xmin[iter];
            vel=av.xvell[iter][ibox]/av.nl[iter][ibox];
            agr_wcoor(fp,x,vel/vel0_d);
        }
    }
    
    agr_begset(fp,gn,2,1,1,2.0,0,0.0,3);
    for(ibox=0;ibox<nbox;ibox++)
    {
        if ((av.nr[iter][ibox]+av.nl[iter][ibox])>0)
        {
            x=ibox*dx+av.xmin[iter];
            vel=(av.xvelr[iter][ibox]+av.xvell[iter][ibox])/
            (av.nl[iter][ibox]+av.nr[iter][ibox]);
            agr_wcoor(fp,x,vel/vel0_d);
        }
    }
    
    fclose(fp);
}

void wav_xvelprofToT ()
{
    /* This function calculates an averaged profile <nl(t,x)>,
     <nr(t,x)> and <ntot(t,x)>; where nr and nl are the number
     of MTs pointing with their minus-end to the right and left.
     */
    
    int iMT,jMT,it,ibox;
    char filename[100]="AvXvelProfTot";
    char title[100]="Average X-Velocity Averaged over all iter";
    float vell[nbox+1],velr[nbox+1];
    FILE *fp;
    
    
    /* These files are denoted graph number 1 == G1 */
    /*gn=1;*/
    
    /*empty bitches*/
    for (ibox=0; ibox<nbox; ibox++) {
        vell[ibox]=0.;
        velr[ibox]=0.;
        for (it=neqsteps+1; it<=niter; it++) {
            if (av.nl[it][ibox]>0) {
               
                vell[ibox]+=av.xvell[it][ibox]/(float)av.nl[it][ibox];
            }
            if (av.nr[it][ibox]>0) {
                velr[ibox]+=av.xvelr[it][ibox]/(float)av.nr[it][ibox];
            }
            if (ibox==50 && vell[ibox]>1000000.){
                printf("it%i xvell=%f nl=%f\n",it,av.xvell[it][ibox],av.nl[it][ibox]);
            }
        }
      
        velr[ibox]/=(float)(niter-neqsteps);
        vell[ibox]/=(float)(niter-neqsteps);

    }
    
    /*fill them*/
    
            
    /* Open a seperate file for each filenum */
    strcat(filename,".dat");
    
    if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
       agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);
    for(ibox=0;ibox<nbox;ibox++)
        agr_wcoor(fp,ibox,velr[ibox]);
    
    agr_begset(fp,gn,1,4,1,2.0,0,0.0,2);
    for(ibox=0;ibox<nbox;ibox++)
        agr_wcoor(fp,ibox,vell[ibox]);
    
    agr_begset(fp,gn,2,1,1,2.0,0,0.0,3);
    for(ibox=0;ibox<nbox;ibox++)
        agr_wcoor(fp,ibox,velr[ibox]+vell[ibox]);
    
    fclose(fp);
}



/* ----------------------------------------------------------------*/

/* ----------------------------------------------------------------*/
/* ----------------------------------------------------------------*/
/*Get tubulin order parameter*/
void wav_tubulin_op(int filenum){
    
    FILE *fp,*fp2;
    float xmin, xmax,dx, xbox, op;
    float av_max=0., av_min=0., bundle_rad=0., dens=0.;
    int ibox,it;
    char filename[100]="AvMTorder_";
    char filename2[100]="AvMTdens_";
    
    /*get average MT density per volume*/
    for (it=1; it<=niter; it++) {
        bundle_rad+= av.ryz[iter];
        for (ibox=1; ibox<=nbox; ibox++) {
            dens+= av.tubulin_dens[it][ibox];
        }
    }
    /*av length of bundle*/
    for (it=neqsteps+1; it<=niter; it++) {
        av_min += av.xmin_real[it];
        av_max += av.xmax_real[it];
    }
    av_max/=(niter-neqsteps);
    av_min/=(niter-neqsteps);
    
    bundle_rad= bundle_rad/niter+exclude;
    dens/= (niter*nbox);
    dens/= bundle_rad*bundle_rad*PI*(av_max-av_min);
    
    
    xmax=av.xmax_real[filenum];
    xmin=av.xmin_real[filenum];
    
    dx= (xmax-xmin)/nbox;
    
    strcat(filename,(const char *)itoa(filenum));
    strcat(filename2,(const char *)itoa(filenum));
    strcat(filename,".dat");
    strcat(filename2,".dat");
    
    if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
    if (remove (filept(filename2))==0) printf("Removed old file:%s\n",filename2);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    if ((fp2=fopen(filept(filename2),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename2);
    
    fprintf(fp,"radius= %f tubulin dens= %f\n",bundle_rad, dens);
    fprintf(fp2,"radius= %f tubulin dens= %f\n",bundle_rad, dens);
    fprintf(fp,"box_it%i [microns]\tOPit%i\n", iter,iter);
    fprintf(fp2,"box_it%i [microns]\tOPit%i\n", iter,iter);
    
    for (ibox=1; ibox<=nbox; ibox++) {
        xbox=(2*ibox-1)*dx/2+xmin;
        if(av.tubulin_r[filenum][ibox]<EPSI && av.tubulin_l[filenum][ibox]<EPSI)
            op=0.;
        else
        /*commented so that the order parameter is equal to ther percentage of + end right*/
            op= (av.tubulin_l[filenum][ibox]/*-av.tubulin_l[filenum][ibox]*/)/(av.tubulin_r[filenum][ibox]+av.tubulin_l[filenum][ibox]);
        fprintf(fp,"%f\t%f\n", xbox, op);
        fprintf(fp2,"%f\t%f\t%f\n", xbox, av.tubulin_r[filenum][ibox],av.tubulin_l[filenum][ibox]);
    }
    
    fclose(fp);
    fclose(fp2);
    
}

/* ----------------------------------------------------------------*/

/* -------------------------------------------------------------*/
/* -------------------------------------------------------------*/

void wav_totabsveldensity ()
{
  /* This function calculates an averaged profile <v_l(t,x)>, 
     <v_r(t,x)> and <v_t(t,x)>; where v_r and v_l are the number 
     of MTs traveling with velocity v, and whose minus-end points 
     to the right and left, respectively. 
  */

  int iMT,jMT,it;
  char filename[100]="Tot_absVelProf";
  char title[100]="TOT AV VELOCITY PROFILE ";
  float absvmin,zmin,absvmax,zmax,dv,v,pv;
  FILE *fp;
  int ibox;
  
  absvmin=0.0;
  absvmax=0.0;
  /*for (it=1;it<=niter;it++)
    {
      absvmin+=av.absvmin[it];
      absvmax+=av.absvmax[it];
    }
  absvmin/=(float)niter;
  absvmax/=(float)niter;
*/
  /* Uncomment to use fixed values for the min/max values of Vel */ 
  /* these are defined in rdata.c                                */
    absvmin=0;
    absvmax=maxvel;

  zmin=0.0;
  zmax=1.0;
  dv=(absvmax-absvmin)/(float)(nbox-1);
  
  /* Add suffice to filename */
  strcat(filename,".dat");
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"TOT VELOCITY PROFILE"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"v/v0","Averaged Velocity Density"
  ,absvmin,0,absvmax,1);
  
  /* void agr_begset(FILE *fp,
     int gn,
     iset k, 
     int lncolor, 
     int lnstyle, 
     float lnwidth, 
     int symbol, 
     float symbsz, 
     int symbcl) 
  */
  
  agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      pv=0.0;
      for (it=neqsteps+1;it<=niter;it++)
          pv+=av.absvel_nr[it][ibox];///dv;
      pv/=(float)(niter-neqsteps);
      v=ibox*dv+absvmin;

      agr_wcoor(fp,v,pv);
    }

 agr_begset(fp,gn,1,4,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      pv=0.0;
        for (it=neqsteps+1;it<=niter;it++)
          pv+=av.absvel_nl[it][ibox];//dv;
      pv/=(float)(niter-neqsteps);
      v=ibox*dv+absvmin;

      agr_wcoor(fp,v,pv);
    }

 agr_begset(fp,gn,2,1,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      pv=0.0;
      for (it=neqsteps+1;it<=niter;it++)
          pv+=(av.absvel_nr[it][ibox]+av.absvel_nl[it][ibox]);//dv;
      pv/=(float)(niter-neqsteps);
      v=ibox*dv+absvmin;

      agr_wcoor(fp,v,pv);
    }
    
  fclose(fp);
}

/* -------------------------------------------------------------*/


void wav_veldensity (int filenum)
{
  /* This function calculates an averaged profile <v_l(t,x)>, 
     <v_r(t,x)> and <v_t(t,x)>; where v_r and v_l are the number 
     of MTs traveling with velocity v, and whose minus-end points 
     to the right and left, respectively. 
  */

  int iMT,jMT;
  char filename[100]="AvVelProf_";
  char title[100]="AV VELOCITY PROFILE FOR ITERATION No_"; 
  float vmin,zmin,vmax,zmax,dv,v,pv;
  FILE *fp;
  int ibox;

  /* These files are denoted graph number 1 == G1 */
  /*gn=1;*/
  
  vmin=av.vmin[iter];
  vmax=av.vmax[iter];

  /* Uncomment to use fixed values for the min/max values of Vel */ 
  /* these are defined in rdata.c                                */
  
  vmin=-maxvel;
  vmax=maxvel;
  

  zmin=0.0;
  zmax=1.0;
  dv=(vmax-vmin)/(float)(nbox-1);
  
  /* Open a seperate file for each filenum */
  strcat(filename,(const char *)itoa(filenum));
  strcat(filename,".dat");
  strcat(title,(const char *)itoa(filenum));
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"VELOCITY PROFILE"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"v/v0","Averaged Velocity Density"
  ,vmin,0,vmax,1);
  
  /* void agr_begset(FILE *fp,
     int gn,
     iset k, 
     int lncolor, 
     int lnstyle, 
     float lnwidth, 
     int symbol, 
     float symbsz, 
     int symbcl) 
  */
  
  agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      v=ibox*dv+av.vmin[iter];
      pv=av.vel_nr[iter][ibox]/dv;
      agr_wcoor(fp,v/vel0_d,pv);
    }
  
  agr_begset(fp,gn,1,4,1,2.0,0,0.0,2);
  for(ibox=0;ibox<nbox;ibox++)
    {
      v=ibox*dv+av.vmin[iter];
      pv=av.vel_nl[iter][ibox]/dv;
      agr_wcoor(fp,v/vel0_d,pv);
    }
  
  agr_begset(fp,gn,2,1,1,2.0,0,0.0,3);  
  for(ibox=0;ibox<nbox;ibox++)
    {
      v=ibox*dv+av.vmin[iter];
      pv=(av.vel_nl[iter][ibox]+av.vel_nr[iter][ibox])/dv;
      agr_wcoor(fp,v/vel0_d,pv);
    }
  
  fclose(fp);
}

/* -------------------------------------------------------------*/

void wav_absveldensity (int filenum)
{
    /* This function calculates an averaged profile <v_l(t,x)>,
     <v_r(t,x)> and <v_t(t,x)>; where v_r and v_l are the number
     of MTs traveling with velocity v, and whose minus-end points
     to the right and left, respectively.
     */
    
    int iMT,jMT;
    char filename[100]="AvabsVelProf_";
    char title[100]="AV absVELOCITY PROFILE FOR ITERATION No_";
    float absvmin,zmin,absvmax,zmax,dv,v,pv;
    FILE *fp;
    int ibox;
    
    /* These files are denoted graph number 1 == G1 */
    /*gn=1;*/
    
    absvmin=av.absvmin[iter];
    absvmax=av.absvmax[iter];
    zmin=0.0;
    zmax=1.0;
    dv=(absvmax-absvmin)/(float)(nbox-1);
    
    /* Open a seperate file for each filenum */
    strcat(filename,(const char *)itoa(filenum));
    strcat(filename,".dat");
    strcat(title,(const char *)itoa(filenum));
    
    if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    /* Write Standard configuration data For XMGR */
    agr_comment(fp,"\n");
    agr_comment(fp,"-----------------------------------------");
    agr_comment(fp,"VELOCITY PROFILE");
    
    /* void agr_std(FILE *fp,
     int gn,
     const char *Title,
     const char *Xaxis,
     const char *Yaxis,
     float xmin,
     float ymin,
     float xmax,
     float ymax) */
    agr_std(fp,gn,title,"v/v0","Averaged Velocity Density"
            ,absvmin,0,absvmax,1);
    
    /* void agr_begset(FILE *fp,
     int gn,
     iset k,
     int lncolor,
     int lnstyle,
     float lnwidth,
     int symbol,
     float symbsz,
     int symbcl)
     */
    
    agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);
    for(ibox=0;ibox<nbox;ibox++)
    {
        v=ibox*dv+av.absvmin[iter];
        pv=av.absvel_nr[iter][ibox]/dv;
        agr_wcoor(fp,v/vel0_d,pv);
    }
    
    agr_begset(fp,gn,1,4,1,2.0,0,0.0,2);
    for(ibox=0;ibox<nbox;ibox++)
    {
        v=ibox*dv+av.absvmin[iter];
        pv=av.absvel_nl[iter][ibox]/dv;
        agr_wcoor(fp,v/vel0_d,pv);
    }
    
    agr_begset(fp,gn,2,1,1,2.0,0,0.0,3);
    for(ibox=0;ibox<nbox;ibox++)
    {
        v=ibox*dv+av.absvmin[iter];
        pv=(av.absvel_nl[iter][ibox]+av.absvel_nr[iter][ibox])/dv;
        agr_wcoor(fp,v/vel0_d,pv);
    }
    
    fclose(fp);
}


/*Calculates a velocity distribution averaged over one simulation*/
void wav_totveldensity ()
{
  /* This function calculates an averaged profile <v_l(t,x)>, 
     <v_r(t,x)> and <v_t(t,x)>; where v_r and v_l are the number 
     of MTs traveling with velocity v, and whose minus-end points 
     to the right and left, respectively. 
  */

  int iMT,jMT,it;
  char filename[100]="Tot_VelProf";
  char title[100]="TOT AV VELOCITY PROFILE";
  float vmin,zmin,vmax,zmax,dv,v,pv;
  FILE *fp;
  int ibox;
  
  vmin=0.0;
  vmax=0.0;
  /*for (it=1;it<=niter;it++)
    {
      vmin+=av.vmin[it];
      vmax+=av.vmax[it];
    }
  vmin/=(float)niter;
  vmax/=(float)niter;
*/
  /* Uncomment to use fixed values for the min/max values of Vel */ 
  /* these are defined in rdata.c                                */
    vmin=-maxvel;
    vmax=maxvel;

  zmin=0.0;
  zmax=1.0;
  dv=(vmax-vmin)/(float)(nbox-1);
  
  /* Add suffice to filename */
  strcat(filename,".dat");
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"TOT VELOCITY PROFILE"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"v/v0","Averaged Velocity Density"
  ,vmin,0,vmax,1);
  
  /* void agr_begset(FILE *fp,
     int gn,
     iset k, 
     int lncolor, 
     int lnstyle, 
     float lnwidth, 
     int symbol, 
     float symbsz, 
     int symbcl) 
  */
  
  agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      pv=0.0;
      for (it=neqsteps+1;it<=niter;it++)
          pv+=av.vel_nr[it][ibox];
      pv/=(float)(niter-neqsteps);
      v=ibox*dv+vmin;
      agr_wcoor(fp,v/vel0_d,pv);
    }

      
   agr_begset(fp,gn,1,4,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      pv=0.0;
      for (it=neqsteps+1;it<=niter;it++)
          pv+=av.vel_nl[it][ibox];
      pv/=(float)(niter-neqsteps);
      v=ibox*dv+vmin;
      agr_wcoor(fp,v/vel0_d,pv);
    }
      
  agr_begset(fp,gn,2,1,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      pv=0.0;
      for (it=neqsteps+1;it<=niter;it++)
          pv+=(av.vel_nr[it][ibox]+av.vel_nl[it][ibox]);
      pv/=(float)(niter-neqsteps);
      v=ibox*dv+vmin;

      agr_wcoor(fp,v/vel0_d,pv);
    }
    
  fclose(fp);
}


/* -----------endvel--------------------------------------------------*/

void wav_end2end()
{
  /* This function calculates the average distance 
     to the most distant MT in the shaft. 
  */
  
  float dens=0., av_strands=0., av_max=0.,av_min=0.;
  int it, ibox;
  char filename[100]="AvEnd2End";
  char title[100]="AV END TO END DISTANCE AS FUNCTION OF TIME"; 
  FILE *fp;
  

    /* Open a seperate file for each filenum */
  /*strcat(filename,(const char *)itoa(filenum));*/
  strcat(filename,".dat");
  /*  strcat(title,(const char *)itoa(filenum));*/
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"AVERAGE END2END");


  /*agr_std(fp,gn,title,"time (sec)","END to End Distance \\f{12}m\\1m"
	  ,0,0,niter*dt,av.xmax_real[niter]-av.xmin_real[niter]);*/
 

  //agr_begset(fp,gn,0,1,1,2.0,0,0.0,1);
    for(it=1;it<=niter;it++){
        
        agr_wcoor(fp,it*dt,av.xmax_real[it]-av.xmin_real[it]);
    }
  fclose(fp);


}

void wav_avovlpsize_perMT()
{
    /*
     */
    
    int it, ibox;
    char filename[100]="AvOvlpSize_perMT";
    FILE *fp;
    
    
    strcat(filename,".dat");
    
    if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    /* Write Standard configuration data For XMGR */
    agr_comment(fp,"\n");
    agr_comment(fp,"-----------------------------------------");
    agr_comment(fp,"AVERAGE OVLP SIZE per MT AS A FUNCTION OF TIME");
    
    for(it=1;it<=niter;it++){
        agr_wcoor(fp,it*dt,av.ovlpsize_perMT[it]);
    }
    fclose(fp);
    
    
}


void wav_avovlpsize_perOV()
{
    /*
     */
    
    int it, ibox;
    char filename[100]="AvOvlpSize_perOV";
    FILE *fp;
    
    
    strcat(filename,".dat");
    
    if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    /* Write Standard configuration data For XMGR */
    agr_comment(fp,"\n");
    agr_comment(fp,"-----------------------------------------");
    agr_comment(fp,"AVERAGE OVLP SIZE per ovlp AS A FUNCTION OF TIME");
    
    for(it=1;it<=niter;it++){
        agr_wcoor(fp,it*dt,av.ovlpsize_perOV[it]);
    }
    fclose(fp);
    
    
}

void wav_avstrands( int filenum)
{
    /* This function calculates the average number of strands in the bundle as a function of x at iteration it
     */

    int ibox;
    char filename[100]="AvStrands_";
    FILE *fp;
    
    /* Open a seperate file for each filenum */
    strcat(filename,(const char *)itoa(filenum));
    strcat(filename,".dat");
    
    if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    /* Write Standard configuration data For XMGR */
    agr_comment(fp,"\n");
    agr_comment(fp,"-----------------------------------------");
    agr_comment(fp,"AVERAGE NUMBER OF STRANDS ITERATION");
    fprintf(fp,"\n%i",filenum);
    agr_comment(fp,"\n");
    
    for(ibox=1;ibox<=nbox;ibox++){
        fprintf(fp,"%f\n",av.thick[filenum][ibox]);
    }
    fclose(fp);
}


void wav_avstrands_tot()
{
    /* This function calculates the average number of strands in the bundle as a function of x
     */
  
    
    int it, ibox;
    char filename[100]="AvStrands_tot";
    FILE *fp;
    float m[nbox];
    

    for (ibox=0; ibox<nbox; ibox++)
        m[ibox]=0.0;

    /*calculate temporal av*/
    for (ibox=1; ibox<=nbox; ibox++) {
        for (it=neqsteps+1; it<=niter; it++) {
            m[ibox-1]+=av.thick[it][ibox];

        }
        m[ibox-1]/=(niter-neqsteps);
    }

    strcat(filename,".dat");
    
    if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);

    /* Write Standard configuration data For XMGR */
    agr_comment(fp,"\n");
    agr_comment(fp,"-----------------------------------------");
    agr_comment(fp,"AVERAGE NUMBER OF STRANDS TOTAL AVERAGE");
    agr_comment(fp,"\n");
    
    for(ibox=0;ibox<nbox;ibox++){
        fprintf(fp,"%f\n",m[ibox]);
    }

    fclose(fp);
}


void wav_x()
{
  /* This function calculates the average distance 
     to the most distant MT in the shaft. 
  */

  int it;
  char filename[100]="Av_x";
  char title[100]="<xm(t)>, <xp(t)>, <xm(t)+xp(t)>"; 
  FILE *fp;

/* These files are denoted graph number 3 -- G3 */ 
  /*gn=3;*/
  
  /* Open a seperate file for each filenum */
  /*strcat(filename,(const char *)itoa(filenum));*/
  strcat(filename,".dat");
  /*  strcat(title,(const char *)itoa(filenum));*/
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"AVERAGE X-POS"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"time (sec)","<x(t)>\\f{12}m\\1m"
	  ,av.xmin[niter],0,niter*dt,av.xmax[niter]);
  
  /* void agr_begset(FILE *fp,
     iset k, 
     int lncolor, 
     int lnstyle, 
     float lnwidth, 
     int symbol, 
     float symbsz, 
     int symbcl) 
  */
  
  agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);
    if (PERIOD)
      for(it=neqsteps;it<=niter;it++)
	agr_wcoor(fp,(it-neqsteps)*dt,av.rxr[it]);
    else
      for(it=1;it<=niter;it++)
	agr_wcoor(fp,it*dt,av.xr[it]);

  agr_begset(fp,gn,1,4,1,2.0,0,0.0,1);
    if (PERIOD)
      for(it=neqsteps;it<=niter;it++)
	agr_wcoor(fp,(it-neqsteps)*dt,av.rxl[it]);
    else
      for(it=1;it<=niter;it++)
	agr_wcoor(fp,it*dt,av.xl[it]);

  agr_begset(fp,gn,2,1,1,2.0,0,0.0,1);
    if(PERIOD)
      for(it=neqsteps;it<=niter;it++)
	agr_wcoor(fp,(it-neqsteps)*dt,av.rxt[it]);
    else
      for(it=1;it<=niter;it++)
	agr_wcoor(fp,it*dt,av.xt[it]);
  
  fclose(fp);


}

void wav_x2()
{
  /* This function calculates the average distance 
     to the most distant MT in the shaft. 
  */

  int it;
  char filename[100]="Av_x2";
  char title[100]="<xm^2(t)>, <xp^2(t)>, <(xm(t)+xp(t))^2>"; 
  FILE *fp;
  
/* These files are denoted graph number 4 -- G4 */ 
  /*gn=4;*/

  /* Open a seperate file for each filenum */
  /*strcat(filename,(const char *)itoa(filenum));*/
  strcat(filename,".dat");
  /*  strcat(title,(const char *)itoa(filenum));*/
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"AVERAGE X2-POS"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"time (sec)","<x2(t)>\\f{12}m\\1m"
	  ,av.xmin[niter],0,niter*dt,av.xmax[niter]);
  
  /* void agr_begset(FILE *fp,
     int gn,
     iset k, 
     int lncolor, 
     int lnstyle, 
     float lnwidth, 
     int symbol, 
     float symbsz, 
     int symbcl) 
  */
  agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);
  /* if (PERIOD) */
    for(it=neqsteps;it<=niter;it++)
      agr_wcoor(fp,(it-neqsteps)*dt,sqrt(av.rxr2[it]));
  /* else */
/*     for(it=1;it<=niter;it++) */
/*       agr_wcoor(fp,it*dt,sqrt(av.xr2[it])); */
  
  agr_begset(fp,gn,1,4,1,2.0,0,0.0,1);
  /* if (PERIOD) */
    for(it=neqsteps;it<=niter;it++)
      agr_wcoor(fp,(it-neqsteps)*dt,sqrt(av.rxl2[it]));
  /* else */
/*     for(it=1;it<=niter;it++) */
/*       agr_wcoor(fp,it*dt,sqrt(av.xl2[it])); */
  
  agr_begset(fp,gn,2,1,1,2.0,0,0.0,1);
  /* if (PERIOD) */
/*     for(it=neqsteps;it<=niter;it++) */
      agr_wcoor(fp,(it-neqsteps)*dt,sqrt(av.rxt2[it]));
  /* else */
/*     for(it=1;it<=niter;it++) */
/*       agr_wcoor(fp,it*dt,sqrt(av.xt2[it])); */
  
  fclose(fp);
}

/* ----------------------------------------------------------*/
void wav_x2tau()
{
  /* This function calculates the average distance 
     to the most distant MT in the shaft. 
  */

  int tau;
  char filename[100]="Av_x2tau";
  char title[100]="<xm^2(t)>, <xp^2(t)>, <(xm(t)+xp(t))^2>"; 
  FILE *fp;
  
/* These files are denoted graph number 4 -- G4 */ 
  /*gn=4;*/

  /* Open a seperate file for each filenum */
  /*strcat(filename,(const char *)itoa(filenum));*/
  strcat(filename,".dat");
  /*  strcat(title,(const char *)itoa(filenum));*/
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"AVERAGE X2-POS"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"time (sec)","<x2(t)>\\f{12}m\\1m"
	  ,av.xmin[niter],0,niter*dt,av.xmax[niter]);
  
  /* void agr_begset(FILE *fp,
     int gn,
     iset k, 
     int lncolor, 
     int lnstyle, 
     float lnwidth, 
     int symbol, 
     float symbsz, 
     int symbcl) 
  */
  agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);
  for(tau=0;tau<=niter;tau++)
    agr_wcoor(fp,tau*dt,sqrt(av.xr2tau[tau]));
  
  agr_begset(fp,gn,1,4,1,2.0,0,0.0,1);
  for(tau=0;tau<=niter;tau++)
    agr_wcoor(fp,tau*dt,sqrt(av.xl2tau[tau]));
  
  agr_begset(fp,gn,2,1,1,2.0,0,0.0,1);
  for(tau=0;tau<=niter;tau++)
    agr_wcoor(fp,tau*dt,sqrt(av.x2tau[tau]));
  
  fclose(fp);
}

/* -----------------------------------------------------------------*/





/* -------------------------------------------------------- */

void wav_V()
{
  /* This function calculates the average MT velocity */

  int it;
  char filename[100]="Av_V";
  char title[100]="<Vl(t)>, <Vr(t)>, <Vl(t)+Vr(t)>"; 
  FILE *fp;

/* These files are denoted graph number 3 -- G3 */ 
  /*gn=3;*/
  
  /* Open a seperate file for each filenum */
  /*strcat(filename,(const char *)itoa(filenum));*/
  strcat(filename,".dat");
  /*  strcat(title,(const char *)itoa(filenum));*/
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"AVERAGE V/V0"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"time (sec)","<V(t)>/V0"
	  ,av.xmin[niter],0,niter*dt,av.xmax[niter]);
  
  /* void agr_begset(FILE *fp,
     iset k, 
     int lncolor, 
     int lnstyle, 
     float lnwidth, 
     int symbol, 
     float symbsz, 
     int symbcl) 
  */
  
  agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);
    for(it=1;it<=niter;it++){
//        printf("it=%i avvelr=%f\n",it,av.velr[it]);
      agr_wcoor(fp,it*dt,av.velr[it]/vel0_d);
    }
  agr_begset(fp,gn,1,4,1,2.0,0,0.0,1);
  for(it=1;it<=niter;it++)
      agr_wcoor(fp,it*dt,av.vell[it]/vel0_d);

  agr_begset(fp,gn,2,1,1,2.0,0,0.0,1);
    for(it=1;it<=niter;it++){
                //printf("it=%i avvelt=%f\n",it,av.velt[it]);
    agr_wcoor(fp,it*dt,av.velt[it]/vel0_d);
    }
  fclose(fp);


}

void wav_V2()
{
  
  /* This function calculates the average MT velocity square */

  int it;
  char filename[100]="Av_V2";
  char title[100]="SQRT <V\\S2\\N(t)>/V0"; 
  FILE *fp;
  
/* These files are denoted graph number 4 -- G4 */ 
  /*gn=4;*/

  /* Open a seperate file for each filenum */
  /*strcat(filename,(const char *)itoa(filenum));*/
  strcat(filename,".dat");
  /*  strcat(title,(const char *)itoa(filenum));*/
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"AVERAGE V2/V0"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"time (sec)","SQRT<(V-avV\\S2\\N(t))^2>/V0"
	  ,av.xmin[niter],0,niter*dt,av.xmax[niter]);
  
  /* void agr_begset(FILE *fp,
     int gn,
     iset k, 
     int lncolor, 
     int lnstyle, 
     float lnwidth, 
     int symbol, 
     float symbsz, 
     int symbcl) 
  */
  
  agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);
  for(it=1;it<=niter;it++)
      agr_wcoor(fp,it*dt,sqrt(av.velr2[it]-av.velr[it]*av.velr[it])/vel0_d);

  agr_begset(fp,gn,1,4,1,2.0,0,0.0,1);
  for(it=1;it<=niter;it++)
      agr_wcoor(fp,it*dt,sqrt(av.vell2[it]-av.vell[it]*av.vell[it])/vel0_d);

  agr_begset(fp,gn,2,1,1,2.0,0,0.0,1);
  for(it=1;it<=niter;it++)
    agr_wcoor(fp,it*dt,sqrt(av.velt2[it]-av.velt[it]*av.velt[it])/vel0_d);
  fclose(fp);


}


/* -----------------------------------------------------------*/

void wav_neighbors()
{
  /* This function prints the average number of neighbors interacting  
     with each MT; alternatively the average coordination number. 
  */

  int it;
  char filename[100]="Av_neighbors";
  char title[100]="Average Coordination Number; <z(t)>"; 
  FILE *fp;
  
/* These files are denoted graph number 4 -- G4 */ 
  /*gn=4;*/

  /* Open a seperate file for each filenum */
  /*strcat(filename,(const char *)itoa(filenum));*/
  strcat(filename,".dat");
  /*  strcat(title,(const char *)itoa(filenum));*/
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"AVERAGE Coordination Number"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"time (sec)","<z(t)>"
	  ,0,0,niter*dt,10);
  
  /* void agr_begset(FILE *fp,
     int gn,
     iset k, 
     int lncolor, 
     int lnstyle, 
     float lnwidth, 
     int symbol, 
     float symbsz, 
     int symbcl) 
  */
  
  agr_begset(fp,gn,0,1,1,2.0,0,0.0,1);
  for(it=1;it<=niter;it++)
      agr_wcoor(fp,it*dt,av.neighbors[it]);
  
  fclose(fp);
}

/* -----------------------------------------------------------*/

void wav_ovlp()
{
  /* This function prints the average number of neighbors interacting  
     with each MT; alternatively the average coordination number. 
  */

  int it;
  char filename[100]="Av_ovlp";
  char title[100]="Average ovlp distance; <Lij(t)>"; 
  FILE *fp;
  
/* These files are denoted graph number 4 -- G4 */ 
  /*gn=4;*/

  /* Open a seperate file for each filenum */
  /*strcat(filename,(const char *)itoa(filenum));*/
  strcat(filename,".dat");
  /*  strcat(title,(const char *)itoa(filenum));*/
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"AVERAGE Coordination Number"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"time (sec)","<Lij(t)>"
	  ,0,0,niter*dt,10);
  
  /* void agr_begset(FILE *fp,
     int gn,
     iset k, 
     int lncolor, 
     int lnstyle, 
     float lnwidth, 
     int symbol, 
     float symbsz, 
     int symbcl) 
  */
  
  agr_begset(fp,gn,0,1,1,2.0,0,0.0,1);
  for(it=1;it<=niter;it++)
      agr_wcoor(fp,it*dt,av.ovlp[it]);
  
  fclose(fp);
}

/* -------------------------------------------------------------*/

void wav_nMT()
{
  /* This function calculates the average number of MTs; in time. 
  */

  
  int it;
  char filename[100]="Av_nMT";
  char title[100]="Average number of MTs"; 
  FILE *fp;
  
  
  /* These files are denoted graph number 2 -- G2 */ 
  /*gn=2;*/

  /* Open a seperate file for each filenum */
  /*strcat(filename,(const char *)itoa(filenum));*/
  strcat(filename,".dat");
  /*  strcat(title,(const char *)itoa(filenum));*/
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"AVERAGE nMTs"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"time (sec)","<N\\sMT\\N>"
	  ,0,0,niter*dt,av.nMT[niter]);
  
  /* void agr_begset(FILE *fp,
     int gn,
     iset k, 
     int lncolor, 
     int lnstyle, 
     float lnwidth, 
     int symbol, 
     float symbsz, 
     int symbcl) 
  */
  

  agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);
  for(it=1;it<=niter;it++)
      agr_wcoor(fp,it*dt,av.nMTr[it]);

  agr_begset(fp,gn,1,4,1,2.0,0,0.0,1);
  for(it=1;it<=niter;it++)
      agr_wcoor(fp,it*dt,av.nMTl[it]);

  agr_begset(fp,gn,2,1,1,2.0,0,0.0,1);
  for(it=1;it<=niter;it++)
      agr_wcoor(fp,it*dt,av.nMT[it]);
  
  fclose(fp);


}


/* -------------------------------------------------------------*/


/* ----------------------------------------------------*/


/* -------------------------------------------------------------*/
/* -------------------------------------------------------------*/



void wav_xtraject()
{
  /* This function calculates single MT trajectories
  */

  char filename[100]="MT_Traject";
  char title[100]="Trajectories";
  float tmin,xmin,tmax,xmax,xcm,time;
  int i,it,jMT,lc,n,num,t0;
  FILE *fp;

  /* These files are denoted graph number 1 == G1 */
  /*gn=1;*/
  
  tmin=0.0;
  tmax=niter*dt;
  xmin=-1000;
  xmax=1000;

  strcat(filename,".dat");
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL)
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"SINGLE MT TRAJECTORIES");
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title,
     const char *Xaxis,
     const char *Yaxis,
     float xmin,
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"t (sec)","x (\\f{12}m\\1m)",
	  xmin,tmin,xmax,tmax);
  
  /* void agr_begset(FILE *fp,
     int gn,
     iset k,
     int lncolor,
     int lnstyle,
     float lnwidth,
     int symbol,
     float symbsz,
     int symbcl)
  */

  /* number of trajectories to draw */
  if (track.index[niter][MT.number]<50) 
    n=track.index[niter][MT.number];
  else
    n=50;

  num=0;
  i=1;
  while (num<n)
    {
      agr_begset(fp,gn,i,i%15+1,1,2.0,0,0.0,1);
      fprintf(fp,"# Trajectory for iMT No.%d  -  \n",i);
      if (getjMT(i,1)!=0){
          num++;
          for (it=1;it<=niter;it++){
              jMT=getjMT(i,it);
              if(jMT!=0){
                  xcm=track.xcm[it][jMT];
                  t0=track.it0[it][jMT]; /* t0 is the iteration the MT entered
                            the system; we start measuring the 
                            RMSD after neqsteps. */
                  if (t0>neqsteps)
                      time=dt*(it-t0);
                  else
                      time=dt*(it-neqsteps);
                  
                  agr_wcoor(fp,time,xcm);
              }
           }
	  }
      i++;
    }
  fclose(fp);
}

///////////////////////////////////////////////
// get output according to mathematica input///
///////////////////////////////////////////////
void wav_traject_2Dmathematica()
{
    /* This gives the 2D mathematica output
     */
    
    char filename[100];
    float tmin,xmin,tmax,xmax,xcm, ycm,time, length;
    int i,it,jseg,t0,n,q;
    FILE *fp;
    int count[niter];
    
    emptyint(count,0, niter);
    
    /* These files are denoted graph number 1 == G1 */
    /*gn=1;*/
    
    tmin=0.0;
    tmax=niter*dt;
    xmin=-1000;
    xmax=1000;
    
    /*Print parameter file (only niter for now)*/
    sprintf(filename,"parameter.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
    printf("COULDNT OPEN FILE: %s",filename);
    /*niter*/
    fprintf(fp,"%i\n", niter);
    /*max plotrange*/
    fprintf(fp,"%f\n", MT.xmax_real+MT.length[MT.imax_real]/2);
    fprintf(fp, "%f\n", MT.xmin_real-MT.length[MT.imin_real]/2);
    fclose(fp);
    
    
    //Print the color file
    sprintf(filename,"colors.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    for (it=1; it<=niter; it++) {
        for (jseg=1; jseg<=track.nMT[it]*2; jseg++) {
            if (track.length[it][jseg]>0.001) {
                if (track.direct[it][jseg] < -0.1) {
                    
                    
                    if (it-track.it0[it][(int)((float)jseg/2.+0.5)] < 100 && track.it0[it][(int)((float)jseg/2.+0.5)]>1){
                        fprintf(fp,"Purple\t");
                    }
                    else{
                        fprintf(fp,"Blue\t");
                    }
                    //printf("num:%i %f\n",jseg, track.direct[it][jseg]);
                }
                else if (track.direct[it][jseg] > 0.1) {
                    
                    if (it-track.it0[it][(int)((float)jseg/2.+0.5)] < 100 && track.it0[it][(int)((float)jseg/2.+0.5)]>1){
                        fprintf(fp,"Pink\t");
                        
                    }
                    else{
                        fprintf(fp,"Red\t");
                        
                        //printf("num:%i %f\n", jseg, track.direct[it][jseg]);
                    }
                }
                
        
            
                else if (track.direct[it][jseg] < -1.1 || track.direct[it][jseg] > 1.1 || (track.direct[it][jseg] > -0.1 && track.direct[it][jseg] < 0.1) )
                    printf("ERROR tracknmt=%i",track.nMT[it]);
            }
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    

    
    
    //Get total number of Filaments to simulate at each it step
    sprintf(filename,"MTnumbers.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
            printf("COULDNT OPEN FILE: %s",filename);
    for (it=1; it<=niter; it++) {
        i=1;
        for (q=1; q <= 2*track.nMT[it]; q++) {
            if (track.length[it][q]>0.001) {
                fprintf(fp, "%i\t",i++);
            }
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    
    //Print out arrow track data////
    sprintf(filename,"a_num.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    for (it=1; it<=niter; it++) {
        for (i=1; i<=track.arrow_number[it]; i++) {
            fprintf(fp, "%i\t", i);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    
    //Print out arrow color data////
    sprintf(filename,"a_col.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    for (it=1; it<=niter; it++) {
        for (i=1; i<=track.arrow_number[it]; i++) {
            /*Kinesin will be blue*/
            if (track.cross_type[it][i]==-1)
                fprintf(fp, "Blue\t" );
            /*and dynein red*/
            else if (track.cross_type[it][i]==1)
                fprintf(fp, "Red\t" );
            else{
                printf("ovlp type not specified!\n >>mathematica arrows\n");
                exit(0);
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    
    // PRINT ARROW COORDINATES
    sprintf(filename,"arrows.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    for (it=1; it<=niter; it++) {
        
        for (i=1; i<=(int)((float)track.arrow_number[it]*4); i++){
            
            fprintf(fp, "%f\t", track.mathematica_arrows[it][i]);
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);
    
    //PRINT BOUNDARY VISUALS
    sprintf(filename,"bound.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    for (it=1; it<=niter; it++) {
        fprintf(fp, "%f\t%f\t%f\t%f\t", track.lbound[it] -2. , -30., track.lbound[it], 30. );
        fprintf(fp, "%f\t%f\t%f\t%f\n", track.rbound[it], -30., track.rbound[it] + 2. , 30. );
    }
    
    fclose(fp);
    
    //PRINT MT TRACKS
    
    sprintf(filename,"MT_Track.dat");
        if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    
    for (it=1; it<=niter; it++) {
        for (jseg=1; jseg<=2*track.nMT[it]; jseg++) {
            if (track.length[it][jseg]>0.001) {
                xcm=track.xcm_vis[it][jseg];
                ycm=(track.ycm_vis[it][jseg]+track.zcm_vis[it][jseg])/sqrt(2.);//just for dramatic purposes (rotates system around 45 degrees so that we can see more MT in 2D)
                length=track.length[it][jseg];
                time=dt*(it);
                fprintf(fp,"%f\t%f\t%f\t%f\t",time,xcm,100*ycm,length);
            }
            
        }
        fprintf(fp,"\n");
        
    }

    
    fclose(fp);
    
}

void wav_traject_3Dmathematica()
{
    /*This Function gives the mathematica visuals 3D output*/
    
    char filename[100];
    float tmin,xmin,tmax,xmax,xcm, ycm, zcm,time, length;
    int i,it,jseg,t0,n,q;
    FILE *fp;
    int count[niter];
    
    emptyint(count,0, niter);
    
    /* These files are denoted graph number 1 == G1 */
    /*gn=1;*/
    
    tmin=0.0;
    tmax=niter*dt;
    xmin=-1000;
    xmax=1000;
    
    /*Print parameter file (only niter for now)*/
    sprintf(filename,"parameter.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    /*niter*/
    fprintf(fp,"%i\n", niter);
    /*max plotrange*/
    fprintf(fp,"%f\n", MT.xmax_real+MT.length[MT.imax_real]/2);
    fprintf(fp,"%f\n", MT.xmin_real-MT.length[MT.imin_real]/2);
    
    fclose(fp);
    
    
    
    
    sprintf(filename,"colors.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    for (it=1; it<=niter; it++) {
        
        for (jseg=1; jseg<=2*track.nMT[it]; jseg++) {
            
            if (track.length[it][jseg]>0.001) {

                if (track.direct[it][jseg] < -0.1) {

                    
                    if (it-track.it0[it][(int)((float)jseg/2.+0.5)] < 100 && track.it0[it][(int)((float)jseg/2.+0.5)]>1){
                        fprintf(fp,"Purple\t");
                    }
                    else{
                        fprintf(fp,"Blue\t");
                    }
                    //printf("num:%i %f\n",jseg, track.direct[it][jseg]);
                }
                else if (track.direct[it][jseg] > 0.1) {

                    if (it-track.it0[it][(int)((float)jseg/2.+0.5)] < 100 && track.it0[it][(int)((float)jseg/2.+0.5)]>1){
                        fprintf(fp,"Pink\t");
                        
                    }
                    else{
                        fprintf(fp,"Red\t");
                        
                        //printf("num:%i %f\n", jseg, track.direct[it][jseg]);
                    }
                    
                }
                
                if (track.direct[it][jseg] < -1.1  ) printf("ERROR");
                if (track.direct[it][jseg] > 1.1  ) printf("ERROR");
                if (track.direct[it][jseg] > -0.1 && track.direct[it][jseg] < 0.1 ) printf("ERROR");
            }
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    
    
    //Print out arrow color data////
    sprintf(filename,"a_col.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    for (it=1; it<=niter; it++) {
        for (i=1; i<=track.arrow_number[it]; i++) {
            /*Kinesin will be blue*/
            if (track.cross_type[it][i]==-1)
                fprintf(fp, "Blue\t" );
            /*and dynein red*/
            if (track.cross_type[it][i]==1)
                fprintf(fp, "Red\t" );
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    
    
    //Get total number of Filaments to simulate at each it step
    sprintf(filename,"MTnumbers.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    for (it=1; it<=niter; it++) {
        i=1;
        for ( q=1; q <= 2*track.nMT[it]; q++) {
            if (track.length[it][q]>0.001) {
                fprintf(fp, "%i\t",i++);
            }
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    
    //Print out arrow track data////
    sprintf(filename,"a_num.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    for (it=1; it<=niter; it++) {
        for (i=1; i<=track.arrow_number[it]; i++) {
            fprintf(fp, "%i\t", i);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    
    // PRINT ARROW COORDINATES
    sprintf(filename,"arrows3D.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    for (it=1; it<=niter; it++) {
        
        for (i=1; i<=(int)((float)track.arrow_number[it]*6); i++){
            
            fprintf(fp, "%f\t", track.mathematica_arrows[it][i]);

        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);
    
    //PRINT BOUNDARY VISUALS
    sprintf(filename,"bound.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    for (it=1; it<=niter; it++) {
        fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%f\t", track.lbound[it] -2. , 0., 0., track.lbound[it], 0., 0. );
        fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%f\n", track.rbound[it], 0., 0., track.rbound[it] + 2. , 0., 0. );
    }
    
    fclose(fp);
    
    //PRINT BUNDLE LENGTH
    sprintf(filename,"bundle_length.txt");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    for (it=1; it<=niter; it++) {
        fprintf(fp, "%f\t%f\n" , dt*it, track.rbound[it]-track.lbound[it] );
    }
    
     fclose(fp);
    
    
    //PRINT MT TRACKS
    
    sprintf(filename,"MT_Track3D.dat");
    if (remove (filept_mathematica(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept_mathematica(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    
    for (it=1; it<=niter; it++) {
        for (jseg=1; jseg<=2*track.nMT[it]; jseg++) {
            if (track.length[it][jseg]>0.001) {
                xcm=track.xcm_vis[it][jseg];
                ycm=track.ycm_vis[it][jseg];
                zcm=track.zcm_vis[it][jseg];
                length=track.length[it][jseg];
                time=dt*(it);
                fprintf(fp,"%f\t%f\t%f\t%f\t%f\t",time,xcm,100*ycm,100*zcm, length);
            }
            
        }
        fprintf(fp,"\n");
        
    }
    
    
    fclose(fp);
    
}

//This function is needed to plot the maximum force exerted by a bundle in the current av it
// length! SEE INIT.c
void write_force_fstall() {
    
    char filename[100];
    float bound_av=0, f=0;
    int it;
    FILE *fp;
    
    sprintf(filename,"length_time_extf=%f.dat",ext_fr);
    if (iav==1) {
        if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
        if ((fp=fopen(filept(filename),"w"))==NULL)
            printf("COULDNT OPEN FILE: %s",filename);
        fprintf(fp,"time [s]\tlength [microns]\n");
    }
    if ((fp=fopen(filept(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    for ( it=neqsteps+1; it<=niter; it++)
        bound_av+= track.rbound[it];
    bound_av/= niter-neqsteps;
    
    if (bound_av>rbound0)
        f= (bound_av-rbound0)*spring_r;
    else
        f=0.;
    
    fprintf(fp, "%f\t%f\n", fstall_k, f);
    
    
    fclose(fp);
}

//This function is needed to plot force against polarity
/*For this we average all bundle lengths after neqsteps for a given simulation*/
//  SEE INIT.c
void write_force_polarity() {
    
    char filename[100];
    int it;
    float rbound_av=0, lbound_av,f;
    FILE *fp;
    
    
    sprintf(filename,"polarity_exertedforce.dat");
    if (iav==1){
        if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
        if ((fp=fopen(filept(filename),"w"))==NULL)
            printf("COULDNT OPEN FILE: %s",filename);
        fprintf(fp,"#polarity ratio \t exerted force [pN]\n");
        
    }
    else{
        if ((fp=fopen(filept(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    }
    
    for (it=neqsteps+1; it<=niter; it++) {
        rbound_av+=track.rbound[it];
        lbound_av+=track.lbound[it];
    }
    rbound_av /= niter-neqsteps;
    lbound_av /= niter-neqsteps;

    f= spring_r*(rbound_av-rbound0)+spring_l*(lbound0-lbound_av);
	

    printf("%f %f\n",PolarityRatio0,f);
    fprintf(fp, "%f\t%f\n", PolarityRatio0, f );
    
    fclose(fp);
}

//This function is needed to plot force against motor fraction/*For this we average all bundle lengths after neqsteps for a given simulation*/
//  SEE INIT.c
void write_force_motorfraction() {
    
    char filename[100];
    int it;
    float rbound_av=0, lbound_av=0, f_r, f_l;
    FILE *fp;
    
    
    sprintf(filename,"motfrac_exertedforce.dat");
    if (iav==1){
        if (ProbBipolar>0.00001 || ProbKinesin<1.){
            printf("motfrac mode only working if ProbBipolar equal to zero and ProbKinesin equal to 1 initially\n");
            exit(0);
        }
        if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
        if ((fp=fopen(filept(filename),"w"))==NULL)
            printf("COULDNT OPEN FILE: %s",filename);
        fprintf(fp,"#bipolar frac \t exerted force [pN]\n");

    }
    else{
        if ((fp=fopen(filept(filename),"a"))==NULL)
            printf("COULDNT OPEN FILE: %s",filename);
    }
    
    for (it=neqsteps+1; it<=niter; it++) {
        rbound_av+=track.rbound[it];
        lbound_av+=track.lbound[it];
    }
    rbound_av /= niter-neqsteps;
    lbound_av /= niter-neqsteps;
    f_r= spring_r*(rbound_av-rbound0);
    f_l= spring_l*(lbound0-lbound_av);

    fprintf(fp, "%f\t%f\n", ProbBipolar, f_r+f_l );
    
    fclose(fp);
}

void write_force_spring_r() {
    
    char filename[100];
    int it;
    float dens=0., bundle_rad=0., av_max=0.,av_min=0.;
    float bound_av=0, f;
    
    FILE *fp;
    
    /*get average bundle thickness from the first few its
     because at the beginning thickness should be the same independent
     of k*/
    for (it=1; it<=10; it++) {
        bundle_rad+= av.ryz[iter];
    }
    bundle_rad= bundle_rad/10.+exclude;



    
    sprintf(filename,"k_exertedforce.dat");
    if (iav==1){
        if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
        if ((fp=fopen(filept(filename),"w"))==NULL)
            printf("COULDNT OPEN FILE: %s",filename);
        fprintf(fp,"#k \t exerted force [pN]\t initial radius\n");
        fclose(fp);
    }

    if ((fp=fopen(filept(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    for (it=neqsteps+1; it<=niter; it++) {
        bound_av+=track.rbound[it];
    }
    bound_av /= niter-neqsteps;
	if (bound_av> rbound0 )
		f= spring_r*(bound_av-rbound0);
	else
		f=0.;
    fprintf(fp, "%f\t%f\t%f\n", spring_r, f, bundle_rad );
    
    fclose(fp);
}
