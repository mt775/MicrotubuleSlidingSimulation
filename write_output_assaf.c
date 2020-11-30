/* #################################################################### */
/* write_output  --  Collection of output files                                */

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
#include "global_var.h"
#define PERIOD (bound.type[0]==PERIODIC && bound.type[1]==PERIODIC)
#define SI12 (iseg-2*iMT+2)   /* seg = {1,2} */


int gn=0; /* graph number */

/* ------------------------------------------------------------ */
/* The two functions below generate a map of the MT array in the 
   XZ and YZ planes.                                            */
/* ------------------------------------------------------------ */
void wmap_xz(int filenum)
{
  int iMT,jMT,imin,imax,icolor,lc,gn,fc,fpt,arrowtype;
  int iseg,jseg,iseg1,iseg2;
  char filename[100]="MapXZ_";
  char title[100]="MT ARRAY FOR ITERATION No_"; 
  float xmin,zmin,xmax,zmax;
  float lw,sum;
  double ov;
  FILE *fmap;
  struct point pmhead, pmleg;

  /* These files are writen as Graph number 0 */
  gn=0;

  /* Obtain the range of the graph */
  if (bound.type[0]==PERIODIC || bound.type[0]==ABSORB 
      || bound.type[0]==POPUP)
    xmin=lbound0;
  else
  xmin=MT.xmin-0.5*MT.length[MT.imin];
  if (bound.type[1]==PERIODIC || bound.type[1]==ABSORB 
      || bound.type[1]==POPUP)
    xmax=rbound0;
  else
    xmax=MT.xmax+0.5*MT.length[MT.imax];
  zmin=-0.5;
  zmax=0.5;

  /* Open a seperate file for each filenum */
  strcat(filename,(const char *)itoa(filenum));
  strcat(filename,".dat");
  strcat(title,(const char *)itoa(filenum));

  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fmap=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);

 /* Write Standard configuration data For XMGR */
  agr_comment(fmap,"\n");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"MAP OF THE MT ARRAY"); 

  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis) */
  agr_std(fmap,gn,title,"x/\\f{12}m\\1m","z/\\f{12}m\\1m",xmin,zmin,xmax,zmax);
  
  /* Generate a plot of the MTs */
  /*gn=0;*/ /* graph number */
  lc=1; /* line color */
  lw=2; /* line width */
  arrowtype=2; /* at end */

  /* agr_comment(fmap,"MT ARROWS"); */
/*   agr_comment(fmap,"---------");     */
/*   for (iseg=1;iseg<=2*MT.number;iseg++) */
/*     { */
/*       if (seg.length[iseg]>0) */
/* 	if (iseg%2!=0) */
/* 	  agr_arrowline(fmap,gn,lw,lc,seg.pend[iseg].x, seg.pend[iseg].z,  */
/* 			seg.mend[iseg].x, seg.mend[iseg].z, arrowtype); */
/* 	else */
/* 	  agr_arrowline(fmap,gn,lw,lc,seg.pend[iseg].x,  */
/* 			seg.pend[iseg].z-0.1*exclude,  */
/* 			seg.mend[iseg].x, seg.mend[iseg].z-0.1*exclude, */
/* 			arrowtype); */
/*     } */

  /* Generate a filled box around all MT segments */
  agr_comment(fmap,"FILLED BOXES");
  agr_comment(fmap,"------------");    
  lw=1;  
  fpt=1;
  for (iseg=1;iseg<=2*MT.number;iseg++)
    {
      iMT=segindex(iseg);
      if(seg.length[iseg]>0)
	{
	  if (seg.direct[iseg].x>0)
	    {
	      if (SI12==1)
		{
		  fc=10; /*fc=2;*/
		  lc=10; /*lc=2;*/
		}
	      else
		{
		  if (seg.length[iseg]>0) fc=11;
		  if (seg.length[iseg]>0) lc=11;
		}
	    }
	  else
	    {
	      if (SI12==1)
		{
		  fc=14;/*fc=4;*/
		  lc=14;/*lc=4;*/
		}
	      else
		{
		  if (seg.length[iseg]>0) fc=9;
		  if (seg.length[iseg]>0) lc=9;
		}
	    }	
	  if (iseg%2!=0)
	    agr_box(fmap,seg.cm[iseg].x,seg.cm[iseg].z,seg.length[iseg],
		    2*exclude,lw,lc,fpt,fc);
	  else
	   agr_box(fmap,seg.cm[iseg].x,seg.cm[iseg].z-0.1*exclude,
		   seg.length[iseg],
		    2*exclude,lw,lc,fpt,fc); 
	}  
    }
  /* write the index number of each MT above all segments */ 
  agr_comment(fmap,"INDEX NUMBERS");
  agr_comment(fmap,"-------------");
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      iseg1=2*iMT-1;
      iseg2=2*iMT;
      if (seg.length[iseg1]>0)
	agr_string(fmap,gn,seg.cm[iseg1].x,seg.cm[iseg1].z,MT.index[iMT]);
      if (seg.length[iseg2]>0)
	agr_string(fmap,gn,seg.cm[iseg2].x,seg.cm[iseg2].z,MT.index[iMT]);
    }
  /* Generate a plot of motors at overlap regions     */
  /* Do not plot the motors in the overlap between two 
     segments of the same MT.                         */
  /*gn=0;*/
  agr_comment(fmap,"MOTORS");
  agr_comment(fmap,"------");
  lw=2.0;
  for (iseg=1;iseg<=2*MT.number-1;iseg++)
    {  
      iMT=segindex(iseg);
      for (jseg=iseg+1;jseg<=2*MT.number;jseg++)
	{
	  jMT=segindex(jseg);
	  /* segindex(iseg)=iMT */
	  if (soverlap(iseg,jseg)>0 && segindex(iseg)!=segindex(jseg) 
	      && ovlp.type[iseg][jseg]!=ZERO)
	    {
	      if (ovlp.type[iseg][jseg]==BIPOLAR)
		{
		  lc=11; /* orange */
		  lw=3.0;
		  arrowtype=3;
		  motorvec(iseg,jseg,&pmhead,&pmleg);
		  agr_arrowline(fmap,gn,lw,lc,pmhead.x,pmhead.z,
				pmleg.x,pmleg.z,arrowtype); 
		}
	      else if (ovlp.type[iseg][jseg]==TWOMOTORS)
		{
		  lc=3; /* green */
		  lw=2.0;
		  arrowtype=2;
		  if (seg.direct[iseg].x*seg.direct[jseg].x>0)
		    {
		      /* parallel pair */
		      motorvec(iseg,jseg,&pmhead,&pmleg);
		      agr_arrowline(fmap,gn,lw,lc,pmhead.x,pmhead.z,
				pmleg.x,pmleg.z,arrowtype);
		      agr_arrowline(fmap,gn,lw,lc,pmhead.x,pmleg.z,
				pmleg.x,pmhead.z,arrowtype);
		    }
		  else
		    {
		      float dx=0.1;
		      /* anti-parallel pair */
		      motorvec(iseg,jseg,&pmhead,&pmleg);
		      agr_arrowline(fmap,gn,lw,lc,pmhead.x+dx,pmhead.z,
				pmleg.x+dx,pmleg.z,arrowtype);
		      agr_arrowline(fmap,gn,lw,lc,pmhead.x-dx,pmleg.z,
				pmleg.x-dx,pmhead.z,arrowtype);
		    }
		}
	      else if (ovlp.type[iseg][jseg]==BUNDLING 
		       /*&& 
		       seg.direct[iseg].x*seg.direct[jseg].x<0 
		       && fabs(MT.pend[iMT].x-MT.pend[jMT].x)<3.0*/)
		{
		  lc=10; /* magenta */
		  lw=4;
		  arrowtype=0;
		  motorvec(iseg,jseg,&pmhead,&pmleg);
		  agr_arrowline(fmap,gn,lw,lc,pmhead.x,pmhead.z,
				pmleg.x,pmleg.z,arrowtype);
		}
	      else /* Unipolar motor, direction determined in motorvec */ 
		{
		  lc=4;
		  lw=3;
		  arrowtype=2;
		  motorvec(iseg,jseg,&pmhead,&pmleg);
		  agr_arrowline(fmap,gn,lw,lc,pmhead.x,pmhead.z,
				pmleg.x,pmleg.z,arrowtype);
		}
	    }
	}
    }
  /* Print the overlap matrix at the bottom of the AGR file */
  
  agr_comment(fmap,"");
  agr_comment(fmap,"OVERLAP TYPE MATRIX");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"0. = no ovlp  1 = leg up     2 = leg down"); 
  agr_comment(fmap,"3. = bipolar  4 = twomotors  5 = bundling");
  agr_comment(fmap,"-----------------------------------------");
  for (iseg=1;iseg<=2*MT.number-1;iseg++)
    {
      iMT=segindex(iseg);
      if (seg.length[iseg]>0)
	{
	  fprintf(fmap,"# iMT.iseg %d.%d  -  ",MT.index[iMT],
		  iseg-2*iMT+2);
	  for (jseg=iseg+1;jseg<=2*MT.number;jseg++)
	    if (seg.length[jseg]>0)
	      fprintf(fmap,"%d ",ovlp.type[iseg][jseg]);
	  fprintf(fmap,"\n");
	}
    }
  /* Print the overlap matrix at the bottom of the AGR file */
  agr_comment(fmap,"");
  agr_comment(fmap,"OVERLAP MATRIX");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"-----------------------------------------");
  for (iseg=1;iseg<=(2*MT.number-1);iseg++)
    {
      if(seg.length[iseg]>0)
	{
	  iMT=segindex(iseg);
	  fprintf(fmap,"# iMT.iseg %d.%d  -  ",MT.index[iMT],
		  iseg-2*iMT+2);
	  for (jseg=iseg+1;jseg<=(2*MT.number);jseg++)
	    {
	      if(seg.length[jseg]>0)
		{
		  fprintf(fmap,"%6.4f ",soverlap(iseg,jseg));
		}
	    }
	  fprintf(fmap,"\n");
	}
    }
  /* List MT velocities at the bottom of the AGR file */
  agr_comment(fmap,"");
  agr_comment(fmap,"MT Velocities");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"-----------------------------------------");
  sum=0.0;
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      fprintf(fmap,"# v[%d] = %f\n",MT.index[iMT],MT.vel[iMT].x);
      sum+=MT.vel[iMT].x;
    }
  fprintf(fmap,"# TOTAL Velocity: %f\n",sum);

/* List MT x-coordinate at the bottom of the AGR file */
  agr_comment(fmap,"");
  agr_comment(fmap,"MT x-Coordinates");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"-----------------------------------------");
  for (iseg=1;iseg<=2*MT.number;iseg++)
    {
      iMT=segindex(iseg);
      if (seg.length[iseg]>0)
	fprintf(fmap,"# iMT.iseg %d.%d %f\n",MT.index[iMT],iseg-2*iMT+2,
		seg.cm[iseg].x);
    }
  /* List segment lengths at the bottom of the AGR file */
  agr_comment(fmap,"");
  agr_comment(fmap,"MT Segment Lengths");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"-----------------------------------------");
  for (iseg=1;iseg<=2*MT.number;iseg++)
    {
      iMT=segindex(iseg);
      if (seg.length[iseg]>0)
	fprintf(fmap,"# iMT.iseg %d.%d %f\n",MT.index[iMT],iseg-2*iMT+2,
		seg.length[iseg]);
    }

   /* List segment polarity at the bottom of the AGR file */
  agr_comment(fmap,"");
  agr_comment(fmap,"MT Segment Polarity");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"-----------------------------------------");
  for (iseg=1;iseg<=2*MT.number;iseg++)
    {
      iMT=segindex(iseg);
      if (seg.length[iseg]>0)
	fprintf(fmap,"# iMT.iseg %d.%d %f\n",MT.index[iMT],iseg-2*iMT+2,
		seg.direct[iseg].x);
    }
    
  fclose(fmap);
}

void wmap_yz(int filenum)
{
  int iMT,jMT,iseg,jseg,imin,imax,icolor,lc,gn,arrowtype;
  char filename[100]="MapYZ_";
  char title[100]="MT ARRAY FOR ITERATION No_"; 
  float ymin,zmin,ymax,zmax,MTrad,ycm,zcm;
  float lw,fillc,fillpt;
  FILE *fmap;
  struct point pmhead, pmleg;

  /* These files are writen as Graph number 0 */
  gn=0;

  /* Obtain the range of the graph */
  ymin=-0.5;
  ymax=0.5;
  zmin=-0.5;
  zmax=0.5;

  /* Open a seperate file for each filenum */
  strcat(filename,(const char *)itoa(filenum));
  strcat(filename,".dat");
  strcat(title,(const char *)itoa(filenum));

  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fmap=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);

 /* Write Standard configuration data For XMGR */
  agr_comment(fmap,"\n");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"MAP OF THE MT ARRAY"); 

  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis) */
  agr_std(fmap,gn,title,"y/\\f{12}m\\1m","z/\\f{12}m\\1m",ymin,zmin,ymax,zmax);
  
  /* Generate a plot of the MTs */
  /* Draw an outward circle that bounds the MT bundle */
  /* In that circle plot red and blue circles 
     for all MTs. */
 
  /* void agr_cyl(FILE *fp,
     float xcm,
     float zcm,
     float r,
     int lw,
     int lc,
     int fpt,
     int fc)*/
  
  lc=1; /* outer circle line color */
  lw=1; /* outer circle line width */
  fillpt=0; /* fill pattern */
  fillc=1;  /* fill color */
  agr_cyl(fmap,0,0,bundle.rad,lw,lc,fillpt,fillc);
  
  lw=1; /* line width */
  fillpt=4; /* fill pattern */
  /* determine radius of MT */ 
  if (exclude>0.0) 
    MTrad=exclude;
  else
    MTrad=sqrtf(bundle.rad/(float)MT.number);
  
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      ycm=MT.cm[iMT].y;
      zcm=MT.cm[iMT].z;
      if (MT.direct[iMT].x>0)
	{
	  lc=2;
	  fillc=2;
	  agr_cyl(fmap,ycm,zcm,MTrad,lw,lc,fillpt,fillc);
	}
      else
	{
	  lc=4;
	  fillc=4;
	  agr_cyl(fmap,ycm,zcm,MTrad,lw,lc,fillpt,fillc);
	}
    }
  /* Mark a number at the cm of each MT */

  for (iMT=1;iMT<=MT.number;iMT++)
    /*void agr_cyl(FILE *fp,float xcm,float zcm,
      float r,int lw,int lc,int fpt,int fc)*/
    /*agr_cyl(fmap,MT.cm[iMT].x,MT.cm[iMT].z,0.01,0,1,1,2);*/
    agr_string(fmap,gn,MT.cm[iMT].y,MT.cm[iMT].z,MT.index[iMT]);
  
  /* Generate a plot of motors at overlap regions */
  /* Do not plot motors in overlap region between segments of same MT */
  /*gn=0;*/
  lw=2.0;
  for (iseg=1;iseg<=2*MT.number-1;iseg++)
    {
      iMT=segindex(iseg);
    for (jseg=iseg+1;jseg<=2*MT.number;jseg++)
      {
	jMT=segindex(jseg);
	if (soverlap(iseg,jseg)>0 && segindex(iseg)!=segindex(jseg)
	    && ovlp.type[iseg][jseg]!=ZERO)
	  {
	    if (ovlp.type[iseg][jseg]==BIPOLAR)
	      {
		lc=11;
		lw=2.0;
		arrowtype=3;
		motorvec(iseg,jseg,&pmhead,&pmleg);
		agr_arrowline(fmap,gn,lw,lc,pmhead.y,pmhead.z,
			      pmleg.y,pmleg.z,arrowtype); 
	      }
	    else if (ovlp.type[iseg][jseg]==TWOMOTORS)
		{
		  lc=3; /* green */
		  lw=2.0;
		  arrowtype=2;
		  if (seg.direct[iseg].x*seg.direct[jseg].x>0)
		    {
		      /* parallel pair */
		      motorvec(iseg,jseg,&pmhead,&pmleg);
		      agr_arrowline(fmap,gn,lw,lc,pmhead.y,pmhead.z,
				    pmleg.y,pmleg.z,arrowtype);
		      agr_arrowline(fmap,gn,lw,lc,pmhead.y,pmleg.z,
				    pmleg.y,pmhead.z,arrowtype);
		    }
		  else
		    {
		      float dy=0.1;
		      /* anti-parallel pair */
		      motorvec(iseg,jseg,&pmhead,&pmleg);
		      agr_arrowline(fmap,gn,lw,lc,pmhead.y+dy,pmhead.z,
				    pmleg.y+dy,pmleg.z,arrowtype);
		      agr_arrowline(fmap,gn,lw,lc,pmhead.y-dy,pmleg.z,
				    pmleg.y-dy,pmhead.z,arrowtype);
		    }
		}
	    else if (ovlp.type[iseg][jseg]==BUNDLING /*&& 
		     seg.direct[iseg].x*seg.direct[jseg].x<0 &&
		     fabs(MT.pend[iMT].x-MT.pend[jMT].x)<3.0*/)
	      {
		lc=10; /* magenta */
		lw=4;
		arrowtype=0;
		motorvec(iseg,jseg,&pmhead,&pmleg);
		agr_arrowline(fmap,gn,lw,lc,pmhead.y,pmhead.z,
			      pmleg.y,pmleg.z,arrowtype);
	      }
	    else /* Unipolar motor, direction determined in motorvec */  
	      {
		lc=2;
		arrowtype=2;
		motorvec(iseg,jseg,&pmhead,&pmleg);
		agr_arrowline(fmap,gn,lw,lc,pmhead.y,pmhead.z,
			      pmleg.y,pmleg.z,arrowtype);
	      }
	  }
      }
    }
  fclose(fmap);
}

/* ----------------------------------------------------------------*/
/*The following function prints the av motor density at iter*/
void wav_motordens(int filenum){
    
    float xmax,xmin;
    
    FILE *fp;
    int ibox;
    float boxlength;
    
    xmin=av.xmin_real[iter];
    xmax=av.xmax_real[iter];
    boxlength=(xmax-xmin)/nbox;
    
    /*if kinesin initially there print output*/
    if (ProbKinesin>0) {
        
        /*open file with iter number in name*/
        char filename[100]="AvKindens_";
        strcat(filename,(const char *)itoa(filenum));
        strcat(filename,".dat");
        if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
        if ((fp=fopen(filept(filename),"w"))==NULL)
            printf("COULDNT OPEN FILE: %s",filename);
        
        fprintf(fp,"#box_cm [microns]\tKiProb%i\n",filenum);
        for (ibox=1; ibox<=nbox; ibox++) {
            fprintf(fp,"%f\t%f\n",xmin+(boxlength*ibox)/2, av.ovlp_act[iter][ibox]*av.ovlp_kin[iter][ibox]);
        }
        fclose(fp);
    }
    
    /*if Dynein initially there print output*/
    if (ProbKinesin<1) {
        
        /*open file with iter number in name*/
        char filename[100]="AvDyndens_";
        strcat(filename,(const char *)itoa(filenum));
        strcat(filename,".dat");
        if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
        if ((fp=fopen(filept(filename),"w"))==NULL)
            printf("COULDNT OPEN FILE: %s",filename);
        
        fprintf(fp,"#box_cm [microns]\tDyProb%i\n",filenum);
        for (ibox=1; ibox<=nbox; ibox++) {
            fprintf(fp,"%f\t%f\n",xmin+(boxlength*ibox)/2,av.ovlp_act[iter][ibox]*av.ovlp_dyn[iter][ibox]);
        }
        fclose(fp);
    }
    
    /*if Dynein initially there print output*/
    if (ProbBipolar>0) {
        
        /*open file with iter number in name*/
        char filename[100]="AvBidens_";
        strcat(filename,(const char *)itoa(filenum));
        strcat(filename,".dat");
        if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
        if ((fp=fopen(filept(filename),"w"))==NULL)
            printf("COULDNT OPEN FILE: %s",filename);
        
        fprintf(fp,"#box_cm [microns]\tBiProb%i\n",filenum);
        for (ibox=1; ibox<=nbox; ibox++) {
            fprintf(fp,"%f\t%f\n",xmin+(boxlength*ibox)/2, av.ovlp_act[iter][ibox]*av.ovlp_bi[iter][ibox]);
        }
        fclose(fp);
    }
}


/* ----------------------------------------------------------------*/

void wav_density (int filenum)
{
  /* This function calculates an averaged profile <nl(t,x)>, 
     <nr(t,x)> and <ntot(t,x)>; where nr and nl are the number 
     of MTs pointing with their minus-end to the right and left.
     Note that this happens completely oblivious about the MTs lengths
  */

  int iMT,jMT;
  char filename[100]="AvDensProf_";
  char title[100]="AV DENSITY PROFILE FOR ITERATION No_"; 
  float xmin,zmin,xmax,zmax,xrange,zrange,dx,x,px;
  FILE *fp;
  int ibox;

  /* These files are denoted graph number 1 == G1 */
  /*gn=1;*/


  xmin=av.xmin[iter];
  xmax=av.xmax[iter];

  zmin=-bundle.rad;
  zmax=bundle.rad;
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
  agr_std(fp,gn,title,"x/\\f{12}m\\1m","Averaged Density",xmin,zmin,xmax,zmax);
  
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
      x=ibox*dx+av.xmin[iter];
      px=av.nr[iter][ibox]/dx;
      agr_wcoor(fp,x,px);
    }
  
  agr_begset(fp,gn,1,4,1,2.0,0,0.0,2);
  for(ibox=0;ibox<nbox;ibox++)
    {
      x=ibox*dx+av.xmin[iter];
      px=av.nl[iter][ibox]/dx;
      agr_wcoor(fp,x,px);
    }
  
  agr_begset(fp,gn,2,1,1,2.0,0,0.0,3);  
  for(ibox=0;ibox<nbox;ibox++)
    {
      x=ibox*dx+av.xmin[iter];
      px=(av.nl[iter][ibox]+av.nr[iter][ibox])/dx;
      agr_wcoor(fp,x,px);
    }
  
  fclose(fp);
}

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
    float absvmin,absvmax,dv,norm[31];
    FILE *fp;
    
    if (remove (filept_max(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    fprintf(fp,"# %i averaging iterations\n ",nav);
    fprintf(fp,"# BLA \n");
    fprintf(fp,"#vel[micron/s]\t");
    for (i=2; i<=30;i++){
        fprintf(fp,"Avfractsize%i\t",i);
        norm[i]=0;
    }
    fprintf(fp, "\n");
  
    dv=maxvel/(float)nbox;

    /*normalize*/
    for (i=2; i<=30; i++) {
        for (ibox=1; ibox<=nbox; ibox++) {
            norm[i]+= av.nv_cluster[i][ibox];
        }
    }

    
    for (ibox=1; ibox<=nbox; ibox++) {
        
        fprintf(fp,"%f\t", (ibox-0.5)*dv);
        for (i=2; i<=30; i++) {
            fprintf(fp,"%f\t", av.nv_cluster[i][ibox]/norm[i]);
        }
        fprintf(fp,"\n");
    }
    
    fclose(fp);
}
/* ----------------------------------------------------------------*/

/*This function calculates the av abs vel*/
void wav_absvel(){
    int ibox,it,i;
    char filename[100]="Av_clustervel.dat";
    float absvmin,absvmax,dv, norm=0,avbox=0;
    FILE *fp;
    
    if (remove (filept_max(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    fprintf(fp,"# %i averaging iterations\n ",nav);
    fprintf(fp,"# BLA \n");
    fprintf(fp,"#vel[micron/s]\t");

    
    fprintf(fp, "\n");
    
    dv=maxvel/(float)nbox;
    
    /*normalize*/
    for (i=2; i<=30; i++) {
        for (ibox=1; ibox<=nbox; ibox++) {
            norm+= av.nv_cluster[i][ibox];
        }
    }
    
    
    for (ibox=1; ibox<=nbox; ibox++) {
        
        fprintf(fp,"%f\t", (ibox-0.5)*dv);
        avbox=0.;
        for (i=2; i<=30; i++) {
            avbox+= av.nv_cluster[i][ibox];
        }
        fprintf(fp,"%f\n", avbox/norm);
    }
    
    fclose(fp);
}
/* ----------------------------------------------------------------*/



void wav_percprob ()
{
  /* This function calculates the percolation probability of a give chi
  */

  int iMT,jMT,it, ibox;
  float dens=0., bundle_rad=0., av_max=0.,av_min=0.;
  char filename[100]="percprob_vs_chi.dat";
  FILE *fp;

  /*get average MT density per volume*/
    for (it=1; it<=niter; it++) {
       
        bundle_rad+= av.thick[it];

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
   
    
  if(iav==1){
      if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
      if ((fp=fopen(filept(filename),"a"))==NULL)
          printf("COULDNT OPEN FILE: %s",filename);
      fprintf(fp,"average radius= \t%f\n", bundle_rad);
  }
  else{
      if ((fp=fopen(filept(filename),"a"))==NULL)
          printf("COULDNT OPEN FILE: %s",filename);
  }
    
  fprintf(fp,"%f\t%f\n", ProbActive, perc_prob/niter);
    
    
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
    float f=0.;
    char filename[100]="Stallforce_vs_Chi.dat";
    FILE *fp;
    
    if(iav==1)
        if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
    
    if ((fp=fopen(filept(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    for (it=neqsteps+1; it<=niter; it++) {
        /*only count positive force*/
        if (track.rbound[it]>rbound0)
            f += track.rbound[it]-rbound0;
    }
    
    f /= (niter-neqsteps);
    
    fprintf(fp,"%f\t%f\n", ProbActive, f*spring_r);
    
    
    fclose(fp);
}


/* ----------------------------------------------------------------*/
void wav_MTorder(int filenum)
{
  /* This function calculates averaged MT order profile 
     s(x)=<nr(t,x)-nl(t,x)>/<nr(t,x)+nl(t,x)>; 
     s2(x)=s(x)*s(x);
     where nr and nl 
     are the number of MTs pointing with their minus-end 
     to the right and left. 
  */

  int iMT,jMT;
  char filename[100]="AvMTorder_";
  char title[100]="AV MT ORDER PROFILE FOR ITERATION No_"; 
  float xmin,zmin,xmax,zmax,xrange,zrange,dx,x,s;
  FILE *fp;
  int ibox;

  /* These files are denoted graph number 1 == G1 */
  /*gn=1;*/


  xmin=av.xmin[iter];
  xmax=av.xmax[iter];

  zmin=-bundle.rad;
  zmax=bundle.rad;
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
  agr_comment(fp,"ORDER PROFILE"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"x/\\f{12}m\\1m","Averaged Order",xmin,zmin,xmax,zmax);
  
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
      if ((av.nr[iter][ibox]+av.nl[iter][ibox])>0)
	{
	  x=ibox*dx+av.xmin[iter];
	  s=(av.nr[iter][ibox])/
	    (av.nr[iter][ibox]+av.nl[iter][ibox]);
	  agr_wcoor(fp,x,s);
	}
    }

  agr_begset(fp,gn,1,4,1,2.0,0,0.0,3);  
  for(ibox=0;ibox<nbox;ibox++)
    {
      if ((av.nr[iter][ibox]+av.nl[iter][ibox])>0)
	{
	  x=ibox*dx+av.xmin[iter];
	  s=(av.nl[iter][ibox])/
	    (av.nr[iter][ibox]+av.nl[iter][ibox]);
	  agr_wcoor(fp,x,s);
	}
    }

  agr_begset(fp,gn,1,1,1,2.0,0,0.0,3);  
  for(ibox=0;ibox<nbox;ibox++)
    {
      if ((av.nr[iter][ibox]+av.nl[iter][ibox])>0)
	{
	  x=ibox*dx+av.xmin[iter];
	  s=(av.nr[iter][ibox]-av.nl[iter][ibox])/
	    (av.nr[iter][ibox]+av.nl[iter][ibox]);
	  agr_wcoor(fp,x,sqrt(s*s));
	}
    }
  
  fclose(fp);
}

/* ----------------------------------------------------------------*/
void wav_OrderParam()
{
  /* This function calculates averaged order parameter 
     profile 
     S(t)=Int(s2(x,t)dx); where,
     s(x)=<nr(t,x)-nl(t,x)>/<nr(t,x)+nl(t,x)>; 
     s2(x)=s(x)*s(x);
     and 
     nr and nl are the number of MTs pointing with their minus-end 
     to the right and left. 
  */

  int iMT,jMT,it,icourse,i;
  char filename[100]="Av_OrderPram";
  char title[100]="AV MT ORDER PARAMETER PROFILE"; 
  float xmin,zmin,xmax,zmax,xrange,zrange,dx,x,s,op,t,nr,nl;
  FILE *fp;
  int ibox;

  /* These files are denoted graph number 1 == G1 */
  /*gn=1;*/

  zmin=0;
  zmax=1;

  /* Open output file */
  strcat(filename,".dat");
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"ORDER PROFILE"); 
  

  agr_std(fp,gn,title,"x/\\f{12}m\\1m","Averaged Order",0,zmin,niter*dt,zmax);
  

  
  agr_begset(fp,gn,0,2,1,2.0,0,0.0,1);

  for (it=1; it<=niter; it++)
    {
      t=it*dt;
      agr_wcoor(fp,t,av.op[it]);
    }
 
  fclose(fp);
}
/* ----------------------------------------------------------------*/
/*Get tubulin order parameter*/
void wav_tubulin_op(int filenum){
    
    FILE *fp;
    float xmin, xmax,dx, xbox, op;
    float av_max=0., av_min=0., bundle_rad=0., dens=0.;
    int ibox,it;
    char filename[100]="AvMTorder_";
    
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
    strcat(filename,".dat");
    
    if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
    if ((fp=fopen(filept(filename),"w"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    
    fprintf(fp,"radius= %f tubulin dens= %f\n",bundle_rad, dens);
    fprintf(fp,"box_it%i [microns]\tOPit%i\n", iter,iter);
    
    for (ibox=1; ibox<=nbox; ibox++) {
        xbox=(2*ibox-1)*dx/2+xmin;
        
        op= (av.tubulin_r[filenum][ibox]-av.tubulin_l[filenum][ibox])/(av.tubulin_r[filenum][ibox]+av.tubulin_l[filenum][ibox]);
        printf("op=%f\n",op);
        fprintf(fp,"%f\t%f\n", xbox, op);
    }
    
    fclose(fp);
    
}

/* -------------------------------------------------------------*/

void wav_velcor (int filenum)
{
  /* This function calculates a velocity/velocity correlation function
     <v(t,0)v(t,x> 
  */

  int iMT,jMT;
  char filename[100]="VelCorr_";
  char title[100]="VELOCITY CORRELATION FUNCTION FOR ITERATION No_"; 
  float xmin,zmin,xmax,zmax,xrange,zrange,dx,x,px;
  FILE *fp;
  int ibox;

  /* These files are denoted graph number 1 == G1 */
  /*gn=1;*/
  
  xmin=av.xmin[iter];
  xmax=av.xmax[iter];
  zmin=-bundle.rad;
  zmax=bundle.rad;
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
  agr_comment(fp,"VELOCITY CORRELATION FUNCTION"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"x/\\f{12}m\\1m","<V(0)V(x)>/<V2>",xmin,zmin,xmax,zmax);
  
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
  
  agr_begset(fp,gn,0,1,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      x=ibox*dx;
      px=av.velcor[iter][ibox]/av.velt2[iter];
      agr_wcoor(fp,x,px);
    }
    
  fclose(fp);
}

/* -------------------------------------------------------------*/

void wav_polarcor (int filenum)
{
  /* This function calculates a polarity correlation function
     <n(t,0)n(t,x> 
  */

  int iMT,jMT;
  char filename[100]="PolarCorr_";
  char title[100]="POLARITY CORRELATION FUNCTION FOR ITERATION No_"; 
  float xmin,zmin,xmax,zmax,xrange,zrange,dx,x,px;
  FILE *fp;
  int ibox;

  /* These files are denoted graph number 1 == G1 */
  /*gn=1;*/
  
  xmin=av.xmin[iter];
  xmax=av.xmax[iter];
  zmin=-1;
  zmax=1;
  dx=(xmax-xmin)/(float)(nbox-1);
  
  /* Open a seperate file for each filenum */
  strcat(filename,(const char *)itoa(filenum));
  /*strcat(filename,"iav_");
    strcat(filename,(const char *)itoa(iav));*/
  strcat(filename,".dat");
  strcat(title,(const char *)itoa(filenum));
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"POLARITY CORRELATION FUNCTION"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"x/\\f{12}m\\1m","<n(0)n(x)>",0,zmin,xmax,zmax);
  
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
  
  agr_begset(fp,gn,0,1,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      x=ibox*dx;
      px=av.polarcor[iter][ibox];
      agr_wcoor(fp,x,px);
    }
    
  fclose(fp);
}


/* -------------------------------------------------------------*/


/* -------------------------------------------------------------*/

void wav_polarcor_ryz (int filenum)
{
  /* This function calculates a radial polarity correlation function
     <n(t,0)n(t,ryz> 
  */

  int iMT,jMT;
  char filename[100]="RYZ_PolarCorr_";
  char title[100]="RADIAL(YZ) POLARITY CORRELATION FUNCTION FOR ITERATION No_"; 
  float ymin,zmin,ymax,zmax,rmin,rmax,xrange,zrange,dr,r,pr;
  FILE *fp;
  int ibox;

  /* These files are denoted graph number 1 == G1 */
  /*gn=1;*/
  
  ymin=av.ymin[iter];
  ymax=av.ymax[iter];
  
  zmin=av.zmin[iter];
  zmax=av.zmax[iter];
  
  rmin=0.0;
  rmax=sqrt((zmax-zmin)*(zmax-zmin)+(ymax-ymin)*(ymax-ymin));
  rmax=10*2*exclude;
  dr=(rmax-rmin)/(float)(nbox2-1);
  
  /* Open a seperate file for each filenum */
  strcat(filename,(const char *)itoa(filenum));
  /*strcat(filename,"iav_");
    strcat(filename,(const char *)itoa(iav));*/
  strcat(filename,".dat");
  strcat(title,(const char *)itoa(filenum));
  
  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fp=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);
  
  /* Write Standard configuration data For XMGR */
  agr_comment(fp,"\n");
  agr_comment(fp,"-----------------------------------------");
  agr_comment(fp,"POLARITY CORRELATION FUNCTION"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"r/\\f{12}m\\1m","<n(0)n(r)>",rmin,-1,rmax,1);
  
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
  for(ibox=0;ibox<nbox2;ibox++)
    {
      r=ibox*dr;
      pr=av.polarcor_ryz[iter][ibox];
      agr_wcoor(fp,r,pr);
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

void wav_totabsveldensity ()
{
  /* This function calculates an averaged profile <v_l(t,x)>, 
     <v_r(t,x)> and <v_t(t,x)>; where v_r and v_l are the number 
     of MTs traveling with velocity v, and whose minus-end points 
     to the right and left, respectively. 
  */

  int iMT,jMT,it;
  char filename[100]="Tot_absVelProf";
  char title[100]="TOT AV VELOCITY PROFILE FOR ITERATION No_"; 
  float absvmin,zmin,absvmax,zmax,dv,v,pv;
  FILE *fp;
  int ibox;
  
  absvmin=0.0;
  absvmax=0.0;
  for (it=1;it<=niter;it++)
    {
      absvmin+=av.absvmin[it];
      absvmax+=av.absvmax[it];
    }
  absvmin/=(float)niter;
  absvmax/=(float)niter;

  /* Uncomment to use fixed values for the min/max values of Vel */ 
  /* these are defined in rdata.c                                */
  /*
  absvmin=hist_absvel_min;
  absvmax=hist_absvel_max;
  */

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
      for (it=1;it<=niter;it++)
          pv+=av.absvel_nr[it][ibox];///dv;
      pv/=(float)niter;
      v=ibox*dv+absvmin;

      agr_wcoor(fp,v,pv);
    }

 agr_begset(fp,gn,1,4,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      pv=0.0;
      for (it=1;it<=niter;it++)
          pv+=av.absvel_nl[it][ibox];//dv;
      pv/=(float)niter;
      v=ibox*dv+absvmin;

      agr_wcoor(fp,v,pv);
    }

 agr_begset(fp,gn,2,1,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      pv=0.0;
      for (it=1;it<=niter;it++)
          pv+=(av.absvel_nr[it][ibox]+av.absvel_nl[it][ibox]);//dv;
      pv/=(float)niter;
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
  /*
  vmin=hist_vel_min;
  vmax=hist_vel_max;
  */

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

void wav_totveldensity ()
{
  /* This function calculates an averaged profile <v_l(t,x)>, 
     <v_r(t,x)> and <v_t(t,x)>; where v_r and v_l are the number 
     of MTs traveling with velocity v, and whose minus-end points 
     to the right and left, respectively. 
  */

  int iMT,jMT,it;
  char filename[100]="Tot_VelProf";
  char title[100]="TOT AV VELOCITY PROFILE FOR ITERATION No_"; 
  float vmin,zmin,vmax,zmax,dv,v,pv;
  FILE *fp;
  int ibox;
  
  vmin=0.0;
  vmax=0.0;
  for (it=1;it<=niter;it++)
    {
      vmin+=av.vmin[it];
      vmax+=av.vmax[it];
    }
  vmin/=(float)niter;
  vmax/=(float)niter;

  /* Uncomment to use fixed values for the min/max values of Vel */ 
  /* these are defined in rdata.c                                */
  /* vmin=hist_vel_min;
     vmax=hist_vel_max;*/

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
      for (it=1;it<=niter;it++)
	  pv+=av.vel_nr[it][ibox]/dv;
      pv/=(float)niter;
      v=ibox*dv+vmin;
      agr_wcoor(fp,v/vel0_d,pv);
    }

      
   agr_begset(fp,gn,1,4,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      pv=0.0;
      for (it=1;it<=niter;it++)
	  pv+=av.vel_nl[it][ibox]/dv;
      pv/=(float)niter;
      v=ibox*dv+vmin;
      agr_wcoor(fp,v/vel0_d,pv);
    }
      
  agr_begset(fp,gn,2,1,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      pv=0.0;
      for (it=1;it<=niter;it++)
	  pv+=(av.vel_nr[it][ibox]+av.vel_nl[it][ibox])/dv;
      pv/=(float)niter;
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
  
  float dens=0., bundle_rad=0., av_max=0.,av_min=0.;
  int it, ibox;
  char filename[100]="AvEnd2End";
  char title[100]="AV END TO END DISTANCE AS FUNCTION OF TIME"; 
  FILE *fp;
  
    
    /*get average MT density per volume*/
    for (it=1; it<=niter; it++) {
        //printf("thick=%f\n", av.thick[it]);
        bundle_rad+= av.thick[it];
        
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
  agr_comment(fp,"AVERAGE END2END");
    fprintf(fp,"dens=%f\tradius=%f \n",dens, bundle_rad);

  /*agr_std(fp,gn,title,"time (sec)","END to End Distance \\f{12}m\\1m"
	  ,0,0,niter*dt,av.xmax_real[niter]-av.xmin_real[niter]);*/
 

  //agr_begset(fp,gn,0,1,1,2.0,0,0.0,1);
    for(it=1;it<=niter;it++){
        
        agr_wcoor(fp,it*dt,av.xmax_real[it]-av.xmin_real[it]);
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

void wav_Rg()
{
  /* This function calculates the average distance 
     to the most distant MT in the shaft. 
  */

  int it;
  char filename[100]="Av_Rg";
  char title[100]="Radius of Gyration as a function of time"; 
  FILE *fp;
  
  /* These files are denoted graph number 5 -- G5 */ 
  /*gn=5;*/

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
  agr_comment(fp,"AVERAGE Rg-POS"); 
  
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
    if (PERIOD)
      for(it=neqsteps;it<=niter;it++)
	agr_wcoor(fp,(it-neqsteps)*dt,sqrt(av.rxr2[it]-av.rxr[it]*av.rxr[it]));
    else
      for(it=1;it<=niter;it++)
	agr_wcoor(fp,it*dt,sqrt(av.xr2[it]-av.xr[it]*av.xr[it]));
  
  agr_begset(fp,gn,1,4,1,2.0,0,0.0,1);
  if (PERIOD)
    for(it=neqsteps;it<=niter;it++)
      agr_wcoor(fp,(it-neqsteps)*dt,sqrt(av.rxl2[it]-av.rxl[it]*av.rxl[it]));
  else
    for(it=1;it<=niter;it++)
      agr_wcoor(fp,it*dt,sqrt(av.xl2[it]-av.xl[it]*av.xl[it]));
  
  agr_begset(fp,gn,2,1,1,2.0,0,0.0,1);
  if (PERIOD)
    for(it=neqsteps;it<=niter;it++)
      agr_wcoor(fp,(it-neqsteps)*dt,sqrt(av.rxt2[it]-av.rxt[it]*av.rxt[it]));
  else
    for(it=1;it<=niter;it++)
      agr_wcoor(fp,it*dt,sqrt(av.xt2[it]-av.xt[it]*av.xt[it]));
  fclose(fp);
  

}

void wav_RgYZ()
{
  /* This function calculates the average distance 
     to the most distant MT in the shaft. 
  */

  int it;
  char filename[100]="Av_RgYZ";
  char title[100]="Bundle Thickness as a function of time"; 
  FILE *fp;
  
  /* These files are denoted graph number 5 -- G5 */ 
  /*gn=5;*/

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
  agr_comment(fp,"AVERAGE RgYZ-POS"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"time (sec)","<ryz2(t)>\\f{12}m\\1m"
	  ,0,0,niter*dt,bundle.rad);
  
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
      agr_wcoor(fp,it*dt,av.ryz[it]);

 
  fclose(fp);


}

/* WRITE XMGR FILE: ARRAY OF MTs */
 
void wvelmap(int filenum)
{
  int iMT,jMT,icolor,lc,gn,arrowtype;
  int ivelmin,ivelmax,nvel;
  float fend,bend,size,dv;
  float velmin,velmax;
  char filename[100]="VelFld_";
  char title[100]="VELOCITY FIELD FOR ITERATION No_"; 
  float xmin,zmin,xmax,zmax,scale;
  float lw;
  FILE *fmap;
  struct point pmhead, pmleg;

  /* These files are writen as Graph number 0 */
  gn=0;

  /* find first and last MT to obtain the range of the graph */
  xmin=MT.xmin-0.5*MT.length[MT.imin];
  xmax=MT.xmax+0.5*MT.length[MT.imax];
  zmin=-0.5;
  zmax=0.5;

  /* Open a seperate file for each filenum */
  strcat(filename,(const char *)itoa(filenum));
  strcat(filename,".dat");
  strcat(title,(const char *)itoa(filenum));

  if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
  if ((fmap=fopen(filept(filename),"w"))==NULL) 
    printf("COULDNT OPEN FILE: %s",filename);

 /* Write Standard configuration data For XMGR */
  agr_comment(fmap,"\n");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"MT VELOCITY FIELD"); 

  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis) */
  agr_std(fmap,gn,title,"x/\\f{12}m\\1m","y/\\f{12}m\\1m",xmin,zmin,xmax,zmax);


  /* Sort all MT.vel in nvel levels; the sloweset is 1, 
     the highest is nvel                               */
  /* Find MTs with min and max velocities              */
  
  nvel=100;
  
/* ivelmin=1; */
/*   for (iMT=1;iMT<=MT.number;iMT++) */
/*     if (fabs(MT.vel[iMT].x)<fabs(MT.vel[ivelmin].x)) ivelmin=iMT; */
/*   ivelmax=1; */
/*  for (iMT=1;iMT<=MT.number;iMT++)*/
/*    if (fabs(MT.vel[iMT].x)>fabs(MT.vel[ivelmax].x)) ivelmax=iMT;*/
/* velmin= fabs(MT.vel[ivelmin].x); */
/* velmax= fabs(MT.vel[ivelmax].x); */  
  

  /*velmin=0.0;
  velmax=10*vel0_d;
  dv=(velmax-velmin)/(nvel-1);*/

  scale=0.03*(xmax-xmin);
  dv=0.5*vel0_d; /* 1microM*scale in the plot is dv */

  
  /* Generate a plot of the MTs */
  /*gn=0;*/ /* graph number */
  lw=2; /* line width */
  arrowtype=2; /* arrow at end */

  for (iMT=1;iMT<=MT.number;iMT++)
    {
      /*size=fabs(fabs(MT.vel[iMT].x)-velmin)/dv;*/
      size=fabs(MT.vel[iMT].x)/dv;
      if (MT.vel[iMT].x<0) 
	{
	  bend=MT.cm[iMT].x+0.5*size*scale;
	  fend=MT.cm[iMT].x-0.5*size*scale;
	}
      else
	{
	  bend=MT.cm[iMT].x-0.5*size*scale;
	  fend=MT.cm[iMT].x+0.5*size*scale;
	}
      
      if (MT.direct[iMT].x>0) 
	lc=2;
      else
	lc=4;

      agr_arrowline(fmap,gn,lw,lc,bend, MT.cm[iMT].z, 
		    fend, MT.cm[iMT].z,arrowtype);      
    }

  {
    float xscale,zscale;
    xscale=xmax+(xmax-xmin)/10;
    zscale=zmax;
    agr_arrowline(fmap,gn,2,6,xscale-0.5*scale,zscale, 
		    xscale+0.5*scale,zscale,arrowtype);

    fprintf(fmap,"@with string\n");
    fprintf(fmap,"@    string on\n");
    fprintf(fmap,"@    string loctype world\n");
    fprintf(fmap,"@    string g%d\n",gn);
    fprintf(fmap,"@    string %f, %f\n",xscale-0.5*scale,0.95*zscale);
    fprintf(fmap,"@    string color 6\n");
    fprintf(fmap,"@    string rot 0\n");
    fprintf(fmap,"@    string font 0\n");
    fprintf(fmap,"@    string just 0\n");
    fprintf(fmap,"@    string char size 1.000000\n");
    fprintf(fmap,"@    string def \"%6.4f \\f{12}m\\1m/sec\"\n",dv);
  }

  /* mark circle at the CM */
  /*gn=0;*/ /* graph number */
  for (iMT=1;iMT<=MT.number;iMT++)
    /*void agr_cyl(FILE *fp,float xcm,float zcm,
      float r,int lw,int lc,int fpt,int fc)*/
    /*agr_cyl(fmap,MT.cm[iMT].x,MT.cm[iMT].z,0.01,0,1,1,2);*/
    agr_string(fmap,gn,MT.cm[iMT].x,MT.cm[iMT].z,MT.index[iMT]);
       

  /* List MT velocities at the bottom of the AGR file */
  agr_comment(fmap,"");
  agr_comment(fmap,"MT Velocities");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"-----------------------------------------");
  for (iMT=1;iMT<=MT.number;iMT++)
    fprintf(fmap,"# v[%d]/v0 = %f\n",MT.index[iMT],MT.vel[iMT].x/vel0_d);
  
  fclose(fmap);
}




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
  for(it=1;it<=niter;it++)
      agr_wcoor(fp,it*dt,av.velr[it]/vel0_d);

  agr_begset(fp,gn,1,4,1,2.0,0,0.0,1);
  for(it=1;it<=niter;it++)
      agr_wcoor(fp,it*dt,av.vell[it]/vel0_d);

  agr_begset(fp,gn,2,1,1,2.0,0,0.0,1);
  for(it=1;it<=niter;it++)
    agr_wcoor(fp,it*dt,av.velt[it]/vel0_d);
  
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

void wav_nzeros()
{
  /* This function writes the average fraction of MTs with zero neighbors
     as a function of time. 
  */

  
  int it;
  char filename[100]="Av_nzeros";
  char title[100]="Average fraction of MTs with zero neighbors"; 
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
  agr_comment(fp,"Fraction of MTs"); 
  
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
	  ,0,0,niter*dt,av.nzeros[niter]);
  
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
      agr_wcoor(fp,it*dt,av.nzeros[it]/av.nMT[it]);
  
  fclose(fp);


}

/* ----------------------------------------------------*/

void wbundle_shape (int filenum)
{
  /* This function prints the mean bundle radius as a function of X
  */

  int iMT,jMT;
  char filename[100]="BundleShape_";
  char title[100]="BUNDLE SHAPE --  ITERATION No_"; 
  float xmin,zmin,xmax,zmax,xrange,zrange,dx,x,px;
  FILE *fp;
  int ibox;

  /* These files are denoted graph number 1 == G1 */
  /*gn=1;*/
  
  xmin=av.xmin[iter];
  xmax=av.xmax[iter];
  zmin=-bundle.rad;
  zmax=bundle.rad;
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
  agr_comment(fp,"BUNDLE SHAPE"); 
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title, 
     const char *Xaxis, 
     const char *Yaxis,
     float xmin, 
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"x/\\f{12}m\\1m","<r>/\\f{12}m\\1m"
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
  
  agr_begset(fp,gn,0,1,1,2.0,0,0.0,1);
  for(ibox=0;ibox<nbox;ibox++)
    {
      x=ibox*dx+av.xmin[iter];
      if (av.rmin[iter][ibox] != bundle.rad)
	agr_wcoor(fp,x,av.rmin[iter][ibox]);
    }
  
  agr_begset(fp,gn,1,2,1,2.0,0,0.0,2);
  for(ibox=0;ibox<nbox;ibox++)
    {
      x=ibox*dx+av.xmin[iter];
      if (av.rmax[iter][ibox] != -bundle.rad)
	agr_wcoor(fp,x,av.rmax[iter][ibox]);
    }
  
  fclose(fp);
}

/* -------------------------------------------------------------*/
/* -------------------------------------------------------------*/


void wav_timevelcor (int filenum)
{
  /* This function calculates a velocity/velocity correlation function
     C(t,tau)=<v(t)v(t+tau)>
  */

  char filename[100]="TimeVelCorr_";
  char title[200]="TIME VELOCITY CORRELATION FUNCTIONS\\n(red ++,purple +-, blue -- black total)\\nFOR ITERATION No_\\n(Five last sets are corresponding <vv> products and average velocities\\n";
  float tmin,cmin,tmax,cmax,c,c0,tau;
  int it,idt;
  FILE *fp;

  /* These files are denoted graph number 1 == G1 */
  /*gn=1;*/
  
  tmin=0.0;
  tmax=niter*dt;
  cmin=0;
  cmax=16*vel0_d*vel0_d;
  it=filenum;

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
  agr_comment(fp,"TIME VELOCITY CORRELATION FUNCTION");
  
  /* void agr_std(FILE *fp,
     int gn,
     const char *Title,
     const char *Xaxis,
     const char *Yaxis,
     float xmin,
     float ymin,
     float xmax,
     float ymax) */
  agr_std(fp,gn,title,"\\f{12}t\\1 (sec)","<(V(t)V(t+\\f{12}t\\1)>/<V\\S2\\N>",
	  tmin,cmin,tmax,cmax);
  
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

  agr_begset(fp,gn,1,2,1,2.0,0,0.0,1);
  c0=av.velcort_pp[it][0]-av.velr[it]*av.velr[it];
  for(idt=0;idt<=(corr_tf-it);idt++)
    {
      tau=idt*dt;
      c=av.velcort_pp[it][idt]-av.velr[it+idt]*av.velr[it+idt];
      if (c0!=0) agr_wcoor(fp,tau,c/c0);
    }
  
  agr_begset(fp,gn,2,4,1,2.0,0,0.0,1);
  c0=av.velcort_mm[it][0]-av.vell[it]*av.vell[it];
  for(idt=0;idt<=(corr_tf-it);idt++)
    {
      tau=idt*dt;
      c=av.velcort_mm[it][idt]-av.vell[it+idt]*av.vell[it+idt];
      if (c0!=0) agr_wcoor(fp,tau,c/c0);
    }

  // Cross correlation.
  // Since its magnitude is found to be zero at all times 
  // we do not normalize it by the value, c0 at idt=0. 
  agr_begset(fp,gn,3,8,1,2.0,0,0.0,1);
  c0=av.velcort_pm[it][0]-av.vell[it]*av.velr[it];
  for(idt=0;idt<=(corr_tf-it);idt++)
    {
      tau=idt*dt;
      c=av.velcort_pm[it][idt]-av.vell[it+idt]*av.velr[it+idt];
      agr_wcoor(fp,tau,c);
    }

  // Plot the velocity products only
  agr_begset(fp,gn,4,2,1,1.0,0,0.0,1);
  for(idt=0;idt<=(corr_tf-it);idt++)
    {
      tau=idt*dt;
      c=av.velcort_pp[it][idt];
      agr_wcoor(fp,tau,c);
    }
  
  agr_begset(fp,gn,5,4,1,1.0,0,0.0,1);
  for(idt=0;idt<=(corr_tf-it);idt++)
    {
      tau=idt*dt;
      c=av.velcort_mm[it][idt];
      agr_wcoor(fp,tau,c);
    }

  agr_begset(fp,gn,6,8,1,1.0,0,0.0,1);
  for(idt=0;idt<=(corr_tf-it);idt++)
    {
      tau=idt*dt;
      c=av.velcort_pm[it][idt];
      agr_wcoor(fp,tau,c);
    }

// Plot the average velocity only
  agr_begset(fp,gn,7,11,1,1.0,0,0.0,1);
  for(idt=0;idt<=(corr_tf-it);idt++)
    {
     agr_wcoor(fp,idt*dt,av.velr[it+idt]); 
    }
  
  agr_begset(fp,gn,8,9,1,1.0,0,0.0,1);
  for(idt=0;idt<=(corr_tf-it);idt++)
    {
     agr_wcoor(fp,idt*dt,av.vell[it+idt]); 
    }
  
  fclose(fp);
}


/* -------------------------------------------------------------*/
/* void wav_xtraject() */
/* { */
/*   /\* This function calculates single MT trajectories */
/*   *\/ */

/*   char filename[100]="MT_Traject"; */
/*   char title[100]="Trajectories"; */
/*   float tmin,xmin,tmax,xmax,xcm,time; */
/*   int it,iMT,lc,n; */
/*   FILE *fp; */

/*   /\* These files are denoted graph number 1 == G1 *\/ */
/*   /\*gn=1;*\/ */
  
/*   tmin=0.0; */
/*   tmax=niter*dt; */
/*   xmin=-1000; */
/*   xmax=1000; */

/*   strcat(filename,".dat"); */
  
/*   if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename); */
/*   if ((fp=fopen(filept(filename),"w"))==NULL) */
/*     printf("COULDNT OPEN FILE: %s",filename); */
  
/*   /\* Write Standard configuration data For XMGR *\/ */
/*   agr_comment(fp,"\n"); */
/*   agr_comment(fp,"-----------------------------------------"); */
/*   agr_comment(fp,"SINGLE MT TRAJECTORIES"); */
  
/*   /\* void agr_std(FILE *fp, */
/*      int gn, */
/*      const char *Title, */
/*      const char *Xaxis, */
/*      const char *Yaxis, */
/*      float xmin, */
/*      float ymin, */
/*      float xmax, */
/*      float ymax) *\/ */
/*   agr_std(fp,gn,title,"t (sec)","x (\\f{12}m\\1m)", */
/* 	  xmin,tmin,xmax,tmax); */
  
/*   /\* void agr_begset(FILE *fp, */
/*      int gn, */
/*      iset k, */
/*      int lncolor, */
/*      int lnstyle, */
/*      float lnwidth, */
/*      int symbol, */
/*      float symbsz, */
/*      int symbcl) */
/*   *\/ */

/*   /\* number of trajectories to draw *\/ */
/*   if (MT.number<50)  */
/*     n=MT.number; */
/*   else */
/*     n=50; */

/*   for (iMT=1;iMT<=n;iMT++) */
/*     { */
/*       if (MT.direct[iMT].x>0)  */
/* 	lc=2; */
/*       else */
/* 	lc=4; */
				
/*       agr_begset(fp,gn,iMT,lc,1,2.0,0,0.0,1); */
/*       for(it=0;it<=niter;it++) */
/* 	{ */
/* 	  xcm=track.xcm[it][iMT]; */
/* 	  time=dt*it; */
/* 	  agr_wcoor(fp,time,xcm); */
/* 	} */
/*     } */
/*   fclose(fp); */
/* } */


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
    int i,it,jMT,t0,n,q;
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
        for (jMT=1; jMT<=track.nMT[it]*2; jMT++) {
            if (track.length[it][jMT]>0.001) {
                if (track.direct[it][jMT] < -0.1) {fprintf(fp,"Blue\t");
                    //printf("num:%i %f\n",jMT, track.direct[it][jMT]);
                }
                if (track.direct[it][jMT] > 0.1) {fprintf(fp,"Red\t");
                    //printf("num:%i %f\n", jMT, track.direct[it][jMT]);
                }
                
                if (track.direct[it][jMT] < -1.1  ) printf("ERROR");
                if (track.direct[it][jMT] > 1.1  ) printf("ERROR");
                if (track.direct[it][jMT] > -0.1 && track.direct[it][jMT] < 0.1 ) printf("ERROR");
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
            if (track.cross_type[it][i]==1)
                fprintf(fp, "Red\t" );
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
        for (jMT=1; jMT<=2*track.nMT[it]; jMT++) {
            if (track.length[it][jMT]>0.001) {
                xcm=track.xcm_vis[it][jMT];
                ycm=(track.ycm_vis[it][jMT]+track.zcm_vis[it][jMT])/sqrt(2.);//just for dramatic purposes (rotates system around 45 degrees so that we can see more MT in 2D)
                length=track.length[it][jMT];
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
    int i,it,jMT,t0,n,q;
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
        for (jMT=1; jMT<=track.nMT[it]*2; jMT++) {
            if (track.length[it][jMT]>0.001) {
                if (track.direct[it][jMT] < -0.1) {fprintf(fp,"Blue\t");
                    //printf("num:%i %f\n",jMT, track.direct[it][jMT]);
                }
                if (track.direct[it][jMT] > 0.1) {fprintf(fp,"Red\t");
                    //printf("num:%i %f\n", jMT, track.direct[it][jMT]);
                }
                
                if (track.direct[it][jMT] < -1.1  ) printf("ERROR");
                if (track.direct[it][jMT] > 1.1  ) printf("ERROR");
                if (track.direct[it][jMT] > -0.1 && track.direct[it][jMT] < 0.1 ) printf("ERROR");
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
        for (jMT=1; jMT<=2*track.nMT[it]; jMT++) {
            if (track.length[it][jMT]>0.001) {
                xcm=track.xcm_vis[it][jMT];
                ycm=track.ycm_vis[it][jMT];//just for dramatic purposes (rotates system around 45 degrees)
                zcm=track.zcm_vis[it][jMT];//that we can see more MT in 2D)
                length=track.length[it][jMT];
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
    float bound_av=0, f;
    FILE *fp;
    
    
    sprintf(filename,"polarity_exertedforce.dat");
    if (iav==1){
        if (remove (filept(filename))==0) printf("Removed old file:%s\n",filename);
        if ((fp=fopen(filept(filename),"w"))==NULL)
            printf("COULDNT OPEN FILE: %s",filename);
        fprintf(fp,"#polarity ratio \t exerted force [pN]\n");
        fclose(fp);
    }
    else{
        if ((fp=fopen(filept(filename),"a"))==NULL)
        printf("COULDNT OPEN FILE: %s",filename);
    }
    
    for (it=neqsteps+1; it<=niter; it++) {
        bound_av+=track.rbound[it];
    }
    bound_av /= niter-neqsteps;
	if (bound_av> rbound0 )
		f= spring_r*(bound_av-rbound0);
	else
		f=0.;
    printf("%f %f\n",PolarityRatio0,f);
    fprintf(fp, "%f\t%f\n", PolarityRatio0, f );
    
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
