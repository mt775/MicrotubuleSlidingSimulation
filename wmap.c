#include <string.h>
#include <stdio.h>
#include <math.h>
#include "global_var.h"
/* #include "func_declaration.h" */


/* WRITE XMGR FILE: ARRAY OF MTs */
 
void wmap_xz(int filenum)
{
  int iMT,jMT,imin,imax,icolor,lc,gn,fc,fpt;
  char filename[100]="MapXZ_";
  char title[100]="MT ARRAY FOR ITERATION No_"; 
  float xmin,zmin,xmax,zmax;
  float lw;
  FILE *fmap;
  struct point pmhead, pmleg;

  /* These files are writen as Graph number 0 */
  gn=0;

  /* Obtain the range of the graph */
  xmin=MT.xmin-0.5*MT.length[MT.imin];
  xmax=MT.xmax+0.5*MT.length[MT.imax];
  zmin=-bundle.rad;
  zmax=bundle.rad;

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
  
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      agr_arrowline(fmap,gn,lw,lc,MT.pend[iMT].x, MT.pend[iMT].z, 
		    MT.mend[iMT].x, MT.mend[iMT].z);
    }

  /* Generate a filled box around all MTs */     
  lw=1;  
  fpt=5;
  for (iMT=1;iMT<=MT.number;iMT++)
    {
      if (MT.direct[iMT].x>0)
	{
	  fc=2;
	  lc=2;
	}
      else
	{
	  fc=4;
	  lc=4;
	}	
      agr_box(fmap,MT.cm[iMT].x,MT.cm[iMT].z,MT.length[iMT],2*exclude,lw,
	      lc,fpt,fc);      
    }  
  
  /* mark circle at the CM */
  /*gn=0;*/ /* graph number */
  for (iMT=1;iMT<=MT.number;iMT++)
    /*void agr_cyl(FILE *fp,float xcm,float zcm,
      float r,int lw,int lc,int fpt,int fc)*/
    /*agr_cyl(fmap,MT.cm[iMT].x,MT.cm[iMT].z,0.01,0,1,1,2);*/
    agr_string(fmap,gn,MT.cm[iMT].x,MT.cm[iMT].z,MT.index[iMT]);
  
  /* Generate a plot of motors at overlap regions */
  /*gn=0;*/
  lc=2; 
  lw=2.0;
  for (iMT=1;iMT<=MT.number-1;iMT++)
    for (jMT=iMT+1;jMT<=MT.number;jMT++)
      {
	if (overlap(iMT,jMT)>0)
	  {
	    motorvec(iMT,jMT,&pmhead,&pmleg);
	    agr_arrowline(fmap,gn,lw,lc,pmhead.x,pmhead.z,
			  pmleg.x,pmleg.z);
	  }
      }
  /* Print the overlap matrix at the bottom of the AGR file */
  agr_comment(fmap,"");
  agr_comment(fmap,"OVERLAP TYPE MATRIX");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"0. = no ovlp, 1 = leg up, 2 = leg down");
  agr_comment(fmap,"-----------------------------------------");
  for (iMT=1;iMT<=MT.number-1;iMT++)
    {
      fprintf(fmap,"# %d  -  ",iMT);
      for (jMT=iMT+1;jMT<=MT.number;jMT++)
	fprintf(fmap,"%d ",ovlp.type[iMT][jMT]);
      fprintf(fmap,"\n");
    }
  
  /* Print the overlap matrix at the bottom of the AGR file */
  agr_comment(fmap,"");
  agr_comment(fmap,"OVERLAP MATRIX");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"-----------------------------------------");
  for (iMT=1;iMT<=MT.number-1;iMT++)
    {
      fprintf(fmap,"# %d  -  ",iMT);
      for (jMT=iMT+1;jMT<=MT.number;jMT++)
	fprintf(fmap,"%6.4f ",overlap(iMT,jMT));
      fprintf(fmap,"\n");
    }

  /* List MT velocities at the bottom of the AGR file */
  agr_comment(fmap,"");
  agr_comment(fmap,"MT Velocities");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"-----------------------------------------");
  for (iMT=1;iMT<=MT.number;iMT++)
      fprintf(fmap,"# v[%d] = %f\n",iMT,MT.vel[iMT].x);

/* List MT x-coordinate at the bottom of the AGR file */
  agr_comment(fmap,"");
  agr_comment(fmap,"MT x-Coordinates");
  agr_comment(fmap,"-----------------------------------------");
  agr_comment(fmap,"-----------------------------------------");
  for (iMT=1;iMT<=MT.number;iMT++)
      fprintf(fmap,"# %d %f\n",iMT,MT.cm[iMT].x);

     
  fclose(fmap);
}

void wmap_yz(int filenum)
{
  int iMT,jMT,imin,imax,icolor,lc,gn;
  char filename[100]="MapYZ_";
  char title[100]="MT ARRAY FOR ITERATION No_"; 
  float ymin,zmin,ymax,zmax,MTrad,ycm,zcm;
  float lw,fillc,fillpt;
  FILE *fmap;
  struct point pmhead, pmleg;

  /* These files are writen as Graph number 0 */
  gn=0;

  /* Obtain the range of the graph */
  ymin=-bundle.rad;
  ymax=bundle.rad;
  zmin=-bundle.rad;
  zmax=bundle.rad;

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
  /*gn=0;*/
  lc=2; 
  lw=2.0;
  for (iMT=1;iMT<=MT.number-1;iMT++)
    for (jMT=iMT+1;jMT<=MT.number;jMT++)
      {
	if (overlap(iMT,jMT)>0)
	  {
	    motorvec(iMT,jMT,&pmhead,&pmleg);
	    agr_arrowline(fmap,gn,lw,lc,pmhead.y,pmhead.z,
			  pmleg.y,pmleg.z);
	  }
      }     
  fclose(fmap);
}

