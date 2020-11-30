#include <sys/types.h>
#include <stdio.h>
#include <math.h>
#include "global_var.h" 


/* #################################################################### */
/* SET OF FUNCTIONS FOR WRITING XMGR FILES                              */
void agr_mathematica(FILE *fp){
    fprintf(fp,"t\t\tx\t\ty\t\tlength\n");
}

void agr_std(FILE *fp, int gn, const char *Title, const char *Xaxis, const char *Yaxis, float minX, float minY, float maxX, float maxY)
{
  
  float dx=*Xaxis/10.0;
  fprintf(fp,"# Grace project file \n");
  fprintf(fp,"@version 50119\n");
  fprintf(fp,"@ title  \"%s\"\n",Title);
  fprintf(fp,"@ title font 0\n");
  fprintf(fp,"@ title size 0.8\n");
  fprintf(fp,"@ title color 1\n");
  /*  fprintf(fp,"@ subtitle \"CN=%d NX=%d NZ=%d Rpep=%3.1f Pent=%3.1f Wall=%d Relax=%d\" \n",
      CN,NX,NZ,pep1.radi,pep1.pent,wall,relax);*/
  fprintf(fp,"@ subtitle font 0\n");
  fprintf(fp,"@ subtitle size 0.7\n");
  fprintf(fp,"@ subtitle color 1 \n");
  fprintf(fp,"@ xaxis  label \"%s\" \n",Xaxis);
  fprintf(fp,"@ xaxis label char size 1.0\n");
  fprintf(fp,"@ xaxis ticklabel char size 1.0\n");
  fprintf(fp,"@ xaxis tick size 1.0\n");
  fprintf(fp,"@ xaxis tick minor size 2.0\n");
  fprintf(fp,"@ xaxis tick major %f\n",dx);
  fprintf(fp,"@ yaxis  label \"%s\" \n",Yaxis);
  fprintf(fp,"@ yaxis label char size 1.0\n");
  fprintf(fp,"@ yaxis ticklabel char size 1.0\n");
  fprintf(fp,"@ yaxis tick size 2.0\n");
  fprintf(fp,"@ yaxis tick minor size 1.0\n");
  fprintf(fp,"@ yaxis tick major 1\n");
  fprintf(fp,"@ world xmin %f\n",minX);
  fprintf(fp,"@ world xmax %f\n",maxX);
  fprintf(fp,"@ world ymin %f\n",minY);
  fprintf(fp,"@ world ymax %f\n",maxY);
  fprintf(fp,"@ view xmin 0.115909 \n");
  fprintf(fp,"@ view xmax 0.888636 \n");
  fprintf(fp,"@ view ymin 0.153934 \n");
  fprintf(fp,"@ view ymax 0.872295 \n");
  fprintf(fp,"@ g%d on\n",gn);
  fprintf(fp,"@ g%d hidden false\n",gn);
}

void agr_legend(FILE *fp, int k, const char *name)
{
  if (k==1) 
    {
      fprintf(fp,"@ legend on\n");
      /*fprintf(fp,"@ legend loctype view\n");*/
      fprintf(fp,"@ legend 0.6, 0.8\n");
    } 
  
  fprintf(fp,"@ s%d legend  \"%s\" \n",k,name);
}
  
void agr_begset(FILE *fp, int gn, int k, int lncolor, int lnstyle, float lnwidth, int symbol, float symbsz, int symbcl)
{
  /*fprintf(fp,"@ s%d line color %d \n",k,lncolor);
  fprintf(fp,"@ s%d line linestyle %d \n",k,lnstyle);
  fprintf(fp,"@ s%d line linewidth %f \n",k,lnwidth);*/

  fprintf(fp,"@ s%d color %d \n",k,lncolor);
  fprintf(fp,"@ s%d linestyle %d \n",k,lnstyle);
  fprintf(fp,"@ s%d linewidth %f \n",k,lnwidth);
  fprintf(fp,"@ s%d symbol %d \n",k,symbol);
  fprintf(fp,"@ s%d symbol size %f \n",k,symbsz);
  fprintf(fp,"@ s%d symbol color %d \n",k,symbcl);
  fprintf(fp,"@target G%d.S%d \n",gn,k);
  fprintf(fp,"@type xy \n");
}

void agr_endset(FILE *fp)
{
  fprintf(fp,"&\n"); 
}

void agr_comment(FILE *fp, const char *comment)
{
  fprintf(fp,"# %s\n",comment);
}

void agr_wcoor(FILE *fp, float x, float y)
{
  /*fprintf(fp,"%16.4f %16.4f\n",x,y);*/
  fprintf(fp,"%16.9f %16.9f\n",x,y);
}

void agr_math_coord(FILE *fp,float time, float x, float y, float length)
{
    /*fprintf(fp,"%16.4f %16.4f\n",x,y);*/
    fprintf(fp,"%16.9f %16.9f %16.9f %16.9f\n",time,x,y,length);
}


void agr_wcoor3D(FILE *fp, float x, float y, float z)
{
  fprintf(fp,"%16.4f %16.4f %16.4f\n",x,y,z);
}

void agr_string(FILE *fp, int gn, float xcm, float zcm, int num)
{
  fprintf(fp,"@with string\n");
  fprintf(fp,"@    string on\n");
  fprintf(fp,"@    string loctype world\n");
  fprintf(fp,"@    string g%d\n",gn);
  fprintf(fp,"@    string %f, %f\n",xcm,zcm);
  fprintf(fp,"@    string color 4\n");
  fprintf(fp,"@    string rot 0\n");
  fprintf(fp,"@    string font 0\n");
  fprintf(fp,"@    string just 0\n");
  fprintf(fp,"@    string char size 1.000000\n");
  fprintf(fp,"@    string def \"%d\"\n",num);
}


void agr_cyl(FILE *fp,float xcm,float zcm,float r,int lw,int lc,int fpt,int fc)
{
  fprintf(fp,"@ with ellipse\n");
  fprintf(fp,"@ ellipse on\n"); 
  fprintf(fp,"@ ellipse loctype world\n");
  fprintf(fp,"@ ellipse g0\n");
  fprintf(fp,"@ ellipse %f, %f, %f, %f\n"
	  ,xcm-r,zcm-r,
	  xcm+r,zcm+r);
  fprintf(fp,"@ ellipse linestyle 1\n");
  fprintf(fp,"@ ellipse linewidth %d\n",lw);
  fprintf(fp,"@ ellipse color %d\n",lc);
  fprintf(fp,"@ ellipse fill color %d\n",fc);
  fprintf(fp,"@ ellipse fill pattern %d\n",fpt);
  fprintf(fp,"@ ellipse def\n");
}

void agr_box(FILE *fp,float xcm,float zcm,float MTlen, float MTwidth,
	     float lw,int lc,int fpt,int fc)
{ 
  fprintf(fp,"@ with box\n");
  fprintf(fp,"@ box on\n"); 
  fprintf(fp,"@ box loctype world\n");
  fprintf(fp,"@ box g0\n");
  fprintf(fp,"@ box %f, %f, %f, %f\n"
	  ,xcm-0.5*MTlen,zcm-0.5*MTwidth,
	  xcm+0.5*MTlen,zcm+0.5*MTwidth);
  fprintf(fp,"@ box linestyle 0\n");
  fprintf(fp,"@ box linewidth %f\n",lw);
  fprintf(fp,"@ box color %d\n",lc);
  fprintf(fp,"@ box fill color %d\n",fc);
  fprintf(fp,"@ box fill pattern %d\n",fpt);
  fprintf(fp,"@ box def\n");
}


void agr_arrowline(FILE *fp, int gn, float lw, int lc, 
		   float pendx, float pendz, float mendx, float mendz, 
		   int arrow)
{
  /* arrow=2 -- arrow at end */
  /* arrow=3 -- arrow at both ends */

  int linestyle=1;
  
  int arrowtype=0;
  float arrowlength=1.0;

  fprintf(fp,"@with line\n");
  fprintf(fp,"@    line on\n");
  fprintf(fp,"@    line loctype world\n");
  fprintf(fp,"@    line g%d\n",gn);
  fprintf(fp,"@    line %9.6f,%9.6f,%9.6f,%9.6f\n",
	  pendx,pendz,mendx,mendz);
  fprintf(fp,"@    line linewidth %4.1f\n",lw);
  fprintf(fp,"@    line linestyle %d\n",linestyle);
  fprintf(fp,"@    line color %d\n",lc);
  fprintf(fp,"@    line arrow %d\n",arrow);
  fprintf(fp,"@    line arrow type %d\n",arrowtype);
  fprintf(fp,"@    line arrow length %f\n",arrowlength);
  fprintf(fp,"@    line arrow layout 1.000000, 1.000000\n");
  fprintf(fp,"@line def\n");
}

