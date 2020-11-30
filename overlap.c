/* ######################################################################## */
/* Overlap finds the overlap region between two MTs                         */
/* We assume the two MTs are in the x-z plane.                              */

#include <stdio.h>
#include <math.h>
#include "global_var.h"
/* #include "func_declaration.h" */

double overlap(int ia, int ib)
/* Returns the total overlap between two MTs */
{
  double ov,ov1,ov2,ov3,ov4;
  int iseg1,iseg2,jseg1,jseg2;

  /* First check the separation between the MTs     */
  /* Then check that both MTs are in the domain     */
  /* If they are in the doamin and are close enough, 
     calculate the overlapping length */

  if (ia==ib) return 0.0;
  iseg1=2*ia-1;
  iseg2=2*ia;
  
  jseg1=2*ib-1;
  jseg2=2*ib;

  
  ov1=soverlap(iseg1,jseg1);
  ov2=soverlap(iseg1,jseg2);
  ov3=soverlap(iseg2,jseg1);
  ov4=soverlap(iseg2,jseg2);

  if (ovlp.type[iseg1][jseg1]==ZERO) ov1=0.0;
  if (ovlp.type[iseg1][jseg2]==ZERO) ov2=0.0;
  if (ovlp.type[iseg2][jseg1]==ZERO) ov3=0.0;
  if (ovlp.type[iseg2][jseg2]==ZERO) ov4=0.0;

  ov=ov1+ov2+ov3+ov4;
  return ov;
}


double xoverlap(int ia, int ib)
/* Returns the total overlap between two MTs irrespective */
/* of the MTs position in the y-z plane.                  */
{
  double ov,ov1,ov2,ov3,ov4;
  int iseg1,iseg2,jseg1,jseg2;

  /* First check the separation between the MTs     */
  /* Then check that both MTs are in the domain     */
  /* If they are in the doamin and are close enough, 
     calculate the overlapping length */

  if (ia==ib) return 0.0;

  iseg1=2*ia-1;
  iseg2=2*ia;
  
  jseg1=2*ib-1;
  jseg2=2*ib;

  ov1=sxoverlap(seg.cm[iseg1].x,seg.cm[jseg1].x,
		seg.length[iseg1],seg.length[jseg1]);
  ov2=sxoverlap(seg.cm[iseg1].x,seg.cm[jseg2].x,
		seg.length[iseg1],seg.length[jseg2]);
  ov3=sxoverlap(seg.cm[iseg2].x,seg.cm[jseg1].x,
		seg.length[iseg2],seg.length[jseg1]);
  ov4=sxoverlap(seg.cm[iseg2].x,seg.cm[jseg2].x,
		seg.length[iseg2],seg.length[jseg2]);

  ov=ov1+ov2+ov3+ov4;
  return ov;
}


/* Overlap between MT segments */
double soverlap(int ia, int ib)
{
  double ov,dist2,ar,al,br,bl,mxr,minl;
  int iMT,jMT;

  /* Check the separation between the MT segments  */
  /* Then check that both MT segments are in the domain  */
  /* If they are in the doamin and are close enough, but not too close, 
     calculate the overlapping length */
  /* Interactions between two segments of the same MT are not included.*/

  if (ia==ib) return 0.0;

  /*if (segindex(ia)==segindex(ib)) return 0.0;*/ /* segindex(ia)=iMT */
  if (seg.length[ia]==0.0 || seg.length[ib]==0.0) return 0.0;

  ar=seg.cm[ia].x+0.5*seg.length[ia]; /* right end of MT seg a */
  al=seg.cm[ia].x-0.5*seg.length[ia]; /* left end of MT seg a */
  
  br=seg.cm[ib].x+0.5*seg.length[ib];
  bl=seg.cm[ib].x-0.5*seg.length[ib];

  dist2=(seg.cm[ia].z-seg.cm[ib].z)*(seg.cm[ia].z-seg.cm[ib].z)+
    (seg.cm[ia].y-seg.cm[ib].y)*(seg.cm[ia].y-seg.cm[ib].y);
    
    //printf("%f\t%f\n",dist2,MinOvlpDist);
  
  if ( sqrtf(dist2)>MinOvlpDist || ((sqrtf(dist2)<2*MT.radi) && 
				    (segindex(ia)!=segindex(ib)))
       || (bound.type[0]==ABSORB && (ar<lbound0 || br<lbound0))
       || (bound.type[1]==ABSORB && (al>rbound0 || bl>rbound0))
       || (bound.type[1]==POPUP && (al>rbound0 || bl>rbound0))
       || (bound.type[0]==POPUP && (ar<lbound0 || br<lbound0))
       ) ov=0.0;
  else 
    {
      mxr=maxf(ar,br);        
      minl=minf(al,bl);       	  
      ov=seg.length[ia]+seg.length[ib]-(mxr-minl);
    }

  if (ov>EPSI)
    return (double)ov;
  else
    return (double)0.0;
}




double sxoverlap(double xa, double xb, double lena, double lenb)
{
  /* This function returns the overlap between any two MT segments 
     along the x coordinate, irrespective of the z and y coordinates.
     It is similar to soverlap(), only the latter returns a zero number 
     if the two MTs are too far in the yz plane */

  /* xa, xb -- are the centers of mass of the segments and,*/
  /* lena, lenb -- their length */
  /* the segments could also be the whole MTs */

  double ov,ar,al,br,bl,mxr,minl;

  ar=xa+0.5*lena; /* right end of MT a */
  al=xa-0.5*lena; /* left end of MT a */
  
  br=xb+0.5*lenb;
  bl=xb-0.5*lenb;
    
  mxr=maxf(ar,br);        
  minl=minf(al,bl);       	  
  ov=lena+lenb-(mxr-minl);
  
  if (ov>EPSI)
    return (ov);
  else
    return (0.0);
}
