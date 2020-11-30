/* #################################################################### */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "global_var.h"

void motorvec(int ia, int ib, struct point *motorhead, 
		 struct point *motorleg)
{
  /* Generates a vector of motor direction in given 
     overlap region between MT segments ia and ib */

  int iup,ilw;
  float upl,upr,lwl,lwr,xmid;
  float dx=0.2;

  /*if (overlap(ia,ib)==0)
    return;*/


  /* determine which MT segment is the upper one 
     and which is the lower one */

  if (seg.cm[ia].z>seg.cm[ib].z) 
    {iup=ia; ilw=ib;}
  else
    {iup=ib; ilw=ia;}
  
  /* determine the end points of each MT segment */

  upr=seg.cm[iup].x+0.5*seg.length[iup]; /* right end of uppper MT seg */
  upl=seg.cm[iup].x-0.5*seg.length[iup]; /* left end of upper MT seg */
  
  lwr=seg.cm[ilw].x+0.5*seg.length[ilw];
  lwl=seg.cm[ilw].x-0.5*seg.length[ilw];

  /* Distingwish four overlapping configurations */
  /* and determine the mid-point of the overlap region */
 
  if ((upl<=lwr && lwr<=upr) && (lwl<=upl && upl<=lwr))
    /* case 1. */ 
    /*             -----|---------------   */
    /* -----------------|------            */
    /*                 xmid                */
    xmid=lwr-0.5*(lwr-upl);
  
  else if ((lwl<=upr && upr<=lwr) && (upl<=lwl && lwl<=upr))
    /* case 2. */ 
    /* -----------------|---                       */
    /*               ---|---------------------     */
    /*                 xmid                        */
    xmid=upr-0.5*(upr-lwl);
  
  else if ((upl<=lwl && lwl<=upr) && (upl<=lwr && lwr<=upr))
    /* case 3. */ 
    /* -----------------|---------------------     */
    /*             -----|-----                     */
    /*                 xmid                        */
    xmid=lwl+0.5*(lwr-lwl);
  
  else if ((lwl<=upl && upl<=lwr) && (lwl<=upr && upr<=lwr))
    /* case 4. */ 
    /*               ---|---                       */
    /* -----------------|---------------------     */
    /*                 xmid                        */
    xmid=upl+0.5*(upr-upl);

  /* Determine a vector for the motor from the lower MT
     to the upper MT (or vise versa) at the overlap region 
     the direction of this motor is from the head to the leg. 
     The function calculates the head and leg coordinates 
     of the motor. */

 
  
  if (ovlp.type[iup][ilw]==LEGUP)
    {
      if (seg.direct[iup].x<0) dx=-dx;
      motorhead->x=xmid-dx;
      motorhead->y=seg.cm[ilw].y;
      motorhead->z=seg.cm[ilw].z;
      motorleg->x=xmid+dx;
      motorleg->y=seg.cm[iup].y;
      motorleg->z=seg.cm[iup].z;
      
    }
  if (ovlp.type[iup][ilw]==LEGDOWN)
    {
      if (seg.direct[ilw].x<0) dx=-dx;
      motorhead->x=xmid-dx;
      motorhead->y=seg.cm[iup].y;
      motorhead->z=seg.cm[iup].z;
      motorleg->x=xmid+dx;
      motorleg->y=seg.cm[ilw].y;
      motorleg->z=seg.cm[ilw].z;
      
    }

  if (ovlp.type[iup][ilw]==BIPOLAR)
    {
      /* Parallel pair */
      if (seg.direct[ilw].x==seg.direct[iup].x) 
	{
	  if (MT.vel[segindex(iup)].x>MT.vel[segindex(ilw)].x)
	    {
	      motorhead->x=xmid-dx;
	      motorhead->y=seg.cm[ilw].y;
	      motorhead->z=seg.cm[ilw].z;
	      motorleg->x=xmid+dx;
	      motorleg->y=seg.cm[iup].y;
	      motorleg->z=seg.cm[iup].z;
	    } 
	  else
	    {
	      motorhead->x=xmid+dx;
	      motorhead->y=seg.cm[ilw].y;
	      motorhead->z=seg.cm[ilw].z;
	      motorleg->x=xmid-dx;
	      motorleg->y=seg.cm[iup].y;
	      motorleg->z=seg.cm[iup].z;
	    }
	}
      else
	/* Anti-Parallel pair */
	{
	  if (seg.direct[ilw].x>0)
	    {
	      motorhead->x=xmid-dx;
	      motorhead->y=seg.cm[ilw].y;
	      motorhead->z=seg.cm[ilw].z;
	      motorleg->x=xmid+dx;
	      motorleg->y=seg.cm[iup].y;
	      motorleg->z=seg.cm[iup].z;
	    }
	  else
	    {
	      motorhead->x=xmid+dx;
	      motorhead->y=seg.cm[ilw].y;
	      motorhead->z=seg.cm[ilw].z;
	      motorleg->x=xmid-dx;
	      motorleg->y=seg.cm[iup].y;
	      motorleg->z=seg.cm[iup].z;
	    }
	}
    }
  if (ovlp.type[iup][ilw]==TWOMOTORS)
    {
      motorhead->x=xmid-dx;
      motorhead->y=seg.cm[ilw].y;
      motorhead->z=seg.cm[ilw].z;
      motorleg->x=xmid+dx;
      motorleg->y=seg.cm[iup].y;
      motorleg->z=seg.cm[iup].z;
    }
  
  if (ovlp.type[iup][ilw]==BUNDLING)
    {
      motorhead->x=xmid;
      motorhead->y=seg.cm[ilw].y;
      motorhead->z=seg.cm[ilw].z;
      motorleg->x=xmid;
      motorleg->y=seg.cm[iup].y;
      motorleg->z=seg.cm[iup].z;
    }
}
