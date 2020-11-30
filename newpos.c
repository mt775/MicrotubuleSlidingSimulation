/* ######################################################################### */
/* newpos -- moves the MTs to a new position in accordence with MT.vel       */

#include <stdio.h>
#include <stdlib.h>
#include "global_var.h"
#include "math.h"

void newpos()
{
  int iMT,flag=0;
  float dx;
  	
  /* store new configuration */
  for (iMT=1; iMT<=MT.number; iMT++)
    {
      dx=dt*MT.vel[iMT].x;
      /*if (fabs(dx)>MT.length[iMT])
      {
      	 if (flag==0)
         {
	     //printf("iter=%i MT moves too fast (dx>length)! readjust to dx=MT.length!\n",iter);
             flag=1;
	 }
	 //dx=dx/fabs(dx)*0.1*MT.length[iMT];
      }	*/
      MT.cm[iMT].x += dx;
      MT.mend[iMT].x += dx;
      MT.pend[iMT].x += dx;
      MT.direct[iMT].x=(MT.mend[iMT].x-MT.pend[iMT].x)/MT.length[iMT];
      if (iter>neqsteps) MT.rcm[iMT].x += dx; 
      /* real x position of iMT; not affected by b.c. */
        /*if (bound.type[0]==WALL && MT.cm[iMT].x<lbound0 ) {
            printf("MT moved beyond WALL boundary!\n");
            exit(0);
        }*/
    }
    

  /* apply boundary conditions */
  boundary();
}
