/*-----------------------------------------------*/
/* find the jMT whose index at iteration it, is i */

#include <stdio.h>
#include <math.h>
#include "global_var.h"

int getjMT(int i,int it)
{
  int jMT,j,nMT;

  jMT=0;
  nMT=track.nMT[it]; //number of MT in iteration it
  do
    {
      jMT++;
      j=track.index[it][jMT]; // track the index of MT number jMT in iteration it
      if (jMT>nMT)
	{
	  jMT=0; // gives zero if microtubule not there
	  break;
	}
    }
  while (j != i); // if the index j is equal to given i the current indexing gets returned
  
  /* printf("i=%d j=%d index[%d][%d]=%d",i,j,it2,jMT,track.index[it2][jMT]); */
/*   getchar(); */

  return jMT;
}
