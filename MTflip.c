/*###############################################################*/
/* Flip MTs according to a constant flipping rate, e.g. +end out
becomes -end out and vice versa */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include "global_var.h"
#include "ran1.h"
#include "nrutil.h"

void flipMTs()
{
	int iMT,j;
	float tmp;
	
	for (iMT=1;iMT<=MT.number; iMT++)
	{
		if (ran1(&idum)<pFlp)
		{
			MT.direct[iMT].x=MT.direct[iMT].x*(-1);
			tmp=MT.pend[iMT].x;
			MT.pend[iMT].x=MT.mend[iMT].x;	
			MT.mend[iMT].x=tmp;	
		}
	}	
}

