/* #################################################################### */
/* wdata  --  Generate output files                                     */

#include <stdio.h>
#include <math.h>
#include "global_var.h"

void wdata()
{
  /*int nbox;*/

  /* MT Map */ 
  if ((iter+wfreq-1)%wfreq==0 /*|| (iter>2417 && iter<2474)*/ && iav==1) wmap_xz(iter);
  if ((iter+wfreq-1)%wfreq==0 && iav==1) wmap_yz(iter);
  if ((iter+wfreq-1)%wfreq==0 /*|| (iter>2417 && iter<2674)*/ && iav==1) wvelmap(iter);

  /* Density profile */
  /*nbox=15;
    if ((iter+wfreq-1)%wfreq==0) wdensity(iter, nbox);*/
}
