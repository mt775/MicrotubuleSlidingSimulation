/* ######################################################################### */
/* grow -- makes MTs grow at their plus end and may initiate instability     */

#include <stdio.h>
#include "global_var.h"

void grow(){
    /*Should add those parameters to data*/
    float instab_p=0.001;
    float rescue_p=0.5;
    float rand;
    int iseg,iMT;
    
    for (iMT=1; iMT<=MT.number; iMT++) {
        /*If MT is currently growing*/
        if (MT.grow_direct[iMT]==1) {
            MT.cm[iMT].x+= 0.5*velpol*dt;
            if ( ran1(&idum) < instab_p ){
                MT.grow_direct[iMT]=-1;
            }
        }
        /*If MT is shrinking*/
        else if (MT.grow_direct[iMT]==-1) {
            MT.cm[iMT].x+= 0.5*veldepol*dt;
            if ( ran1(&idum) < rescue_p ) {
                MT.grow_direct[iMT]=1;
            }
        }
        /*IF MT is idle*/
        else if (MT.grow_direct[iMT]==0) {
            if ( ran1(&idum) < instab_p ) {
                MT.grow_direct[iMT]=-1;
            }
            else{
                MT.grow_direct[iMT]=1;
            }
        }
        else{
            printf("Fatal Error in grow.c!\n MT growth not specified!\n")
        }
    }
}
