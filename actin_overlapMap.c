/* #################################################################### */
/*If actin flag is set to TRUE this function generates the connections  */
/*of MT close to the boundary to the myosin/actin network               */
/*They can associate with probility actin.attach and deassociate with   */
/*actin.deattach                                                        */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ran1.h"
#include "global_var.h"



void actin_overlapMap(){

    float actin_rad= bundle_rad()+exclude, rand;

    float MT_rad, active_ovlps ,active,notactive;
    int iMT;
    
    active=0.0;
    notactive=0.0;

    for ( iMT=1; iMT<=MT.number; iMT++) {
        /*get distance from 0*/
        MT_rad= sqrtf(MT.cm[iMT].y*MT.cm[iMT].y+MT.cm[iMT].z*MT.cm[iMT].z);
        
        /* Calculate the percentage of already active overlaps */

        /*First check for actin overlaps*/
        /*If there is a actin overlap*/
        if ( (actin_rad-MT_rad)< MinOvlpDist ){
            if (ovlp.actin[iMT]==-1 || ovlp.actin[iMT]==1)
                active++;
            else if (ovlp.actin[iMT]==0 )
                notactive++;
            else{
                printf("ERROR in actin motors!\n");
                printf("iMT%i has not defined overlap\n", iMT);
                exit(0);
            }
        }
    }
    if ( (active+ notactive) >0 )
        active_ovlps=active/(active+notactive);
    else
        active_ovlps=0.;
    

    
    for (iMT=1; iMT<=MT.number; iMT++) {
        /*If not connected check if MT is close to boundary 
         given by MT.zmax and MT.ymax*/
        /*Distribute new actin ovlps according to active fraction*/
        /*get radius of iMT and look if overlap is possible*/
        MT_rad= sqrtf(MT.cm[iMT].y*MT.cm[iMT].y+MT.cm[iMT].z*MT.cm[iMT].z);
        
        /*if ((actin_rad-MT_rad)< MinOvlpDist) {

            printf("%i\n", ovlp.actin[iMT]);
            printf("MT_rad=%f sctrad=%f\n", MT_rad, actin_rad);
            printf("%f\n",active_ovlps);
        }*/

        rand=ran1(&idum);
        
        if ( ovlp.actin[iMT]==0 && (actin_rad-MT_rad)< MinOvlpDist
            && rand< (actin.attach-active_ovlps) ){

            if (ran1(&idum)>ProbBipolar) {
                if ( ran1(&idum)< ProbKinesin) {
                    //printf("%f %f CONNECT! KINESIN MT%i\n", rand, actin.attach-active_ovlps, iMT);
                    ovlp.actin[iMT]=-1;
                    ovlp.actin_it0[iMT]=iter;
                }
                else{
                    //printf("%f %f CONNECT! DYNEIN MT%i\n",rand,actin.attach-active_ovlps , iMT);
                    ovlp.actin[iMT]=1;
                    ovlp.actin_it0[iMT]=iter;
                }
            }
            else{
                printf("ERROR: actin network together with bipolar motors not implemented!\n");
                exit(0);
            }
        }
        /*There is a certain runtime after which a filament detaches from the actin
         */
        if (ovlp.actin[iMT]==1 || ovlp.actin[iMT]==-1 ) {
            
            if (ran1(&idum)< actin.deattach) {
                ovlp.actin[iMT]=0;
            }
        }
      
        
    }
    
}
