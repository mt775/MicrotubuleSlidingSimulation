#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "global_var.h"
#define UNIPOLAR (type==LEGUP || type==LEGDOWN)

void prob_update(){

    int ibox, type, iseg, jseg,i;
    float sz, min_ovlp;
    float iseg_l,jseg_l;

    int active_kin[nbox+1],active_dyn[nbox+1], active_bi[nbox+1], notactive[nbox+1];

    float dx=0.00000001+(MT.xmax_real+MT.length[MT.imax_real]/2-MT.xmin_real+MT.length[MT.imin_real]/2)/nbox;
    int active;
    
    for (i=1; i<=nbox; i++) {
        active_kin[i]=0;
        active_dyn[i]=0;
        active_bi[i]=0;
        notactive[i]=0;
        
    }

    for (iseg=1; iseg<=2*MT.number; iseg++) {
        for (jseg=iseg+1; jseg<=2*MT.number; jseg++) {
            sz=soverlap(iseg,jseg);
            type=ovlp.type[iseg][jseg];

            if (sz>EPSI) {
                iseg_l=seg.cm[iseg].x-seg.length[iseg]/2;
                jseg_l=seg.cm[jseg].x-seg.length[jseg]/2;
                min_ovlp=max(iseg_l, jseg_l);
                /*get ibox which contains the left end of the overlap*/
                ibox=(int)floor((min_ovlp-MT.xmin_real+MT.length[MT.imin_real]/2)/dx)+1;
                /*a floating point error may lead to ibox=0,nbox has to be corrected*/
                if (ibox<=0 || ibox>nbox) {
                    if (ibox==0) {
                        ibox=1;
                    }
                    else if (ibox==nbox){
                        ibox=nbox-1;
                    }
                    else{
                        printf("FATAL ERROR in overlapmap!!\n");
                        printf("minovlp=%f\n", min_ovlp);
                        printf("dx=%f\n", dx);
                        printf("box=%i\n", ibox);
                        printf("MT.xmin=%f\n", MT.xmin_real-MT.length[MT.imax_real]/2);
                        printf("MT.xmax=%f\n\n", MT.xmax_real+MT.length[MT.imin_real]/2);
                        printf("MT.length=%f\n", MT.length[MT.imin_real]);
                        printf("isegl=%f\n", iseg_l);
                        printf("isegl=%f\n", jseg_l);
                        printf("iseg=%i\n", iseg);
                        printf("MTcm=%f\n",MT.cm[segindex(iseg)].x);
                        printf("MT left end=%f\n",MT.cm[segindex(iseg)].x-MT.length[segindex(iseg)]/2);
                        printf("checksum=%f\n",(min_ovlp-MT.xmin_real+MT.length[MT.imin_real]/2)/dx);
                        exit(0);
                    }
                }


                
                while (Interval_overlap(min_ovlp, (ibox-1)*dx+MT.xmin_real-MT.length[MT.imin_real]/2, sz, dx )> 0.001) {
                    
                    if (type==ZERO) {
                        notactive[ibox]+=1;
                    }
                    else if (UNIPOLAR && ovlp.motor_direction[iseg][jseg]==-1) {
                        active_kin[ibox]+=1;
                    }
                    else if (UNIPOLAR && ovlp.motor_direction[iseg][jseg]==1) {
                        active_dyn[ibox]+=1;
                    }
                    else if (type==BIPOLAR){
                        active_bi[ibox]+=1;
                        printf("BIPOLAR DETECTED! IS THIS OK?\n");
                        exit(0);
                    }
                    else{
                        printf("FATAL ERROR: unexpected motor type!\n");
                        exit(0);
                    }
                    ibox++;
                }
            }
        }
    }

    
    /*calculate new distribution*/
    for (ibox=1; ibox<=nbox; ibox++) {
        active=active_kin[ibox]+active_dyn[ibox]+active_bi[ibox];
        
        
        
        /*If there are no active overlaps set active to zero but keep 
         Ki and Dy dist the same as we cannot see any changes in them*/
        if (active==0 && notactive[ibox]==0) {

            //ProbAct_dist[ibox]=0.;
            /*ProbBi_dist[ibox]=0;
            ProbKi_dist[ibox]=0;
             
            ProbDy_dist[ibox]=0;*/
           
        }
        else if (active==0){
            //ProbAct_dist[ibox]=0.;
            /*ProbBi_dist[ibox]=0;
            ProbKi_dist[ibox]=0;
            ProbDy_dist[ibox]=0;*/
         
            
        }
        else{
            
            //ProbAct_dist[ibox]=(float)active/(float)(active+notactive[ibox]);
            ProbBi_dist[ibox]=(float)active_bi[ibox]/(float)active;
            ProbKi_dist[ibox]=(float)active_kin[ibox]/(float)active;
            ProbDy_dist[ibox]=(float)active_dyn[ibox]/(float)active;
            if(ProbBi_dist[ibox]>0){
                printf("wild BIPOLAR MOTOR appears out of nowhere! Serious error anticipated! If not intended\n");
                exit(0);
            }
            
        
        }
        
    }


}
