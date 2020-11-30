#include <string.h>
#include <stdio.h>
#include <math.h>
#include "global_var.h"

/*#############################################################*/
//THIS FUNCTION DIVIDES THE ZONE INTO NBOX REGIONS AND OBTAINS AV MOTOR VELOCITIES PER REGION
/*#############################################################*/
void Box_av_vel (float boxlength, float *av_vel){
    
    struct BOXsegmentdata {
        float length_i,length_j, leftend_i, leftend_j;
    } boxseg;
    int high,iMT,jMT;
    float sz=0, tot_seg_l=0,tot_seg_r=0; // tot_seg is the total length of all segments in the ibox
                                 // free_seg is the length of segment surface unoccupied by motors
    
    for (int ibox=1; ibox <= nbox; ibox++) {
        
        // SEARCH THROUGH ALL FILAMENT SECTIONS IN BOX
        /*#############################################################*/
        for (int iseg=1; iseg <= 2*MT.number; iseg++) {
            // only carry on in this loop if the segment is existing
            if (seg.length[iseg]<0.0001) continue;

            iMT=segindex(iseg);
            
            if (seg.direct[iseg].x < 0)  tot_seg_l+=boxseg.length_i; // update the total segment length for minus end to left in the box
            if (seg.direct[iseg].x > 0)  tot_seg_r+=boxseg.length_i; // update the total segment length for minus end to right in the box
            
        
            //LOOK FOR ACTIVE OVERLAPS
            /*#############################################################*/
            for (int jseg=1; jseg <= 2*MT.number; jseg++) {
                
                if (seg.length[jseg]<0.0001) continue;
                

                jMT=segindex(jseg);
                
                sz=sbox_overlap(iseg,jseg,ibox,boxlength);
                // IF THERE IS AN OVERLAP BETWEEN 2 SEGMENTS IN THE BOX CHECK FOR MOTORS THERE
                if (  sz > 0 && ovlp.type[iseg][jseg]!=ZERO ) {
                    
                    //IF UNIPOLAR MOTOR CONNECTION
                    if (ovlp.type[iseg][jseg] == LEGUP || ovlp.type[iseg][jseg] == LEGDOWN) {
                        
                        // GET high filament
                        /*#############################################################*/
                        if (MT.cm[iMT].z>MT.cm[jMT].z
                            || ((MT.cm[iMT].z==MT.cm[jMT].z) &&
                                (MT.cm[iMT].y>MT.cm[jMT].y))
                            /* or, arrbitrary choice (iMT>jMT) if they overlap */
                            || ((MT.cm[iMT].z==MT.cm[jMT].z) &&
                                (MT.cm[iMT].y==MT.cm[jMT].y) && iseg>jseg))
                            high=iMT;
                        else
                            high=jMT;
                        /*#############################################################*/
                        /////////////////////////////////
                        av_vel[ibox]+= sz*MT.vel[high].x;
                        /////////////////////////////////
                    }
                    
                    //IF BIPOLAR MOTORS
                    else if (ovlp.type[iseg][jseg] == BIPOLAR ){
                        if (ovlp.motor_direction[iseg][jseg] == 1 )
                            av_vel[ibox]+= sz*0.5*(MT.vel[iMT].x+MT.vel[jMT].x+vel0_d*(seg.direct[iseg].x+seg.direct[jseg].x));
                        else if (ovlp.motor_direction[iseg][jseg] == -1 )
                            av_vel[ibox]+= sz*0.5*(MT.vel[iMT].x+MT.vel[jMT].x+vel0_k*(seg.direct[iseg].x+seg.direct[jseg].x));
                        
                    }
                }
                
            }
            /*#############################################################*/
        }
        /*#############################################################*/
        
        //LET MOTORS WALK ON THE FILAMENT WITH THEIR FREE VELOCITY
        /*#############################################################*/
        
        
        
        /*#############################################################*/
    }
    

}
/*#############################################################*/

void update_motor_dist () {
    
    // Divide x axis into boxes of length boxlength 
    float boxlength= (MT.xmax+0.5*MT.length[MT.imax]-MT.xmin+0.5*MT.length[MT.imin])/(float)nbox;
    
    float *av_vel = vector (1,nbox);
    empty(av_vel,nbox);
    
    /*#############################################################*/
    //Obtain average velocities in every box
    Box_av_vel(boxlength, av_vel);
    /*#############################################################*/
    
    
}
