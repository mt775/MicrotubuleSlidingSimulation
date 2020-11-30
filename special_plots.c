/* #################################################################### */
/* special plot file be carful changing something here                  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global_var.h"


int special_plots(){
    
    /*The following reinterpret the av loop to raising the
     stall force one av iteration to the other and plot the stall force of the bundle*/
    
    if (forcelength){
        if (spring_r<0.000001) {
            printf("ERROR: force length (stallmode) mode only makes sense for spring!\n");
            exit(0);
        }
        
        printf("STEADY LENGTH CALC MODE ON\n");
        write_force_fstall();
        fstall_k+=0.1;
        fstall_d+=0.1;
        return(0);
    }
    
    else if (0){
        if (spring_r>0.000001) {
            printf("ERROR: force length mode only makes sense for spring!\n");
            exit(0);
        }
        
        printf("STEADY LENGTH CALC MODE ON\n");
        wav_stalllength_extf();
        ext_fr-=0.1;
        return(0);
    }
    
    /*The following 3 lines reinterpret the av loop to raising the
     polarity ratio from one av iteration to the other*/
     else if (forcepol) {
        printf("POLARITY AGAINST FORCE PLOT MODE ON\n");
        write_force_polarity();
        
        PolarityRatio0+= 1/(float)nav;
        PolarityRatioNew+= 1/(float)nav;
        if (PolarityRatio0!=PolarityRatioNew) {
            printf("Warning: The FOrce(Polarity) plot only functions if PolarityRatio0==PolarityRatioNew!!\n");
            exit(0);
        }
         return(0);
     }
    /*The following reinterprets the av loop to raising the
     motor fraction from one av iteration to the other*/
     else if (0) {
         
         printf("Motor Fraction AGAINST FORCE PLOT MODE ON\n");
         write_force_motorfraction();
         wav_percprob_bi();
         ProbBipolar+= 1/(float)nav;
         ProbKinesin-= 1/(float)nav;
         return(0);
     }
    
    
    
    
            /*The following raises the spring constant at each av step*/
     else if (0) {
        printf("spring const against FORCE PLOT MODE ON\n");
        write_force_spring_r();
        
        spring_r+= 1;

        return(0);
     }
     
    
    /*If the following switch is on the averaging iterations make place 
     for iterations raising the chi step by step to look for a percolation transition*/
     else if (0){
         
         /*check if chi is smallest in first interation*/
         if (iav==1 && ProbActive>0.0) {
             printf("ERROR! chi>0");
             exit(0);
         }
         printf("Percolation Probability (+stationary force) vs CHI on!\n");
         printf("No av data will be produced!\n");
         wav_percprob();
         wav_stallforce();
         ProbActive+=1./(float)nav;
         
        return(0);
     }
    
    
    return(1);
    
}
